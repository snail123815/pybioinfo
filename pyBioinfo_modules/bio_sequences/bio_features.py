from concurrent.futures import ThreadPoolExecutor, as_completed
from copy import deepcopy
from typing import Literal

from Bio.Seq import Seq
from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    ExactPosition,
    FeatureLocation,
    SeqFeature,
)
from Bio.SeqRecord import SeqRecord


def truncate_feat_translation(
    feat: SeqFeature,
    side: Literal["left", "right", "both_sides"],
    codon_table=11,
    on_seq: Seq | None = None,
    inplace: bool = True,
) -> None | SeqFeature:
    """Truncate the translation of a truncated CDS feature.
    The feature must
        1. Already been truncated to the correct location (length).
        2. Have a "translation" qualifier.
    Only the "translation" qualifier will be modified.
    The feature length will be used to generate new translation.

    Args:
        feat (SeqFeature): Truncated feature.
        side (Literal["left", "right", "both_sides"]): Side to truncate the translation.
        on_seq (Seq): Sequence where the feature lies, only necessary for "both_sides" truncation.

    Returns:
        None: The feature will be modified in place.
        | SeqFeature: If inplace is False, a new feature will be returned.
    """
    if not "translation" in feat.qualifiers:
        raise ValueError("Feature must have a 'translation' qualifier")
    # Reverse translation if the feature is on the -1 strand
    # so that it follows the actual DNA sequence
    if not inplace:
        feat = deepcopy(feat)
    translation = feat.qualifiers["translation"][0]
    if feat.location.strand == -1:
        translation = translation[::-1]

    max_translation_len = (
        len(feat) // 3
    )  # only used for left and right truncation

    if side == "left":
        translation = translation[len(translation) - max_translation_len :]
    elif side == "right":
        translation = translation[:max_translation_len]
    elif side == "both_sides":
        # Now I cannot know where the codon starts, so have to guess
        assert (
            on_seq is not None
        ), "Sequence (on_seq parameter) must be provided for 'both_sides' truncation"
        possible_ts = []
        for i in [0, 1, 2]:
            s = feat.extract(on_seq)[i:]
            t = s.translate(table=codon_table, to_stop=False, stop_symbol="")
            if len(t) >= max_translation_len - 2:
                if str(t) in feat.qualifiers["translation"][0]:
                    possible_ts.append((i, t))
        if len(possible_ts) == 0:
            translation = ""
            codon_start = 0
        else:
            possible_ts = sorted(possible_ts, key=lambda x: len(x[1]))
            possible_t = possible_ts[-1]
            codon_start = possible_t[0]
            translation = possible_t[1]

        if "note" not in feat.qualifiers:
            feat.qualifiers["note"] = []
        feat.qualifiers["note"].append(
            "Truncated translation was ambiguous, "
            "the longest possible translation was chosen. "
            f"Codon start position: {codon_start} (0-based) "
            f"from strand {feat.location.strand}"
        )
    else:
        raise ValueError(
            f"Side must be 'left', 'right' or 'both_sides' not {side}"
        )

    # Reverse translation back if the feature is on the -1 strand
    if feat.location.strand == -1 and side != "both_sides":
        translation = translation[::-1]

    feat.qualifiers["translation"] = [translation]

    if not inplace:
        return feat


def find_truncated_features(
    source_seq: SeqRecord,
    location: tuple[int, int] | FeatureLocation,
    expand: int = 20000,
    include_inner_feats: bool = False,
) -> list[SeqFeature]:
    """Get features from a sequence region with proper handling of truncated
    features. Scanning the region with optional expansion, for features that
    might be truncated by splicing.
    Used to append features back to the sequence after slicing.
    Can find features that span the location bundaries, or are completely
    within the region.

    Args:
        source_seq (SeqRecord): Source sequence record.
        location (tuple[int, int] | FeatureLocation): Start and end positions
            in genome coordinates (0-based, inclusive start, exclusive end) or
            a FeatureLocation object.
        expand (int, optional): Number of bases to expand the region by.
            Defaults to 20000 to include larger features produced by antismash.
        include_inner_feats (bool, optional): Whether to include features that
            are completely within the region. Defaults to False.

    Returns:
        list[SeqFeature]: List of features in the region with proper truncation handling
    """
    # Calculate region boundaries with expansion
    if isinstance(location, FeatureLocation):
        start = int(location.start)
        end = int(location.end)
    else:
        start, end = location
    expand_start = max(0, start - expand)
    expand_end = min(
        len(source_seq), end + expand
    )  # +1 to include end position

    spanFeats = []
    for feat in source_seq[expand_start:expand_end].features:
        all_in = (
            start <= feat.location.start <= end
            and start <= feat.location.end <= end
        )
        if all_in and not include_inner_feats:
            # Feature will be included by splicing
            continue
        left_in = (
            start <= feat.location.start <= end and feat.location.end > end
        )
        right_in = (
            feat.location.start <= start and start <= feat.location.end <= end
        )
        span = feat.location.start <= start and feat.location.end >= end
        if sum([all_in, left_in, right_in, span]) == 0:
            continue
        if sum([all_in, left_in, right_in, span]) > 1:
            assert feat.location.start == feat.location.end, (
                "Feature should be either in one of the categories: "
                "all_in, left_in, right_in, span"
                "or have a length of 1 (start == end)"
            )
            spanFeats.append(deepcopy(feat))

        # Create new feature with copied qualifiers
        newFeat = deepcopy(feat)

        # Feature completely within region
        if all_in:
            newFeat.location = FeatureLocation(
                feat.location.start - start,
                feat.location.end - start,
                feat.location.strand,
            )
        # Feature left side in region, truncated on its right
        elif left_in:
            newFeat.location = FeatureLocation(
                feat.location.start - start,
                AfterPosition(end - start),
                feat.location.strand,
            )
            if feat.type == "CDS":
                truncate_feat_translation(newFeat, side="right")
            newFeat.qualifiers["truncated"] = ["right"]
        # Feature right side in region, truncated on its left
        elif right_in:
            newFeat.location = FeatureLocation(
                BeforePosition(0),
                feat.location.end - start,
                feat.location.strand,
            )
            if feat.type == "CDS":
                truncate_feat_translation(newFeat, side="left")
            newFeat.qualifiers["truncated"] = ["left"]
        # Feature spans the region
        elif span:
            newFeat.location = FeatureLocation(
                BeforePosition(0),
                AfterPosition(end - start),
                feat.location.strand,
            )
            if feat.type == "CDS":
                truncate_feat_translation(
                    newFeat, side="both_sides", on_seq=source_seq[start:end].seq
                )
            newFeat.qualifiers["truncated"] = ["both_sides"]
        else:
            raise ValueError(
                "Feature should be either in one of the categories: "
                "all_in, left_in, right_in, span"
                "or have a length of 1"
            )
        spanFeats.append(newFeat)

    # Sort features by start position
    return sorted(spanFeats, key=lambda f: f.location.start)


def seqFeature_to_tuple(seqFeature) -> tuple:
    """Convert a SeqFeature object to a hashable tuple.
    Note only limited feature info is preserved."""
    qualifier_keys = tuple(sorted(seqFeature.qualifiers.keys()))
    qualifier_values = tuple(
        tuple(seqFeature.qualifiers[key]) for key in qualifier_keys
    )
    return (
        seqFeature.type,
        str(seqFeature.location),
        qualifier_keys,
        qualifier_values,
    )


def slice_sequence_keep_truncated_features(
    source_seqrec: SeqRecord,
    location: FeatureLocation | tuple[int, int, int | None] | tuple[int, int],
    id: str | None = None,
):
    """Slice a SeqRecord with a FeatureLocation object. The sliced sequence
    will have the same features as the source sequence, but they will be
    adjusted to the new sequence coordinates.
    Truncated features will have a "truncated" qualifier added to them,
    indicating which side of this feature is truncated ('left', 'right').

    Args:
        sourceSeq (SeqRecord): Source sequence to slice.
        location (FeatureLocation | Tuple): Location of the slice.
        id (str, optional): ID of the new sequence. Defaults to the origional.

    Returns:
        SeqRecord: Sliced sequence.
                   Attributes 'id' and 'description' will be set to
                   "{sourceSeq.id}_{start}-{end}_rc" if id is None.
    """
    rc = False
    if isinstance(location, FeatureLocation):
        start = location.start
        end = location.end
        if location.strand == -1:
            rc = True
    else:
        start, end = location[:2]
        if len(location) == 3:
            if location[2] == -1:
                rc = True
            elif location[2] != 1 and location[2] is not None:
                raise ValueError(
                    "Strand must be 1, -1 or None\n"
                    f"Current location is {location}"
                )
        # Validate indices
        if start < 0 or end < 0:
            raise ValueError("Negative indices not allowed")
        if start >= end:
            raise ValueError("Start must be less/equal than end")
        if end > len(source_seqrec) + 1:
            raise ValueError("End position exceeds sequence length")
    sliced = source_seqrec[start:end]
    # Expand features if the cut location is inside features
    descrip = []
    # Add features, add gene names to description
    sliced.features.extend(find_truncated_features(source_seqrec, (start, end)))
    if len(sliced.features) > 0:
        for feat in sliced.features:
            if feat.type == "gene":
                try:
                    descrip.append(feat.qualifiers["locus_tag"][0])
                except:
                    pass
    if rc:
        sliced = reverse_complement_SeqRecord_with_features(sliced)
    if id:
        sliced.id = id
    else:
        sliced.id = f"{source_seqrec.id}_{start}-{end}{'_rc' if rc else ''}"
    if descrip:
        descrip.insert(0, sliced.id)
    else:
        descrip = [sliced.id]
    sliced.description = "-".join(descrip).replace(" ", "_")
    return sliced


def slice_seq_concurrent(
    source_seqs: list[SeqRecord], peak_info: dict[str, list[str, list[int]]]
) -> list[SeqRecord]:
    """
    peak_info data structure:
    peak_info -> {peak_id: [chr, [start, end]]}
    """
    extracted_seqs = []
    sources = {}
    for seq in source_seqs:
        sources[seq.id] = seq
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = []
        for peak, info in peak_info.items():
            chr = info[0]
            loc = info[1]
            if chr not in sources:
                raise Exception(f"Chromosome {chr} not found in genome file")
            source_seq = sources[chr]
            futures.append(
                executor.submit(
                    slice_sequence_keep_truncated_features,
                    source_seq,
                    loc,
                    peak,
                )
            )
        for future in as_completed(futures):
            sliceSeq = future.result()
            extracted_seqs.append(sliceSeq)
    return extracted_seqs


def add_seq_to_SeqRecord_as_feature(
    target_seqrec: SeqRecord,
    to_add_seq: Seq,
    feature_type: str = "feature",
    qualifiers: dict = {},
) -> SeqRecord:
    """
    Add sequence to the end of the sequence record. Automatically finds the
    correct strand and location for the target sequence.
    If a "primer_bind" sequence is provided, it will check if there is
    5' overhang, and avoid that in the match. Return a new SeqRecord.

    "overhang" should be a string qualifier in the qualifiers dictionary, can
    a number or sequences.

    Args:
        seq_record (SeqRecord): Sequence record to add the feature to.
        target_seq (Seq): Sequence to add as a feature.
        feature_type (str, optional): Type of feature to add. Defaults to "feature".
        qualifiers (dict, optional): Qualifiers for the feature. Defaults to {}.

    Returns:
        SeqRecord: Sequence record with the added feature.
                   For seq with overhanges, if int was used as "overhang"
                   qualifier, it will be converted to a string.
    """

    if (
        feature_type == "primer_bind" and "overhang" in qualifiers
    ):  # 5' overhang
        try:
            overhang_len = int(qualifiers["overhang"][0])
            overhang_seq = to_add_seq[:overhang_len]
            qualifiers["overhang"] = [str(overhang_seq)]
        except ValueError:
            overhang_len = len(qualifiers["overhang"][0])
            overhang_seq = Seq(qualifiers["overhang"][0])
            assert to_add_seq.lower().startswith(
                str(overhang_seq.lower())
            ), "Overhang not found in target sequence"
        match_seq = Seq(to_add_seq[overhang_len:].lower())
    else:
        match_seq = Seq(to_add_seq.lower())

    seq_seq = target_seqrec.seq.lower()
    seq_len = len(seq_seq)
    match_len = len(match_seq)
    on_strand = 1

    if match_seq in seq_seq:
        pass
    elif match_seq.reverse_complement() in seq_seq:
        on_strand = -1
        match_seq = match_seq.reverse_complement()
    else:
        raise ValueError("Oligo not found in sequence")

    for i in range(seq_len - match_len + 1):
        if seq_seq[i : i + match_len] == match_seq:
            location = FeatureLocation(i, i + match_len, strand=on_strand)
            feature = SeqFeature(
                location=location,
                type=feature_type,
                qualifiers=qualifiers,
            )
            break
    new_seq = target_seqrec[:]
    new_seq.features.append(feature)
    return new_seq


def reverse_complement_SeqRecord_with_features(record: SeqRecord) -> SeqRecord:
    """Reverse complement a SeqRecord, process its features to fit in
    the new SeqRecord."""
    new_seq = record.seq.reverse_complement()
    new_features = []
    seq_length = len(record)

    for f in record.features:
        new_start = reverse_complement_position(seq_length, f.location.end)
        new_end = reverse_complement_position(seq_length, f.location.start)
        new_strand = (
            None
            if f.location.strand is None
            else (-f.location.strand if f.location.strand else 0)
        )
        new_location = FeatureLocation(new_start, new_end, strand=new_strand)
        f.location = new_location
        new_features.append(f)

    new_record = SeqRecord(
        new_seq,
        id=record.id,
        name=record.name,
        features=new_features,
        description=record.description,
        dbxrefs=record.dbxrefs,
        annotations=record.annotations,
        letter_annotations=record.letter_annotations,
    )

    return new_record


def reverse_complement_position(
    seq_length: int, pos: int | BeforePosition | AfterPosition | ExactPosition
) -> int | BeforePosition | AfterPosition | ExactPosition:
    """Reverse complement a single position on a sequence of length
    seq_length. Takes care of fuzzy positions.
    """
    new_pos = seq_length - pos
    if isinstance(pos, BeforePosition):
        return AfterPosition(new_pos)
    elif isinstance(pos, AfterPosition):
        return BeforePosition(new_pos)
    elif isinstance(pos, ExactPosition):
        return ExactPosition(new_pos)
    else:
        return new_pos


def reverse_complement_location(
    seq_length: int,
    location: FeatureLocation | tuple[int, int, int | None] | tuple[int, int],
) -> FeatureLocation:
    """Reverse complement a location (start, end, [strand]) on a sequence
    of length seq_length. Takes care of fuzzy positions.
    ALWAYS returns a FeatureLocation object!

    Args:
        seq_length (int): Length of the sequence.
        location (FeatureLocation | tuple[int, int]):
                          Location to reverse complement. Tuple can contain
                          a third element for the strand [-1,1,None].

    Returns:
        FeatureLocation: Reversed location.
    """
    if isinstance(location, tuple):
        if len(location) == 2:
            location = FeatureLocation(location[0], location[1])
        elif len(location) == 3:
            assert location[2] in [None, 1, -1], "Strand must be 1, -1 or None"
            location = FeatureLocation(location[0], location[1], location[2])
        else:
            raise ValueError("Location must be a tuple of length 2 or 3")
    return FeatureLocation(
        reverse_complement_position(seq_length, location.end),
        reverse_complement_position(seq_length, location.start),
        strand=None if location.strand is None else -location.strand,
    )
