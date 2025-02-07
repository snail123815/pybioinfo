from math import ceil
from copy import deepcopy

from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    ExactPosition,
    FeatureLocation,
    SeqFeature,
)
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def getSpanFetures(
    source_seq: SeqRecord,
    location: tuple[int, int] | FeatureLocation,
    expand: int = 20000,
    include_inner_feats: bool = False,
):
    """Get features from a sequence region with proper handling of truncated
    features. Scanning the region with optional expansion, for features that
    might be truncated by splicing.
    0-based indexing is used for start and end positions.

    Args:
        source_seq (SeqRecord): Source sequence record.
        location (tuple[int, int] | FeatureLocation): Start and end positions in genome coordinates (0-based, inclusive start, exclusive end) or a FeatureLocation object.
        expand (int, optional): Number of bases to expand the region by. Defaults to 20000 to include larger features produced by antismash.
        include_inner_feats (bool, optional): Whether to include features that are completely within the region. Defaults to False.

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
            newFeat.qualifiers["truncated"] = ["right"]
        # Feature right side in region, truncated on its left
        elif right_in:
            newFeat.location = FeatureLocation(
                BeforePosition(0),
                feat.location.end - start,
                feat.location.strand,
            )
            newFeat.qualifiers["truncated"] = ["left"]
        # Feature spans the region
        elif span:
            newFeat.location = FeatureLocation(
                BeforePosition(0),
                AfterPosition(end - start),
                feat.location.strand,
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


# Convert SeqFeature to hashable tuple
def seqFeature_to_tuple(seqFeature):
    qualifier_keys = tuple(sorted(seqFeature.qualifiers.keys()))
    qualifier_values = tuple(
        tuple(seqFeature.qualifiers[key]) for key in qualifier_keys
    )
    return (
        seqFeature.type,
        seqFeature.location.start,
        seqFeature.location.end,
        seqFeature.location.strand,
        qualifier_keys,
        qualifier_values,
    )


def slice_sequence(
    sourceSeq: SeqRecord,
    location: FeatureLocation,
    id: str | None = None,
    with_features=False,
):
    rc = False
    if isinstance(location, FeatureLocation):
        start = location.start
        end = location.end
        if location.strand == -1:
            rc = True
    else:
        assert len(location) == 2
        start, end = location
        # Validate indices
        if start < 0 or end < 0:
            raise ValueError("Negative indices not allowed")
        if start >= end:
            raise ValueError("Start must be less than end")
        if end > len(sourceSeq):
            raise ValueError("End position exceeds sequence length")
    sliced = sourceSeq[start:end]
    # Expand features if the cut location is inside features
    descrip = []
    if with_features:
        # Add features, add gene names to description
        sliced.features.extend(getSpanFetures(sourceSeq, start, end))
        if len(sliced.features) > 0:
            for feat in sliced.features:
                if feat.type == "gene":
                    try:
                        descrip.append(feat.qualifiers["locus_tag"][0])
                    except:
                        pass
    else:
        sliced = SeqRecord(
            sliced.seq,
            name=sliced.name,
            dbxrefs=sliced.dbxrefs,
            annotations=sliced.annotations,
        )
    if rc:
        sliced = reverse_complement_with_features(sliced)
    sliced.id = f"{sourceSeq.id}_{start}-{end}" if id is None else id
    if descrip:
        sliced.description = "-".join(descrip).replace(" ", "_")
    else:
        sliced.description = f"{sourceSeq.id}_{start}-{end}"
    return sliced


def slice_seq_record_preserve_truncated(
    seq_record: SeqRecord, slice_tuple: tuple[int, int]
) -> SeqRecord:
    sliced_rec = seq_record[slice_tuple[0] : slice_tuple[1]]
    truncated_features_left = []
    truncated_features_right = []
    for feature in seq_record.features:
        if feature.location.start < slice_tuple[0] < feature.location.end:
            feat = deepcopy(feature)
            feat.location = (
                FeatureLocation(
                    BeforePosition(slice_tuple[0]),
                    (
                        feat.location.end
                        if feat.location.end <= slice_tuple[1]
                        else AfterPosition(slice_tuple[1])
                    ),
                    feat.location.strand,
                )
                - slice_tuple[0]
            )
            feat.id = f"{feature.id}_truncated"
            feat.qualifiers["truncated"] = ["left"]
            truncated_features_left.append(feat)
            continue
        if feature.location.start < slice_tuple[1] < feature.location.end:
            feat = deepcopy(feature)
            feat.location = (
                FeatureLocation(
                    (
                        feat.location.start
                        if feat.location.start >= slice_tuple[0]
                        else BeforePosition(slice_tuple[0])
                    ),
                    AfterPosition(slice_tuple[1]),
                    feat.location.strand,
                )
                - slice_tuple[0]
            )
            feat.id = f"{feature.id}_truncated"
            feat.qualifiers["truncated"] = ["right"]
            truncated_features_right.append(feat)
    sliced_rec.features = (
        truncated_features_left + sliced_rec.features + truncated_features_right
    )
    return sliced_rec


def add_seq_to_SeqRecord_as_feature(
    seq_record: SeqRecord,
    target_seq: Seq,
    feature_type: str = "feature",
    qualifiers: dict = {},
) -> SeqRecord:
    """
    Add sequence to the end of the sequence record. Automatically finds the
    correct strand and location for the target sequence.
    If a primer bind sequence is provided, it will check if there is
    5' overhang, and avoid that.
    """

    if (
        feature_type == "primer_bind" and "overhang" in qualifiers
    ):  # 5' overhang
        try:
            overhang = int(qualifiers["overhang"][0])
            overhang = Seq("N" * overhang)
        except ValueError:
            overhang = Seq(qualifiers["overhang"][0])
            assert (
                overhang in target_seq
            ), "Overhang not found in target sequence"
        anneal_seq = target_seq[len(overhang) :]
    else:
        anneal_seq = target_seq

    seq_seq = seq_record.seq
    seq_len = len(seq_record)
    anneal_len = len(anneal_seq)
    on_strand = 1
    if anneal_seq.reverse_complement() in seq_seq:
        # located on -1 strand
        on_strand = -1
    else:
        raise ValueError("Oligo not found in sequence")
    match_seq = (
        anneal_seq if on_strand == 1 else anneal_seq.reverse_complement()
    )
    for i in range(seq_len - anneal_len + 1):
        if seq_seq[i : i + anneal_len] == match_seq:
            location = FeatureLocation(i, i + anneal_len, strand=on_strand)
            feature = SeqFeature(
                location=location,
                type=feature_type,
                qualifiers=qualifiers,
            )
            seq_record.features.append(feature)
            return seq_record
    raise ValueError("Oligo not found in sequence")


def reverse_complement_with_features(record: SeqRecord):
    new_seq = record.seq.reverse_complement()
    new_features = []
    seq_length = len(record)

    for f in record.features:
        new_start = rev_pos(seq_length, f.location.end)
        new_end = rev_pos(seq_length, f.location.start)
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


def rev_pos(seq_length, pos):
    new_pos = seq_length - pos
    if isinstance(pos, BeforePosition):
        return AfterPosition(new_pos)
    elif isinstance(pos, AfterPosition):
        return BeforePosition(new_pos)
    elif isinstance(pos, ExactPosition):
        return ExactPosition(new_pos)
    else:
        return new_pos
