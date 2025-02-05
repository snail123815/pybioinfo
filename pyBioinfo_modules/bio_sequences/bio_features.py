from math import ceil

from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    ExactPosition,
    FeatureLocation,
    SeqFeature,
)
from Bio.SeqRecord import SeqRecord


def getSpanFetures(genome, startOri, endOri, expand=20000):
    newStart = max(0, startOri - expand)
    newEnd = min(len(genome), endOri + expand)
    sourceSeq = genome[newStart:newEnd]

    start = startOri - newStart
    end = start + (endOri - startOri)

    spanFeats = []
    for feat in sourceSeq.features:
        spanStart = start in feat
        spanEnd = end in feat
        # feat.__contains__(self, value)
        # Check if an integer position is within the feature.

        spanFeat = SeqFeature(type=feat.type)

        if spanStart and spanEnd:
            # Target position is inside feature
            if feat.type == "CDS":
                # calculate correct start and end location to make it inframe
                newStart = (3 - abs(start - feat.location.start)) % 3
                newEnd = end - start - abs(end - feat.location.start) % 3
            else:
                newStart = 0
                newEnd = end - start
            spanFeat.location = FeatureLocation(
                newStart, newEnd, strand=feat.location.strand
            )
            for key in feat.qualifiers:
                if key in [
                    "gene_synonym",
                    "locus_tag",
                    "product",
                    "protein_id",
                    "db_xref",
                    "mol_type",
                ]:
                    spanFeat.qualifiers[key] = [
                        f"{keyStr} (slice)" for keyStr in feat.qualifiers[key]
                    ]
                elif key == "translation":
                    spanFeat.qualifiers[key] = list(feat.qualifiers[key])
                    if feat.location.strand == 1:
                        cutPointA = ceil((start - feat.location.start) / 3)
                        cutPointB = (end - feat.location.start) // 3
                    else:
                        len(spanFeat.qualifiers[key][0])
                        cutPointA = (
                            len(spanFeat.qualifiers[key][0])
                            + 1
                            - (end - feat.location.start) // 3
                        )
                        cutPointB = (
                            len(spanFeat.qualifiers[key][0])
                            + 1
                            - ceil((start - feat.location.start) / 3)
                        )
                    spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][0][
                        cutPointA:cutPointB
                    ]
                else:
                    spanFeat.qualifiers[key] = feat.qualifiers[key]
            spanFeats.append(spanFeat)

        elif spanStart:
            # Start position inside feature, feature ends in this range
            if feat.type == "CDS":
                newStart = (3 - abs(start - feat.location.start)) % 3
            else:
                newStart = 0
            newEnd = feat.location.end - start
            spanFeat.location = FeatureLocation(
                newStart, newEnd, strand=feat.location.strand
            )
            for key in feat.qualifiers:
                if key in [
                    "gene_synonym",
                    "locus_tag",
                    "product",
                    "protein_id",
                    "db_xref",
                    "mol_type",
                ]:
                    spanFeat.qualifiers[key] = [
                        f"{keyStr} (right part)"
                        for keyStr in feat.qualifiers[key]
                    ]
                elif key == "translation":
                    spanFeat.qualifiers[key] = [
                        keyStr for keyStr in feat.qualifiers[key]
                    ]
                    if feat.location.strand == 1:
                        cutPoint = (
                            len(spanFeat.qualifiers[key][0])
                            - ceil(len(spanFeat) / 3)
                            + 1
                        )
                        spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][
                            0
                        ][cutPoint:]
                    else:
                        cutPoint = len(spanFeat) // 3
                        spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][
                            0
                        ][:cutPoint]
                else:
                    spanFeat.qualifiers[key] = feat.qualifiers[key]
            spanFeats.append(spanFeat)

        elif spanEnd:
            # End position inside feature, feature ends in this range
            if feat.type == "CDS":
                newEnd = end - start - abs(end - feat.location.start) % 3
            else:
                newEnd = end - start
            newStart = feat.location.start - start
            spanFeat.location = FeatureLocation(
                newStart, newEnd, strand=feat.location.strand
            )
            for key in feat.qualifiers:
                if key in [
                    "gene_synonym",
                    "locus_tag",
                    "product",
                    "protein_id",
                    "db_xref",
                    "mol_type",
                ]:
                    spanFeat.qualifiers[key] = [
                        f"{keyStr} (left part)"
                        for keyStr in feat.qualifiers[key]
                    ]
                elif key == "translation":
                    spanFeat.qualifiers[key] = list(feat.qualifiers[key])
                    if feat.location.strand == 1:
                        cutPoint = len(spanFeat) // 3
                        spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][
                            0
                        ][:cutPoint]
                    else:
                        cutPoint = ceil((len(feat) - len(spanFeat)) / 3)
                        spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][
                            0
                        ][cutPoint:]
                else:
                    spanFeat.qualifiers[key] = feat.qualifiers[key]
            spanFeats.append(spanFeat)

        else:
            # Not in range, ignore
            continue

    return spanFeats


# getSpanFetures


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
