from collections.abc import Iterable
from pathlib import Path
from typing import Literal, NotRequired, TypedDict

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord

from pyBioinfo_modules.bio_sequences.bio_features import seqFeature_to_tuple


class VarianceData(TypedDict):
    CHROM: str
    POS: int
    ID: str
    REF: str
    ALT: str
    QUAL: float
    FILTER: str
    INFO: str
    FORMAT: NotRequired[str]
    SAMPLES: NotRequired[dict[str, str]]


essentialVcfColumns = [
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
]


def vcfParser(vcf: Path) -> list[VarianceData]:
    varianceDatas: list[VarianceData] = []
    with open(vcf, "r") as vcfFh:
        getFormat: bool
        # Parse header lines
        for line in vcfFh:
            if line.startswith("##"):  # header lines
                continue
            elif line.startswith("#"):  # column names
                columns = line.strip()[1:].split("\t")
                break
        # Assert correct format
        assert "columns" in locals(), f"Header line not found for {vcf}"
        assert (
            len(columns) >= 8
        ), f"Header is not complete for {vcf}, {len(columns)} columns."
        assert (
            columns[:8] == essentialVcfColumns
        ), f"The file {vcf} is not a standard VCF file, \n{columns}"
        getFormat = len(columns) > 8 and columns[8] == "FORMAT"
        if getFormat:
            samplesColumns = columns[9:]
            assert len(samplesColumns) > 0, (
                f"No sample is found in {vcf}" " but FORMAT field present."
            )
        # Parse data
        for line in vcfFh:
            fields = line.strip().split("\t")
            varianceData = VarianceData(
                CHROM=fields[0],
                POS=int(fields[1]),
                ID=fields[2],
                REF=fields[3],
                ALT=fields[4],
                QUAL=(0 if fields[5] == "." else float(fields[5])),
                FILTER=fields[6],
                INFO=fields[7],
            )
            if getFormat:
                varianceData["FORMAT"] = fields[8]
                varianceData["SAMPLES"] = {
                    s: fields[9 + i] for i, s in enumerate(samplesColumns)
                }
            varianceDatas.append(varianceData)
    return varianceDatas


def filterVarianceData(
    varianceDatas: list[VarianceData],
    locationRange: tuple[int, int],
    quality: float,
    filter: bool,  # on FILTER column (PASS or not)
    # on INFO column, specific keys -> values
    info: dict[str, float | int | str],
    additional: dict[str, dict[str, float | int | str]],  # on SAMPLEs columns
    # {sampleName: {filter dict}}
) -> list[VarianceData]:
    """
    TODO
    filter: bool, # on FILTER column (PASS or not)
    info: dict[str, float|int|str], # on INFO column, specific keys -> values
    additional: dict[str, dict[str, float|int|str]], # on SAMPLEs columns
    # {sampleName: {filter dict}}

    Is this function worth making? bcftools actually have filter tool.
    """
    return varianceDatas


def applyVarianceDataOnSeqRecord(
    varianceDatas: list[VarianceData], seqRecord: SeqRecord
) -> SeqRecord:
    """
    Prokaryotes only: Intron and exon are not considered in this function.
    Multiple records with the same POS might break or produce false result
    """
    if not all(
        varianceData["CHROM"] == seqRecord.id for varianceData in varianceDatas
    ):
        noMatch = [
            varianceData["CHROM"]
            for varianceData in varianceDatas
            if varianceData["CHROM"] != seqRecord.id
        ][0]
        raise AssertionError(
            f"Source sequence ID is {seqRecord.id}, it does not match "
            f"in .vcf file: {noMatch}"
        )
    assert all(
        seqRecord.seq[varianceData["POS"] - 1].lower()
        == varianceData["REF"][0].lower()
        for varianceData in varianceDatas
    ), "Reference nucleotide do not match variance data"

    newSeqRecord = seqRecord[:]
    # https://www.insdc.org/submitting-standards/feature-table/#7.3
    excludeQualifiers = [
        "codon_start",
        "altitude",
        "anticodon",
        "artificial_location",
        "bio_material",
        "bound_moiety",
        "cell_line",
        "cell_type",
        "chromosome",
        "circular_RNA",
        "citation",
        "clone_lib",
        "collected_by",
        "collection_date",
        "compare",
        "country",
        "cultivar",
        "culture_collection",
        "db_xref",
        "dev_stage",
        "direction",
        "EC_number",
        "ecotype",
        "environmental_sample",
        "estimated_length",
        "exception",
        "experiment",
        "focus",
        "frequency",
        "gap_type",
        "germline",
        "haplogroup",
        "haplotype",
        "host",
        "identified_by",
        "inference",
        "isolate",
        "isolation_source",
        "lab_host",
        "lat_lon",
        "linkage_evidence",
        "macronuclear",
        "map",
        "mating_type",
        "metagenome_source",
        "mobile_element_type",
        "mod_base",
        "mol_type",
        "number",
        "operon",
        "organelle",
        "pop_variant",
        "ribosomal_slippage",
        "rpt_family",
        "rpt_type",
        "rpt_unit_range",
        "satellite",
        "serotype",
        "serovar",
        "sex",
        "specimen_voucher",
        "strain",
        "sub_clone",
        "sub_species",
        "sub_strain",
        "tag_peptide",
        "tissue_lib",
        "tissue_type",
        "transgenic",
        "translation",
        "transl_except",
        "transl_table",
        "trans_splicing" "type_material",
        "variety",
    ]

    # Valid VCF file do not need this step. But who knows.
    varianceDatas = sorted(varianceDatas, key=lambda vd: vd["POS"])

    for i, varianceData in enumerate(varianceDatas):
        if i > 0:
            assert varianceDatas[i - 1]["POS"] < varianceData["POS"], (
                "Multiple records at the same position is not supported\n"
                f"{varianceData}\n{varianceDatas[i-1]}"
            )
        assert "," not in varianceData["ALT"], (
            "Mixed VCF record is not supported,"
            f" consider remove them before applying.\n"
            f"{varianceData}"
        )
        allFeatures_startSorted = sorted(
            newSeqRecord.features, key=lambda f: f.location.start
        )
        allFeatures_endSorted = sorted(
            newSeqRecord.features, key=lambda f: f.location.end
        )
        feattuples_startSorted = [
            seqFeature_to_tuple(feat) for feat in allFeatures_startSorted
        ]
        feattuples_endSorted = [
            seqFeature_to_tuple(feat) for feat in allFeatures_endSorted
        ]
        start = varianceData["POS"] - 1 + (len(newSeqRecord) - len(seqRecord))
        assert newSeqRecord[start].lower() == varianceData["REF"][0].lower()
        length = len(varianceData["REF"])
        upStreamSeq = newSeqRecord[:start]
        nFeaturesUpStream = len(upStreamSeq.features)
        downStreamSeq = newSeqRecord[start + length :]
        nFeaturesDownStream = len(downStreamSeq.features)
        newSeqRecord = upStreamSeq + Seq(varianceData["ALT"]) + downStreamSeq
        posDiff = len(varianceData["ALT"]) - len(varianceData["REF"])

        # Magick:
        if nFeaturesDownStream == 0:
            upStreamFeattuples = set(feattuples_startSorted)
        else:
            upStreamFeattuples = set(
                feattuples_startSorted[:-nFeaturesDownStream]
            )
        downStreamFeattuples = set(feattuples_endSorted[nFeaturesUpStream:])
        affectedFeature_indexs = [
            feattuples_startSorted.index(feattuple)
            for feattuple in list(
                upStreamFeattuples.intersection(downStreamFeattuples)
            )
        ]
        affectedFeatures = [
            allFeatures_startSorted[i] for i in affectedFeature_indexs
        ]

        # add affected features back
        mutFeatures = []
        for feat in affectedFeatures:
            newLoc = FeatureLocation(
                feat.location.start,
                feat.location.end + posDiff,
                feat.location.strand,
            )
            newFeat = feat
            newFeat.location = newLoc

            for qualifier, contents in newFeat.qualifiers.items():
                if qualifier == "translation" and newFeat.type == "CDS":
                    newDna = newSeqRecord[
                        newFeat.location.start : newFeat.location.end
                    ]
                    if newFeat.location.strand == -1:
                        newDna = newDna.reverse_complement()
                    try:
                        translateTable: int | Literal["Standard"] = int(
                            newFeat.qualifiers["transl_table"][0]
                        )
                    except IndexError:
                        translateTable = "Standard"
                    newTranslate = (
                        newDna[: len(newDna) // 3 * 3]
                        .translate(translateTable, cds=False, stop_symbol="*")
                        .seq
                    )
                    if "*" in newTranslate:
                        newTranslate = newTranslate.split("*")[0]
                    else:
                        newTranslate += "..."
                    if newTranslate[0].lower() == "v" and translateTable == 11:
                        if newTranslate[0].islower():
                            newTranslate = "m" + newTranslate[1:]
                        else:
                            newTranslate = "M" + newTranslate[1:]

                    newFeat.qualifiers[qualifier] = [
                        newTranslate,
                    ]
                # white list
                if qualifier in excludeQualifiers:
                    continue
                if any(c.startswith("MUTATED_") for c in contents):
                    newFeat.qualifiers[qualifier] = contents
                else:
                    newFeat.qualifiers[qualifier] = [
                        "MUTATED_" + c for c in contents
                    ]
            mutFeatures.append(newFeat)

        newSeqRecord.features.extend(mutFeatures)
        assert len(newSeqRecord.features) == len(seqRecord.features)
    newSeqRecord.features = sorted(
        sorted(
            newSeqRecord.features,
            key=lambda feat: feat.type,
            reverse=True,  # gene ahead of CDS
        ),
        key=lambda feat: feat.location.start,
    )
    return newSeqRecord


def applyVariancesOnSeqRecords(
    varianceDatas: list[VarianceData], seqRecords: Iterable[SeqRecord]
) -> dict[str, SeqRecord]:
    """
    Prokaryotes only: Intron and exon are not considered in this function.
    """
    variantSeqRecordDict: dict[str, SeqRecord] = {}
    print(f"There are {len(varianceDatas)} to apply.")
    totalN = 0
    for seqRec in seqRecords:
        seqid = seqRec.id
        varianceDatas_specific = [
            varianceData
            for varianceData in varianceDatas
            if varianceData["CHROM"] == seqid
        ]
        print(
            f"Will apply {len(varianceDatas_specific)} variances to " f"{seqid}"
        )
        variantSeqRecordDict[seqid] = applyVarianceDataOnSeqRecord(
            varianceDatas_specific, seqRec
        )
        totalN += len(varianceDatas_specific)
    print(f"Applied {totalN} variances.")
    return variantSeqRecordDict
