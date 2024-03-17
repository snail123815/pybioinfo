import os
from pathlib import Path
from typing import Literal

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

from pyBioinfo_modules.basic.decompress import decompFileIfCompressed
from pyBioinfo_modules.bio_sequences.bio_seq_file_extensions import (
    GBK_EXTENSIONS,
)


def getCdss(seqObj, codonTable=11, getProteins=0) -> SeqRecord:
    """Extract proteins from a SeqRecord.
    If your input file have multiple contigs, do a loop"""
    proteins = []
    cdss = []

    for i, feat in enumerate(seqObj.features):
        if feat.type == "CDS":
            cds = seqObj.seq[feat.location.start : feat.location.end]

            if "locus_tag" in feat.qualifiers:
                proteinLocustag = feat.qualifiers["locus_tag"][0]
            elif "label" in feat.qualifiers:  # plasmids
                proteinLocustag = feat.qualifiers["label"][0]
            else:
                proteinLocustag = f"CDS_{str(i).zfill(4)}"

            proteinGeneId = ""
            proteinProduct = ""
            proteinTranslation = ""
            if "gene" in feat.qualifiers:
                proteinGeneId = feat.qualifiers["gene"][0]
                proteinGeneId = proteinGeneId[0].upper() + proteinGeneId[1:]
            if "product" in feat.qualifiers:
                proteinProduct = feat.qualifiers["product"][0]

            if getProteins:
                if "translation" in feat.qualifiers:
                    proteinTranslation = Seq(feat.qualifiers["translation"][0])
                else:
                    if len(cds) % 3 == 0:
                        cds = (
                            cds
                            if feat.location.strand == 1
                            else cds.reverse_complement()
                        )
                        try:
                            # Get a reliable tranlsation
                            proteinTranslation = cds.translate(
                                to_stop=True, cds=True, table=codonTable
                            )
                        except TranslationError:
                            continue
                    else:  # Ignore incomplete cds
                        continue

                p = SeqRecord(
                    proteinTranslation,
                    id=proteinLocustag,
                    name=proteinGeneId,
                    description=proteinProduct,
                )
                proteins.append(p)
            else:
                cdsRec = SeqRecord(
                    cds,
                    id=(
                        proteinGeneId
                        if proteinGeneId != ""
                        else proteinLocustag
                    ),
                    name=proteinLocustag,
                    description=proteinProduct,
                )
                cdss.append(cdsRec)

    if getProteins:
        return proteins
    else:
        return cdss


def getProteins(seqObj, codonTable=11):
    """Extract proteins from a SeqRecord.
    If your input file have multiple contigs, do a loop"""
    return getCdss(seqObj, codonTable=codonTable, getProteins=1)


def getFeatureFromGbk(
    gbkPath: Path,
    codonTable=11,
    targetFeature: Literal["cds", "protein"] = "protein",
) -> Path:
    gbkPath, unzip = decompFileIfCompressed(gbkPath)
    try:
        if targetFeature == "protein":
            faaPath = gbkPath.with_suffix(".faa")
            proteins = []
            for s in SeqIO.parse(str(gbkPath), "genbank"):
                proteins.extend(getProteins(s, codonTable=codonTable))
            n = SeqIO.write(proteins, faaPath, "fasta")
            print(f"Successfully wrote {n} proteins")
            outputFile = faaPath
        elif targetFeature == "cds":
            cdss = []
            fnaPath = gbkPath.with_suffix(".cds.fna")
            for s in SeqIO.parse(str(gbkPath), "genbank"):
                cdss.extend(getCdss(s, codonTable=codonTable))
            n = SeqIO.write(cdss, fnaPath, "fasta")
            print(f"Successfully wrote {n} CDSs")
            outputFile = fnaPath

    except Exception as e:
        raise e
    finally:
        if unzip:
            os.remove(str(gbkPath))
    return outputFile


def getFaaFromGbk(gbkPath: Path, codonTable=11) -> Path:
    return getFeatureFromGbk(
        gbkPath, codonTable=codonTable, targetFeature="protein"
    )


def getCdsFromGbk(
    gbkPath: Path,
) -> Path:
    return getFeatureFromGbk(gbkPath, targetFeature="cds")
