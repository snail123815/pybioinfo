import os
from pathlib import Path
from typing import Literal
import logging

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

from pyBioinfo_modules.basic.decompress import decompFileIfCompressed
from pyBioinfo_modules.bio_sequences.bio_seq_file_extensions import (
    GBK_EXTENSIONS,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
log = logging.getLogger(__name__)


def _getCdss(
    seqObj, codonTable=11, getProteins=0, getIdFrom=None
) -> list[SeqRecord]:
    """Extract proteins from a SeqRecord.
    If your input file have multiple contigs, do a loop"""
    proteins = []
    cdss = []
    idFromSequence = ["locus_tag", "label", "protein_id"]

    for i, feat in enumerate(seqObj.features):
        if feat.type == "CDS":
            cds = seqObj.seq[feat.location.start : feat.location.end]

            proteinId = ""
            proteinGeneId = ""
            proteinProduct = ""
            proteinTranslation = ""

            if getIdFrom is not None:
                if getIdFrom in feat.qualifiers:
                    proteinId = feat.qualifiers[getIdFrom][0]
                else:
                    log.warning(
                        f'Feature name "{getIdFrom}" not found in CDS'
                        f'{f"CDS_{str(i).zfill(4)}"}'
                    )
            if proteinId == "":
                proteinId = f"CDS_{str(i).zfill(4)}"
                for idFrom in idFromSequence:
                    if idFrom in feat.qualifiers:
                        proteinId = feat.qualifiers[idFrom][0]
                        break

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
                    id=proteinId,
                    name=proteinGeneId,
                    description=proteinProduct,
                )
                proteins.append(p)
            else:
                cdsRec = SeqRecord(
                    cds,
                    id=(proteinGeneId if proteinGeneId != "" else proteinId),
                    name=proteinId,
                    description=proteinProduct,
                )
                cdss.append(cdsRec)

    if getProteins:
        return proteins
    else:
        return cdss


def _getFeatureFromGbk(
    gbkPath: Path,
    codonTable=11,
    targetFeature: Literal["cds", "protein"] = "protein",
    getIdFrom=None,
) -> Path:
    """
    Extract a feature from a GenBank file and write it to a new file
    This function deals with reading file and deals with possible multiple
    contigs.
    Actual feature extraction is done by getCdss()
    """
    gbkPath, unzip = decompFileIfCompressed(gbkPath)
    try:
        if targetFeature == "protein":
            faaPath = gbkPath.with_suffix(".faa")
            proteins = []
            for s in SeqIO.parse(str(gbkPath), "genbank"):
                proteins.extend(
                    _getCdss(
                        s,
                        codonTable=codonTable,
                        getProteins=1,
                        getIdFrom=getIdFrom,
                    )
                )
            n = SeqIO.write(proteins, faaPath, "fasta")
            log.info(f"Successfully wrote {n} proteins")
            outputFile = faaPath
        elif targetFeature == "cds":
            cdss = []
            fnaPath = gbkPath.with_suffix(".cds.fna")
            for s in SeqIO.parse(str(gbkPath), "genbank"):
                cdss.extend(_getCdss(s, getIdFrom=getIdFrom))
            n = SeqIO.write(cdss, fnaPath, "fasta")
            log.info(f"Successfully wrote {n} CDSs")
            outputFile = fnaPath

    except Exception as e:
        raise e
    finally:
        if unzip:
            os.remove(str(gbkPath))
    return outputFile


def getFaaFromGbk(gbkPath: Path, codonTable=11, getIdFrom=None) -> Path:
    """
    Extract protein sequences from a GenBank file
    """
    return _getFeatureFromGbk(
        gbkPath,
        codonTable=codonTable,
        targetFeature="protein",
        getIdFrom=getIdFrom,
    )


def getCdsFromGbk(gbkPath: Path, getIdFrom=None) -> Path:
    """
    Extract CDS sequences from a GenBank file
    """
    return _getFeatureFromGbk(gbkPath, targetFeature="cds", getIdFrom=getIdFrom)
