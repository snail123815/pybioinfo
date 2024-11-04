import logging
import os
import re
from pathlib import Path
from typing import Literal

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pyBioinfo_modules.basic.decompress import decompFileIfCompressed
from pyBioinfo_modules.bio_sequences.bio_seq_file_extensions import \
    GBK_EXTENSIONS

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


def get_target_region(args, genome_with_annotation):
    """
    Determine the target region based on the provided arguments.
    This function calculates the start and end positions of a target region
    either based on a specified region or a gene with optional flanking regions.
    Args:
        args: An object containing the following attributes:
            - region (str): A string specifying the region in the format "start-end".
            - gene (str): The name of the gene to center the region around.
            The above two arguments are mutually exclusive.
            - flanking (int, optional): The number of base pairs to include on
              either side of the gene.
        genome_with_annotation: A BioRecord object containing annotated features.
    Returns:
        tuple: A tuple containing the start and end positions (0-based) of the
        target region.
    Raises:
        ValueError: If the region format is invalid.
    """

    if args.region and args.flanking:
        log.error("--flanking is only effective when --gene is provided.")
    if args.gene and args.flanking is None:
        args.flanking = 1500
        log.warning(
            f"--flanking is not provided, using default value: {args.flanking}"
        )

    if args.region:
        # Parse the region
        match = re.match(r"(\d+)-(\d+)", args.region.replace(",", ""))
        if not match:
            log.error(
                "Invalid region format. Please provide a region in the format: "
                "start-end"
            )
            raise ValueError("Invalid region format.")
        tr_start, tr_end = match.groups()
        tr_start, tr_end = int(tr_start) - 1, int(tr_end)
        log.info(f"Region to plot: {tr_start + 1}-{tr_end}")
    elif args.gene:
        # Find the start of the gene
        for feature in genome_with_annotation.features:
            if (
                feature.type == "gene"
                and feature.qualifiers["gene"][0] == args.gene
            ):
                if feature.location.strand == -1:
                    gene_start = int(feature.location.end)
                else:
                    gene_start = int(feature.location.start)
                log.info(
                    f"Gene {args.gene} found on strand "
                    f"{feature.location.strand}, position: {gene_start}"
                )
                break
        else:
            log.error(f"Gene {args.gene} not found in the genome file.")
        tr_start = max(0, gene_start - args.flanking - 1)
        tr_end = min(len(genome_with_annotation), gene_start + args.flanking)
        log.info(
            f"Region with flanking region to plot: {tr_start + 1}-{tr_end}"
        )
    return tr_start, tr_end
