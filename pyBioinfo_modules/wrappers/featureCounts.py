import logging
import subprocess
import time
from pathlib import Path

from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from pyBioinfo_modules.basic.basic import getTimeStr, timeDiffStr
from pyBioinfo_modules.wrappers._environment_settings import (
    SHELL, SHORTREADS_ENV, withActivateEnvCmd)

logger = logging.getLogger(__name__)


def featureCounts(
    file: Path,
    gffFile: Path,
    ncpu: int,
    output_path: Path,
    isPe: bool,
    targetFeature: str,
    groupFactor: str,
    fractionCounting: bool,
    peLoose: bool,
    dryRun: bool = False,
):
    """
    Run the featureCounts command with specified parameters.

    Args:
        file (Path): Path to the input BAM/SAM file.
        gffFile (Path): Path to the GFF annotation file.
        ncpu (int): Number of CPUs to use.
        output_path (Path): Path to the output directory.
        isPe (bool): Whether the input is paired-end.
        targetFeature (str): Feature type to count.
        groupFactor (str): Grouping factor for counting.
        fractionCounting (bool): Whether to use fractional counting.
        peLoose (bool): Whether to use loose paired-end criteria.
        dryRun (bool): If True, only log the command without executing.
    """
    ts = time.time()
    b = file.stem
    cmdList = [
        "featureCounts",
        "-T",
        str(ncpu),
        "-a",
        str(gffFile),
        "-F",
        "GTF",
        "-t",
        targetFeature,
        "-g",
        groupFactor,
        "--minOverlap",
        "20",  # Minimum number of overlapping bases in a read that is
        # required for read assignment. 1 by default. Number of
        # overlapping bases is counted from both reads if paired
        # end. If a negative value is provided, then a gap of up
        # to specified size will be allowed between read and the
        # feature that the read is assigned to.
        "--fracOverlap",
        "0.25",  # Minimum fraction of overlapping bases in a read that is
        # required for read assignment. Value should be within range
        # [0,1]. 0 by default. Number of overlapping bases is
        # counted from both reads if paired end. Both this option
        # and '--minOverlap' option need to be satisfied for read
        # assignment.
        "-o",
        str(output_path / f"{b}.txt"),
        str(file),
    ]
    if fractionCounting:
        # Multi-mapping reads will also be counted. For a multi-
        cmdList.insert(3, "-M")
        # mapping read, all its reported alignments will be
        # counted. The 'NH' tag in BAM/SAM input is used to detect
        # multi-mapping reads.
        # Assign fractional counts to features. This option must
        cmdList.insert(4, "--fraction")
        # be used together with '-M' or '-O' or both. When '-M' is
        # specified, each reported alignment from a multi-mapping
        # read (identified via 'NH' tag) will carry a fractional
        # count of 1/x, instead of 1 (one), where x is the total
        # number of alignments reported for the same read. When '-O'
        # is specified, each overlapping feature will receive a
        # fractional count of 1/y, where y is the total number of
        # features overlapping with the read. When both '-M' and
        # '-O' are specified, each alignment will carry a fractional
        # count of 1/(x*y).
        # Assign reads to all their overlapping meta-features (or
        cmdList.insert(5, "-O")
        # features if -f is specified).

    if isPe:
        # If specified, fragments (or templates) will be counted
        cmdList.insert(3, "-p")
        # instead of reads. This option is only applicable for
        # paired-end reads; single-end reads are always counted as
        # reads.
        if not peLoose:
            # Check validity of paired-end distance when counting read
            cmdList.insert(4, "-P")
            # pairs. Use -d and -D to set thresholds.
            cmdList.insert(5, "-B")  # Only count read pairs that have both ends
            # aligned. (must set together with -P)
    logger.info(" ".join(cmdList))
    cmd = withActivateEnvCmd(" ".join(cmdList), SHORTREADS_ENV)
    if dryRun:
        logger.info("Dry run, not executing command")
        logger.info(cmd)
        return
    res = subprocess.run(cmd, shell=True, capture_output=True, executable=SHELL)
    if res.returncode != 0:
        logger.info(res.stdout.decode())
        logger.info(res.stderr.decode())
        raise Exception
    logger.info(res.stderr.decode())
    logger.info(f"Finished in {timeDiffStr(ts)}\n")


def multiple_featureCounts(
    input_path: Path,
    output_path: Path,
    gbk_path: Path,
    ncpu: int,
    isPe: bool,
    targetFeature: str,
    groupFactor: str,
    fractionCounting: bool,
    peLoose: bool,
    dryRun: bool = False,
):
    """
    Run featureCounts on multiple files in a directory.

    Args:
        input_path (Path): Path to the directory containing input BAM/SAM files.
        output_path (Path): Path to the output directory.
        gbk_path (Path): Path to the GenBank file for annotation conversion.
        ncpu (int): Number of CPUs to use.
        isPe (bool): Whether the input is paired-end.
        targetFeature (str): Feature type to count.
        groupFactor (str): Grouping factor for counting.
        fractionCounting (bool): Whether to use fractional counting.
        peLoose (bool): Whether to use loose paired-end criteria.
        dryRun (bool): If True, only log the command without executing.
    """
    output_path.mkdir(exist_ok=True)
    log_file = output_path / "featureCounts.log"
    logger.addHandler(logging.FileHandler(log_file))
    # convert gbk to gff
    gffFile = output_path / "annotation.gff"
    # keep only selected features
    seqs = []
    for seq in SeqIO.parse(gbk_path, "genbank"):
        newSeq = SeqRecord(Seq(""), id=seq.id)
        for feat in seq.features:
            if (
                feat.type in targetFeature.split(",")
                and groupFactor in feat.qualifiers
            ):
                if feat.qualifiers[groupFactor][0] == "none":
                    continue
                start = feat.location.start
                end = feat.location.end
                # Convert start and end to ExactPosition if they are not already
                # This is because featureCounts does not accept otherwise
                if not isinstance(start, ExactPosition):
                    start = ExactPosition(int(start))
                if not isinstance(end, ExactPosition):
                    end = ExactPosition(int(end))
                newFeat = SeqFeature(
                    location=FeatureLocation(
                        start=start, end=end, strand=feat.location.strand
                    ),
                    type=feat.type,
                    id=feat.id,
                )
                newFeat.qualifiers[groupFactor] = feat.qualifiers[groupFactor]
                newSeq.features.append(newFeat)
        if len(newSeq.features) == 0:
            errmsg = (
                f"Did not find any feature with type {targetFeature} "
                f"and have qualifier {groupFactor}."
            )
            logger.error(errmsg)
            raise ValueError(errmsg)
        seqs.append(newSeq)

    with gffFile.open("w+") as handle:
        GFF.write(seqs, handle)
        handle.seek(0)
        logger.info("####### Head of gff file #######")
        for i, l in enumerate(handle.readlines()):
            if i > 30:
                break
            logger.info(l.strip())
        logger.info("####### End head of gff ########")

    finalTs = time.time()
    logger.info("=" * 20 + getTimeStr() + "=" * 20)
    files = [
        f for f in input_path.iterdir() if f.is_file() and f.suffix == ".bam"
    ]
    if len(files) == 0:
        files = [
            f
            for f in input_path.iterdir()
            if f.is_file() and f.suffix == ".sam"
        ]
    if dryRun:
        if len(files) == 0:
            logger.info(
                f"Dry run: No files found in {input_path}, will mimic some."
            )
            files = [input_path / "file1.bam", input_path / "file2.bam"]
    assert len(files) != 0, input_path

    for i, f in enumerate(files):
        logger.info(f"Processing {i+1}/{len(files)}: ")
        featureCounts(
            file=f,
            gffFile=gffFile,
            ncpu=ncpu,
            output_path=output_path,
            isPe=isPe,
            targetFeature=targetFeature,
            groupFactor=groupFactor,
            fractionCounting=fractionCounting,
            peLoose=peLoose,
            dryRun=dryRun,
        )
    logger.info(f"All done, time elapsed {timeDiffStr(finalTs)}")
    logger.info("=" * 20 + getTimeStr() + "=" * 20)
    # Remove logger handler with path "out / align.log"
    for handler in logger.handlers[:]:
        if isinstance(
            handler, logging.FileHandler
        ) and handler.baseFilename == str(log_file):
            logger.removeHandler(handler)
            handler.close()
