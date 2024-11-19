import subprocess
from tempfile import NamedTemporaryFile
import os
import time
import argparse
import logging
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from BCBio import GFF
from pathlib import Path
from pyBioinfo_modules.wrappers._environment_settings import (
    SHORTREADS_ENV,
    SHELL,
    withActivateEnvCmd,
)

logger = logging.getLogger(__name__)


def arg_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="alignment folder (bam or sam files)",
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="gene counts output folder"
    )
    parser.add_argument(
        "--gbk", type=Path, required=True, help="genbank file with annotation"
    )
    parser.add_argument(
        "--ncpu", type=int, default=1, help="number of cpu to use"
    )
    parser.add_argument("--isPe", action="store_true", help="set if pairend")
    parser.add_argument(
        "-t", "--targetFeature", type=str, default="gene", help="target feature"
    )
    parser.add_argument(
        "-g",
        "--groupFactor",
        type=str,
        default="locus_tag",
        help="group factor",
    )
    parser.add_argument(
        "--fractionCounting",
        action="store_true",
        help="will add -M --fraction -O if True",
    )
    parser.add_argument(
        "--peLoose",
        action="store_true",
        help=(
            "Will use loose configuration in pairend counting, "
            "set it if you donot want the features: -P, Check validity of "
            "paired-end distance when counting read pairs. "
            "Use -d and -D to set thresholds; -B, only count read pairs that "
            "have both ends aligned. (must set together with -P) "
        ),
    )
    return parser


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
):
    output_path.mkdir(exist_ok=True)
    log_file = output_path / "!featureCounts.log"
    log_file_handler = logging.FileHandler(log_file)
    # console_handler = logging.StreamHandler(stream=sys.stdout)
    console_handler = logging.StreamHandler()
    log_formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    log_file_handler.setFormatter(log_formatter)
    console_handler.setFormatter(log_formatter)
    logger.addHandler(console_handler)
    logger.addHandler(log_file_handler)
    logger.setLevel(logging.DEBUG)
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
    logger.info("=" * 20 + getTime() + "=" * 20)
    files = [
        f for f in input_path.iterdir() if f.is_file() and f.suffix == ".bam"
    ]
    if len(files) == 0:
        files = [
            f
            for f in input_path.iterdir()
            if f.is_file() and f.suffix == ".sam"
        ]
    assert len(files) != 0

    for i, f in enumerate(files):
        logger.info(f"Processing {i+1}/{len(files)}: ")
        ts = time.time()
        b = os.path.splitext(os.path.split(f)[-1])[0]
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
            str(f),
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
                cmdList.insert(
                    5, "-B"
                )  # Only count read pairs that have both ends
                # aligned. (must set together with -P)
        logger.info(" ".join(cmdList))
        cmd = withActivateEnvCmd(" ".join(cmdList), SHORTREADS_ENV)
        res = subprocess.run(
            cmd, shell=True, capture_output=True, executable=SHELL
        )
        if res.returncode != 0:
            logger.info(res.stdout.decode())
            logger.info(res.stderr.decode())
            raise Exception
        logger.info(res.stderr.decode())
        logger.info(f"Finished in {diffTime(ts)}\n")
    logger.info(f"All done, time elapsed {diffTime(finalTs)}")
    logger.info("=" * 20 + getTime() + "=" * 20)


def getTime():
    return time.strftime("%z, %a, %d %b %Y, %H:%M:%S", time.localtime())


def diffTime(a):
    d = abs(time.time() - a)
    h = int(d // 3600)
    return str(h).zfill(2) + time.strftime(":%M:%S", time.gmtime(d))


def main():
    parser = arg_setup()
    args = parser.parse_args()
    multiple_featureCounts(
        input_path=args.input,
        output_path=args.output,
        gbk_path=args.gbk,
        ncpu=args.ncpu,
        isPe=args.isPe,
        targetFeature=args.targetFeature,
        groupFactor=args.groupFactor,
        fractionCounting=args.fractionCounting,
        peLoose=args.peLoose,
    )


if __name__ == "__main__":
    main()
