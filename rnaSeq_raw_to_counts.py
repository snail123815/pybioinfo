import argparse
import logging
from pathlib import Path

from Bio import SeqIO

from pyBioinfo_modules.wrappers.bowtie2 import multiple_raw_align_bowtie2
from rnaSeq_featureCounts import multiple_featureCounts

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("rnaSeq_raw_to_counts.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert raw RNA-Seq data to counts using featureCounts."
    )
    parser.add_argument(
        "--raw", required=True, type=Path, nargs="+", help="path(s) to raw data"
    )
    parser.add_argument(
        "--out", required=True, type=Path, help="path to output alignment"
    )
    parser.add_argument(
        "--gbk",
        required=True,
        type=Path,
        help="path to genome file, genbank format",
    )
    parser.add_argument(
        "--isPe", action="store_true", help="set if you have pairend"
    )
    parser.add_argument(
        "--ncpu", default=1, type=int, help="number of cpu to use"
    )
    parser.add_argument(
        "--dryRun",
        action="store_true",
        help="perform a dry run without executing Bowtie2",
    )
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
            "set it if you do not want the features: -P, Check validity of "
            "paired-end distance when counting read pairs. "
            "Use -d and -D to set thresholds; -B, only count read pairs that "
            "have both ends aligned. (must set together with -P) "
        ),
    )
    return parser.parse_args()


def main(args):
    args.out.mkdir(exist_ok=True, parents=True)
    alignment_output_path = args.out / "alignment"
    # convert gbk to fasta
    fna_file = args.out / (args.gbk.stem + ".fna")
    if not fna_file.exists():
        logger.info(f"Converting {args.gbk} to fasta")
        SeqIO.write(SeqIO.parse(args.gbk, "genbank"), fna_file, "fasta")
    logger.info(f"Genome fasta file: {fna_file}")
    logger.info(f"Alignment output path: {alignment_output_path}")
    logger.info(f"Raw data: {args.raw}")
    logger.info(f"Number of CPUs: {args.ncpu}")
    logger.info(f"Is pairend: {args.isPe}")
    logger.info("Running bowtie2 alignment")
    multiple_raw_align_bowtie2(
        raw=args.raw,
        sampleNames=None,
        out=alignment_output_path,
        genomes=[fna_file],
        isPe=args.isPe,
        ncpu=args.ncpu,
        dryRun=args.dryRun,
    )
    logger.info("Alignment done")

    counts_output_path = args.out / "counts"
    logger.info(f"Counts output path: {args.out/'counts'}")
    logger.info(f"Target feature: {args.targetFeature}")
    logger.info(f"Group factor: {args.groupFactor}")
    logger.info(f"Fraction counting: {args.fractionCounting}")
    logger.info(f"Pairend loose: {args.peLoose}")
    logger.info("Running featureCounts")
    multiple_featureCounts(
        input_path=alignment_output_path,
        output_path=counts_output_path,
        gbk_path=args.gbk,
        ncpu=args.ncpu,
        isPe=args.isPe,
        targetFeature=args.targetFeature,
        groupFactor=args.groupFactor,
        fractionCounting=args.fractionCounting,
        peLoose=args.peLoose,
        dryRun=args.dryRun,
    )
    logger.info("FeatureCounts done")


if __name__ == "__main__":
    args = parse_args()

    main(args)
