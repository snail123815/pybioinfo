import argparse
import logging
from pathlib import Path

from pyBioinfo_modules.wrappers.bowtie2 import multiple_raw_align_bowtie2

logger = logging.getLogger(__name__)


def arg_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--raw", required=True, type=Path, nargs="+", help="path(s) to raw data"
    )
    parser.add_argument(
        "--out", required=True, type=Path, help="path to output alignment"
    )
    parser.add_argument(
        "--genome",
        nargs="+",
        required=True,
        type=Path,
        help="path to genome file(s), also supports indexed genome (*.bt2)",
    )
    parser.add_argument(
        "--isPe", action="store_true", help="set if you have pairend"
    )
    parser.add_argument(
        "--pesuffix",
        nargs=2,
        help=(
            "suffix to pairend file name. Will impute by default."
            " eg. a_1_.fq.gz and a_2_.fq.gz, --pesuffix _1_ _2_"
        ),
    )
    parser.add_argument(
        "--ncpu", default=1, type=int, help="number of cpu to use"
    )
    parser.add_argument(
        "--sampleNames",
        nargs="+",
        help=(
            "List of sample names."
            " Set when sample names are random string"
            " which will unpair the samples after removing pairend suffix"
        ),
    )
    parser.add_argument(
        "--dryRun",
        action="store_true",
        help="perform a dry run without executing Bowtie2",
    )
    return parser


def main():
    parser = arg_setup()
    args = parser.parse_args()

    # Process loggers
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(args.out / "align.log"),
            logging.StreamHandler(),
        ],
    )

    logger.debug(args)
    multiple_raw_align_bowtie2(
        raw=args.raw,
        genomes=args.genome,
        out=args.out,
        isPe=args.isPe,
        sampleNames=args.sampleNames,
        ncpu=args.ncpu,
        dryRun=args.dryRun,
    )

    logger.handlers.clear()


if __name__ == "__main__":
    main()
