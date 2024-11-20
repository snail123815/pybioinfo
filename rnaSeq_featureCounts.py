import argparse
import logging
from pathlib import Path

from pyBioinfo_modules.wrappers.featureCounts import multiple_featureCounts

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


def main():
    parser = arg_setup()
    args = parser.parse_args()
    console_handler = logging.StreamHandler()
    logger.addHandler(console_handler)
    logger.setLevel(logging.DEBUG)
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
