import time
import argparse
import logging
from pathlib import Path
from rnaSeq_Align_bowtie2 import align_multiple_raw_bowtie2
from rnaSeq_featureCounts import multiple_featureCounts

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("rnaSeq_raw_to_counts.log"),
        logging.StreamHandler()
    ]
)

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
    pass

if __name__ == "__main__":
    args = parse_args()
    align_args = args

    count_args = args
    main(args)

