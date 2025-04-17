import argparse
from pathlib import Path

from pyBioinfo_modules.wrappers.macs import (macs_callPeak, macs_predictd,
                                             parse_comparison_file)

parser = argparse.ArgumentParser(description="ChipSeq MACS3 Peak Calling")

parser.add_argument(
    "-b", "--bamPath", type=Path, help="path to sorted bam files"
)
parser.add_argument("-o", "--output", type=Path, help="path to output files")
parser.add_argument(
    "-c",
    "--comparisonsFile",
    type=Path,
    help='path to a tsv file of comparisons, with headers: "name", "ctr", "exp"',
)
parser.add_argument(
    "-gs",
    "--genomeSize",
    default="8.67e6",
    help='genome size, default to "8.67e6"',
)
parser.add_argument(
    "--pairend", help="set if your data is paired ends", action="store_true"
)


arg_group1 = parser.add_mutually_exclusive_group()
arg_group1.add_argument(
    "--predictd",
    help=(
        "set if you want to use 'macs3 predictd' for fragment size estimation. "
        "This option will run predictd on all the experiments and use the "
        "average fragment size for peak calling. "
        "Imply --nomodel in peak calling. "
        " NOTE pairend reads is not supported for predictd"
    ),
    action="store_true",
)
arg_group1.add_argument(
    "--extsize",
    type=int,
    default=0,
    help=(
        "The arbitrary extension size in bp. When nomodel is true, MACS will "
        "use this value as fragment size to extend each read "
        "towards 3' end, then pile them up."
        "Exclusive with --predictd"
    ),
)
args = parser.parse_args()

if args.pairend and args.predictd:
    parser.error("--pairend and --predictd cannot be used together")

bamPath = args.bamPath
gsize = args.genomeSize
outputDir = args.output
compFile = args.comparisonsFile
doPredictd = args.predictd
isPe = args.pairend
extsize = args.extsize

experimentDict = parse_comparison_file(compFile, bamPath)
if extsize == 0:
    if doPredictd and not isPe:
        fragsize = macs_predictd(experimentDict, outputDir, gsize)
    else:
        fragsize = None
else:
    fragsize = extsize

macs_callPeak(experimentDict, outputDir, gsize, isPe=isPe, fragsize=fragsize)
