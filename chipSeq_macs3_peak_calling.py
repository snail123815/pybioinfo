import argparse
from pathlib import Path
from pyBioinfo_modules.chipseq.macs3 import readComps, predictd, callPeak

parser = argparse.ArgumentParser()
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
    "--predictd",
    help=(
        "set if you want to use 'macs3 predictd' for fragment size estimation."
        " NOTE pairend reads is not supported for predictd"
    ),
    action="store_true",
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


args = parser.parse_args()
bamPath = args.bamPath
gsize = args.genomeSize
outputDir = args.output
compFile = args.comparisonsFile
doPredictd = args.predictd
isPe = args.pairend

experimentDict = readComps(compFile, bamPath)
if doPredictd and not isPe:
    fragsize = predictd(experimentDict, outputDir, gsize)
else:
    fragsize = None
callPeak(experimentDict, outputDir, gsize, isPe=isPe, fragsize=fragsize)
