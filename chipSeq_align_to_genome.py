# This file is licensed under the MIT License

import argparse
from pathlib import Path

from pyBioinfo_modules.wrappers.bowtie2 import multiple_raw_align_bowtie2

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--rawReadsDir", type=Path, help="path to raw reads")
parser.add_argument(
    "-g",
    "--genome",
    type=Path,
    help="path to genome file, can be fasta or genbank",
)
parser.add_argument(
    "-o", "--output", type=Path, help="Target directory for output bam files"
)
parser.add_argument("-p", "--processers", type=int, help="number of cpu to use")

args = parser.parse_args()
output_path = args.output
if not output_path.exists():
    output_path.mkdir()

multiple_raw_align_bowtie2(
    raw=[args.rawReadsDir],
    sampleNames=None,
    genomes=[args.genome],
    out=output_path,
    isPe=False,
    ncpu=args.processers,
    dryRun=False,
)
