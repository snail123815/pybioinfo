import argparse
import logging
import lzma
import re
from pathlib import Path
from typing import TextIO

import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.patches import Arrow, Rectangle
from matplotlib.colors import LinearSegmentedColormap

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

# Parse the input arguments
agrparser = argparse.ArgumentParser(
    description=(
        "Read pileup value (coverage) from macs output pileup file, "
        "and plot the pileup line on the given genome within a given region."
    )
)
agrparser.add_argument(
    "--macsOutput",
    type=Path,
    help="The path to the macs output dir.",
    required=True,
)
agrparser.add_argument(
    "--genome", type=Path, help="The path to the genome file.", required=True
)
arg_group1 = agrparser.add_mutually_exclusive_group(required=True)
arg_group1.add_argument("--region", type=str, help="The region to plot.")
arg_group1.add_argument("--gene", type=str, help="The gene to plot.")
agrparser.add_argument(
    "--flanking",
    type=int,
    default=None,
    help=(
        "When plotting a gene, the flanking region to include. "
        "Only effective when --gene is provided."
    ),
)

# Parse the arguments
args = agrparser.parse_args()
if args.region and args.flanking:
    agrparser.error("--flanking is only effective when --gene is provided.")
if args.gene and args.flanking is None:
    args.flanking = 1500
    log.warning(
        f"--flanking is not provided, using default value: {args.flanking}"
    )

genome_with_annotation = SeqIO.read(args.genome, "genbank")

if args.region:
    # Parse the region
    match = re.match(r"(\d+)-(\d+)", args.region.replace(",", ""))
    if not match:
        agrparser.error(
            "Invalid region format. Please provide a region in the format: "
            "start-end"
        )
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
        agrparser.error(f"Gene {args.gene} not found in the genome file.")
    tr_start = max(0, gene_start - args.flanking - 1)
    tr_end = min(len(genome_with_annotation), gene_start + args.flanking)
    log.info(f"Region with flanking region to plot: {tr_start + 1}-{tr_end}")


# Find all matching files
pileup_files = [
    file
    for file in args.macsOutput.glob("*.*")
    if file.name.endswith(".bdg") or file.name.endswith(".bdg.xz")
]

assert len(pileup_files) == 2, (
    "There should be two and only two .bdg files in the macs output dir."
    f" Found {len(pileup_files)} files." + str(pileup_files)
)

control_lambda_path = None
treat_pileup_path = None

for f in pileup_files:
    if "_control_lambda" in f.name:
        control_lambda_path = f
    elif "_treat_pileup" in f.name:
        treat_pileup_path = f
    else:
        raise ValueError(
            "Unknown file: "
            + f
            + ". Expected either _control_lambda or _treat_pileup"
        )

assert control_lambda_path, (
    "Control pileup file not found. "
    "Please make sure the file name contains '_control_lambda'."
)
assert treat_pileup_path, (
    "Experimental pileup file not found. "
    "Please make sure the file name contains '_treat_pileup'."
)

if control_lambda_path.suffix == ".xz":
    control_lambda = lzma.open(control_lambda_path, "rt")
else:
    control_lambda = open(control_lambda_path, "r")


def read_pileup(
    pileup: TextIO, tr_start: int, tr_end: int
) -> list[tuple[int, float]]:
    # Read the pileup file, a range per line
    tr_range_pileup = []
    for line in pileup:
        _, start, end, value = line.strip().split("\t")
        end = int(end)
        if end <= tr_start:
            continue
        start = int(start)
        if start > tr_end:
            break
        value = float(value)
        tr_range_pileup.append((start, end, value))
    # Fill in the the range per base
    tr_perbase_pileup = []
    for start, end, value in tr_range_pileup:
        effective_range = range(max(start, tr_start), min(end, tr_end + 1))
        for i in effective_range:
            tr_perbase_pileup.append((i, value))
    return tr_perbase_pileup


tr_control_data = np.array(read_pileup(control_lambda, tr_start, tr_end))


if treat_pileup_path.suffix == ".xz":
    treat_pileup = lzma.open(treat_pileup_path, "rt")
else:
    treat_pileup = open(treat_pileup_path, "r")

tr_treat_data = np.array(read_pileup(treat_pileup, tr_start, tr_end))


# Function to plot genes as arrows with gradient fill
def plot_genes(ax, genome_with_annotation, region_start, region_end):
    height = ax.get_ylim()[1] * 0.5
    for feature in genome_with_annotation.features:
        if feature.type == "gene":
            start = feature.location.start
            end = feature.location.end
            strand = feature.location.strand

            partial_right = False
            partial_left = False
            if end < region_start:
                continue
            if start > region_end:
                break
            if start < region_start:
                start = region_start
                partial_left = True
            if end > region_end - 1:
                end = region_end - 1
                partial_right = True

            color = "C3"
            if strand == -1:
                start, end = end, start
                color = "C4"

            # Add arrowhead
            arrow = Arrow(
                start,
                height,
                end - start + 1,
                0,
                width=height,
                # edgecolor=None,
                facecolor=color,
            )
            ax.add_patch(arrow)

            # Divide the arrow into segments to simulate gradient
            num_segments = 3
            segment_length = 0.005 * (ax.get_xlim()[1] - ax.get_xlim()[0])
            for i in range(num_segments):
                segment_start = start + strand * 2 * i * segment_length
                segment_start += (
                    segment_length * strand * (2 if strand == -1 else 1)
                )
                rect = Rectangle(
                    (segment_start, height - 1000),
                    segment_length,
                    2000,
                    linewidth=0,
                    edgecolor=None,
                    facecolor="C0",
                )
                ax.add_patch(rect)

            ax.text(
                (start + end) / 2,
                height + 0.1,
                feature.qualifiers.get("gene", [""])[0],
                ha="center",
                va="bottom",
                fontsize=8,
            )


fig, ax = plt.subplots(figsize=(10, 2))

ax.plot(tr_control_data[:, 0], tr_control_data[:, 1], label="Control")
ax.plot(tr_treat_data[:, 0], tr_treat_data[:, 1], label="Treat")
plot_genes(ax, genome_with_annotation, tr_start, tr_end)

ax.set_xlim(tr_start, tr_end)
ax.set_yticks([])
ax.set_xlabel("Genomic Position")
ax.set_title("Normalized pileup of ChIP-seq reads")

plt.tight_layout()
plt.savefig("__temp.png", dpi=600)

if control_lambda:
    control_lambda.close()
if treat_pileup:
    treat_pileup.close()
