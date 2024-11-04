# Example usage in module:
#
# from chipSeq_plot_pileup_comparisons import (
# arg_parser,
# read_input,
# plot_pileup,
# get_target_region,
# )
# import matplotlib.pyplot as plt
# from Bio import SeqIO
#
#
# fig, axs = plt.subplots(2, 1, figsize=(10, 8))
# args_list = [
# "--macsOutput",
# "~/data/Proj.Community_DAP_to_ChIP/Phase_I_testing_202408/macs3_peakcalling/comm002",
# "--genome",
# "~/data/Proj.Community_DAP_to_ChIP/M145_assembly_AL645882.gb",
# "--region",
# "1,968,245-1,969,261",
# ]
# args = arg_parser().parse_args(args_list)
# genome_with_annotation = SeqIO.read(args.genome.expanduser(), "genbank")
# tr_start, tr_end = get_target_region(args, genome_with_annotation)
# tr_control_data, tr_treat_data = read_input(args, tr_start, tr_end)
#
# plot_pileup(
# axs[0],
# tr_control_data,
# tr_treat_data,
# tr_start,
# tr_end,
# genome_with_annotation,
# do_logscale=False
# )
#
# args_list = [
# "--macsOutput",
# "~/data/Proj.Community_DAP_to_ChIP/Phase_I_testing_202408/macs3_peakcalling/comm002",
# "--genome",
# "~/data/Proj.Community_DAP_to_ChIP/M145_assembly_AL645882.gb",
# "--region",
# "1,968,245-1,969,261",
# "--logscale",
# ]
# args = arg_parser().parse_args(args_list)
# genome_with_annotation = SeqIO.read(args.genome.expanduser(), "genbank")
# plot_pileup(
# axs[1],
# tr_control_data,
# tr_treat_data,
# tr_start,
# tr_end,
# genome_with_annotation,
# do_logscale=True
# )
#
# fig.savefig("multiple_plot.png", dpi=600)

import argparse
import logging
import lzma
import re
from pathlib import Path
from typing import TextIO

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
from Bio import SeqIO
from matplotlib.patches import FancyArrow, Rectangle

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def arg_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Read pileup value (coverage) from macs output pileup file, "
            "and plot the pileup line on the given genome within a given region."
        )
    )
    parser.add_argument(
        "--macsOutput",
        type=Path,
        help="The path to the macs output dir.",
        required=True,
    )
    parser.add_argument(
        "--genome",
        type=Path,
        help="The path to the genome file.",
        required=True,
    )
    arg_group1 = parser.add_mutually_exclusive_group(required=True)
    arg_group1.add_argument("--region", type=str, help="The region to plot.")
    arg_group1.add_argument("--gene", type=str, help="The gene to plot.")
    parser.add_argument(
        "--flanking",
        type=int,
        default=None,
        help=(
            "When plotting a gene, the flanking region to include. "
            "Only effective when --gene is provided."
        ),
    )
    parser.add_argument(
        "--logscale",
        action="store_true",
        help="Whether to plot the y-axis in log scale.",
    )
    parser.add_argument(
        "--savefig",
        type=Path,
        default="./__temp.png",
        help="The path to save the plot.",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="ChIP-seq pileup comparison",
        help="The title of the plot.",
    )

    return parser


def read_macs_pileup(
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


def get_target_region(args, genome_with_annotation):
    if args.region and args.flanking:
        log.error("--flanking is only effective when --gene is provided.")
    if args.gene and args.flanking is None:
        args.flanking = 1500
        log.warning(
            f"--flanking is not provided, using default value: {args.flanking}"
        )

    if args.region:
        # Parse the region
        match = re.match(r"(\d+)-(\d+)", args.region.replace(",", ""))
        if not match:
            log.error(
                "Invalid region format. Please provide a region in the format: "
                "start-end"
            )
            raise ValueError("Invalid region format.")
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
            log.error(f"Gene {args.gene} not found in the genome file.")
        tr_start = max(0, gene_start - args.flanking - 1)
        tr_end = min(len(genome_with_annotation), gene_start + args.flanking)
        log.info(
            f"Region with flanking region to plot: {tr_start + 1}-{tr_end}"
        )
    return tr_start, tr_end


def read_input(args, tr_start, tr_end):
    """
    Reads and processes the input arguments for the chipSeq plot pileup
    comparisons.

    Args:
        args (argparse.Namespace): The parsed command-line arguments.

    Returns:
        (tr_control_data, tr_treat_data, genome_with_annotation,
         tr_start, tr_end, do_logscale)
    Raises:
        ValueError: If the region format is invalid or the gene is not found
        in the genome file.

    The function performs the following steps:
    1. Reads the genome file with annotations.
    2. If a region is specified, it parses the region and extracts the start
       and end positions.
    3. If a gene is specified, it finds the gene in the genome annotations and
       calculates the start and end positions based on the flanking region.
    """
    # Find all matching files
    pileup_files = [
        file
        for file in args.macsOutput.expanduser().glob("*.*")
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
    if treat_pileup_path.suffix == ".xz":
        treat_pileup = lzma.open(treat_pileup_path, "rt")
    else:
        treat_pileup = open(treat_pileup_path, "r")

    tr_control_data = np.array(
        read_macs_pileup(control_lambda, tr_start, tr_end)
    )
    tr_treat_data = np.array(read_macs_pileup(treat_pileup, tr_start, tr_end))
    control_lambda.close()
    treat_pileup.close()

    return tr_control_data, tr_treat_data


# Function to plot genes as arrows
def plot_genes(
    ax, genome_with_annotation, region_start, region_end, color=["C0", "C2"]
):
    xrange = ax.get_xlim()[1] - ax.get_xlim()[0]
    arrow_y_loc = 0.1
    # Gene x-axis position will be defined by gene location itself,
    # while y-axis position will be fixed to axes coordinate
    trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    for feature in genome_with_annotation.features:
        if feature.type == "gene":
            gene_start = feature.location.start
            gene_end = feature.location.end
            gene_strand = feature.location.strand

            partial_right = False
            partial_left = False
            if gene_end < region_start:
                continue
            if gene_start > region_end:
                break
            if gene_start < region_start:
                gene_start = region_start
                partial_left = True
            if gene_end > region_end - 1:
                gene_end = region_end - 1
                partial_right = True

            arrow_color = color[0] if gene_strand == 1 else color[1]

            # Define fixed head width and head length
            head_width = 0.05
            head_length = 0.007 * xrange

            # Add arrowhead
            # Make sure the arrowhead is of the same size
            arrow_start = gene_start if gene_strand == 1 else gene_end
            arrow_vector = gene_strand * (gene_end - gene_start)
            arrow = FancyArrow(
                arrow_start,
                arrow_y_loc,
                arrow_vector,
                0,
                width=head_width,
                head_width=head_width,
                head_length=head_length,
                length_includes_head=True,
                linewidth=0,
                facecolor=arrow_color,
                transform=trans,
            )
            ax.add_patch(arrow)

            def make_gene_appear_truncated(
                truncate_loc,
                direction,
                num_segments=2,
                seg_length_to_arrow_head_ratio=0.4,
            ):
                segment_length = seg_length_to_arrow_head_ratio * head_length
                for i in range(num_segments):
                    # if i == 1: continue
                    segment_start = (
                        truncate_loc + direction * 2 * i * segment_length
                    )
                    segment_start += (
                        segment_length
                        * direction
                        * (2 if direction == -1 else 1)
                    )
                    rect = Rectangle(
                        (segment_start, arrow_y_loc - head_width / 2),
                        segment_length,
                        head_width,
                        linewidth=0,
                        edgecolor=None,
                        facecolor=ax.get_facecolor(),
                        transform=trans,
                    )
                    ax.add_patch(rect)

            if partial_left:
                make_gene_appear_truncated(gene_start, 1)
            if partial_right:
                make_gene_appear_truncated(gene_end, -1)

            ax.text(
                (gene_start + gene_end) / 2,
                arrow_y_loc,
                feature.qualifiers.get("gene", [""])[0],
                ha="center",
                va="center_baseline",
                fontsize=8,
                transform=trans,
            )


def plot_pileup(
    ax,
    tr_control_data,
    tr_treat_data,
    tr_start,
    tr_end,
    genome_with_annotation,
    do_logscale,
):
    # Plot control pileup
    ax.plot(
        tr_control_data[:, 0],
        tr_control_data[:, 1],
        label="Genome-seq",
        color="silver",
    )

    # Plot treat pileup
    peaks = np.ma.masked_where(
        tr_treat_data[:, 1] < tr_control_data[:, 1], tr_treat_data[:, 1]
    )
    base_line = np.ma.masked_where(
        tr_treat_data[:, 1] > tr_control_data[:, 1], tr_treat_data[:, 1]
    )
    x = tr_treat_data[:, 0]

    ax.plot(x, base_line, label="ChIP-seq", color="dimgray")
    ax.plot(x, peaks, label="ChIP-seq peaks", color="C1")
    ax.legend()
    ax.spines["top"].set_visible(False)  # Hide the top spine
    ax.spines["right"].set_visible(False)  # Hide the right spine

    plot_genes(ax, genome_with_annotation, tr_start, tr_end)

    ax.set_xlim(tr_start, tr_end)
    ax.set_xlabel("Genomic Position")

    if do_logscale:
        ax.set_yscale("log")
        min_data = min(tr_control_data[:, 1].min(), tr_treat_data[:, 1].min())
        if min_data > 100:
            ymin = 100
        elif min_data > 10:
            ymin = 10
        else:
            ymin = 1
        ax.set_ylim(ymin, ax.get_ylim()[1])
    else:
        ymin = 0

    ax.set_ylabel("Pileup")
    # Remove y-ticks below zero
    yticks = ax.get_yticks()
    ax.set_yticks([ytick for ytick in yticks if ytick >= ymin])
    ax.spines["bottom"].set_position(("data", ymin))
    # Remove the y-axis spine below x axis
    ax.spines["left"].set_bounds(ymin, ax.get_ylim()[1])
    ax.xaxis.set_ticks_position("bottom")
    ax.xaxis.set_ticks_position("bottom")


def __main__():
    argparser = arg_parser()
    args = argparser.parse_args()

    genome_with_annotation = SeqIO.read(args.genome.expanduser(), "genbank")

    tr_start, tr_end = get_target_region(args, genome_with_annotation)
    (
        tr_control_data,
        tr_treat_data,
    ) = read_input(args, tr_start, tr_end)

    fig, ax = plt.subplots(figsize=(10, 4))

    plot_pileup(
        ax,
        tr_control_data,
        tr_treat_data,
        tr_start,
        tr_end,
        genome_with_annotation,
        args.logscale,
    )
    ax.set_title(args.title)
    fig.savefig(args.savefig, dpi=100)


if __name__ == "__main__":
    __main__()
