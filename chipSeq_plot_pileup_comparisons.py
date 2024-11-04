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
import numpy as np
from Bio import SeqIO
from pyBioinfo_modules.chipseq.coverage import read_macs_pileup
from pyBioinfo_modules.bio_sequences.features_from_gbk import get_target_region
from pyBioinfo_modules.bio_sequences.plot_genes import plot_genes

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
        max_data = max(tr_control_data[:, 1].max(), tr_treat_data[:, 1].max())
        min_data = min(tr_control_data[:, 1].min(), tr_treat_data[:, 1].min())
        ax.set_yscale("log")
        if min_data > 100:
            ymin = 100
        elif min_data > 10:
            ymin = 10
        else:
            ymin = 1
        ymax = np.power(10, np.ceil(np.log10(max_data)))
        ax.set_ylim(ymin, ymax)
    else:
        ymax = ax.get_ylim()[1]
        ymin = 0

    ax.set_ylabel("Pileup")
    ax.set_ylim(ymin, ymax)
    # Remove y-ticks below zero
    yticks = [t for t in ax.get_yticks() if t >= ymin and t <= ymax]
    ax.set_yticks(yticks)
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
