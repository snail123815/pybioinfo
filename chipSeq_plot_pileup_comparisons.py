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
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from pyBioinfo_modules.bio_sequences.features_from_gbk import get_target_region
from pyBioinfo_modules.bio_sequences.plot_genes import plot_genes
from pyBioinfo_modules.chipseq.coverage import read_macs_pileup

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
    arg_group1.add_argument(
        "--region",
        type=str,
        help=(
            "The region to plot."
            "Format: 'start-end', 'start,end', or just one location."
        ),
    )
    arg_group1.add_argument("--gene", type=str, help="The gene to plot.")
    arg_group1.add_argument(
        "--peak_list",
        type=Path,
        help=("The path to the peak list ',xls' (tsv format) output by macs."),
    )
    parser.add_argument(
        "--flanking",
        type=int,
        default=None,
        help=("When plotting a gene, the flanking region to include. "),
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


def plot_macs_pileup(
    ax: plt.Axes,
    tr_control_data: np.ndarray,
    tr_treat_data: np.ndarray,
    tr_start: int,
    tr_end: int,
    do_logscale: bool = False,
    genome_with_annotation: SeqRecord | None = None,
) -> None:
    """
    Plots the MACS pileup data for control and treatment samples on the
    given axis.

    Parameters:
    ax (matplotlib.axes.Axes): The axis to plot on.
    tr_control_data (numpy.ndarray): The control data array with genomic positions and pileup values.
    tr_treat_data (numpy.ndarray): The treatment data array with genomic positions and pileup values.
    tr_start (int): The start position of the genomic region to plot.
    tr_end (int): The end position of the genomic region to plot.
    do_logscale (bool, optional): Whether to use a logarithmic scale for the y-axis. Default is False.
    genome_with_annotation (dict, optional): A dictionary containing genome annotation data. Default is None.

    Returns:
    None
    """
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

    if genome_with_annotation:
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


def read_peak_list(peak_list: Path) -> Iterable[dict]:
    with peak_list.open("rt") as f:
        for line in f:
            if (
                line.startswith("#")
                or line.startswith("chr")
                or len(line.strip()) == 0
            ):
                continue
            try:
                (
                    _,
                    start,
                    end,
                    length,
                    summit,
                    pileup,
                    _,
                    fold_enrichment,
                    _,
                    name,
                ) = line.strip().split("\t")
            except ValueError:
                print(line)
                raise

            yield {
                "start": int(start),
                "end": int(end),
                "length": int(length),
                "summit": int(summit),
                "pileup": float(pileup),
                "fold_enrichment": float(fold_enrichment),
                "name": str(name),
            }


def __main__():
    argparser = arg_parser()
    args = argparser.parse_args()

    if args.peak_list:
        peaks = list(read_peak_list(args.peak_list))

    genome_with_annotation = SeqIO.read(args.genome.expanduser(), "genbank")

    if not args.peak_list:
        tr_start, tr_end = get_target_region(
            args.region, args.gene, args.flanking, genome_with_annotation
        )
        tr_control_data, tr_treat_data = read_macs_pileup(
            args.macsOutput, tr_start, tr_end
        )

        fig, ax = plt.subplots(figsize=(10, 4))

        plot_macs_pileup(
            ax,
            tr_control_data,
            tr_treat_data,
            tr_start,
            tr_end,
            do_logscale=args.logscale,
            genome_with_annotation=genome_with_annotation,
        )
        ax.set_title(args.title)
        fig.savefig(args.savefig, dpi=100)
    else:
        for peak in peaks:
            location = peak["summit"]
            start, end = peak["start"], peak["end"]
            length = peak["length"]
            # Expand both side by length of "summit to edge plus peak length"
            # or 1500 bp if that value < 1500
            flanking = max(1500, location - start-length, end-location -length)
            tr_start, tr_end = get_target_region(
                location,
                None,
                flanking,
                genome_with_annotation,
            )
            tr_control_data, tr_treat_data = read_macs_pileup(
                args.macsOutput, tr_start, tr_end
            )
            fig, ax = plt.subplots(figsize=(10, 4))

            plot_macs_pileup(
                ax,
                tr_control_data,
                tr_treat_data,
                tr_start,
                tr_end,
                do_logscale=args.logscale,
                genome_with_annotation=genome_with_annotation,
            )
            if "_temp" not in args.savefig.name:
                savefig = args.savefig.with_name(
                    args.savefig.stem + f"_{peak['name']}.png"
                )
            else:
                savefig = Path(f"./{peak['name']}.png")
            if args.title:
                title = f"{args.title} - {peak['name']}"
            log.info(f"Saving plot to {savefig}")
            ax.set_title(title)
            fig.savefig(savefig, dpi=100)
            plt.close(fig)



if __name__ == "__main__":
    __main__()
