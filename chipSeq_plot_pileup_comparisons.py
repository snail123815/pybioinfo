"""Plot MACS pileup comparison between control and treatment samples.
This file is licensed under the MIT License

Example usage in module:

from chipSeq_plot_pileup_comparisons import (
arg_parser,
read_input,
plot_pileup,
get_target_region,
)
import matplotlib.pyplot as plt
from Bio import SeqIO


fig, axs = plt.subplots(2, 1, figsize=(10, 8))
args_list = [
"--macsOutput",
"~/data/Proj.Community_DAP_to_ChIP/Phase_I_testing_202408/macs3_peakcalling/comm002",
"--genome",
"~/data/Proj.Community_DAP_to_ChIP/M145_assembly_AL645882.gb",
"--region",
"1,968,245-1,969,261",
]
args = arg_parser().parse_args(args_list)
genome_with_annotation = SeqIO.read(args.genome.expanduser(), "genbank")
tr_start, tr_end = get_target_region(args, genome_with_annotation)
tr_control_data, tr_treat_data = read_input(args, tr_start, tr_end)

plot_pileup(
axs[0],
tr_control_data,
tr_treat_data,
tr_start,
tr_end,
genome_with_annotation,
do_logscale=False
)

args_list = [
"--macsOutput",
"~/data/Proj.Community_DAP_to_ChIP/Phase_I_testing_202408/macs3_peakcalling/comm002",
"--genome",
"~/data/Proj.Community_DAP_to_ChIP/M145_assembly_AL645882.gb",
"--region",
"1,968,245-1,969,261",
"--logscale",
]
args = arg_parser().parse_args(args_list)
genome_with_annotation = SeqIO.read(args.genome.expanduser(), "genbank")
plot_pileup(
axs[1],
tr_control_data,
tr_treat_data,
tr_start,
tr_end,
genome_with_annotation,
do_logscale=True
)

fig.savefig("multiple_plot.png", dpi=600)
"""

import argparse
import logging
import subprocess
import tempfile
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from matplotlib.axes import Axes
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image, UnidentifiedImageError

from pyBioinfo_modules.bio_sequences.features_from_gbk import get_target_region
from pyBioinfo_modules.bio_sequences.plot_genes import plot_genes
from pyBioinfo_modules.chipseq.coverage import read_macs_pileup
from pyBioinfo_modules.chipseq.read_peak_file import read_peak_file

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def arg_parser():
    """Argument parser for the script chipSeq_plot_pileup_comparisons.py"""
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
    arg_group1 = parser.add_mutually_exclusive_group(required=False)
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
        "--ymax",
        type=int,
        default=None,
        help="The maximum y-axis value for the pileup plot.",
    )
    parser.add_argument(
        "--ymin",
        type=int,
        default=None,
        help="The minimum y-axis value for the pileup plot.",
    )
    parser.add_argument(
        "--peak_alpha",
        type=float,
        default=1.0,
        help="The alpha (transparency) value for the peak plot.",
    )
    parser.add_argument(
        "--no_control_mask",
        action="store_true",
        help=(
            "If set, do not mask treatment values below control values. "
            "This is useful for comparing two samples."
        ),
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
    parser.add_argument(
        "--consolidate",
        action="store_true",
        help=(
            "Generate a consolidated PDF file with all plots and a ZIP file "
            "with all PNG files. "
            "Removes intermediate PNG files after consolidation."
        ),
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=100,
        help="The DPI (dots per inch) for the output image.",
    )

    return parser


def plot_macs_pileup(
    ax: Axes,
    tr_control_data: np.ndarray,
    tr_treat_data: np.ndarray,
    tr_start: int,
    tr_end: int,
    do_logscale: bool = False,
    ymax_given: int | None = None,
    ymin_given: int | None = None,
    no_control_mask: bool = False,
    peak_alpha: float = 1.0,
    genome_with_annotation: SeqRecord | None = None,
) -> None:
    """
    Plots the MACS pileup data for control and treatment samples on the
    given axis.

    Parameters:
    ax (matplotlib.axes.Axes):
        The axis to plot on.
    tr_control_data (numpy.ndarray):
        The control data array with genomic positions and pileup values.
    tr_treat_data (numpy.ndarray):
        The treatment data array with genomic positions and pileup values.
    tr_start (int):
        The start position of the genomic region to plot.
    tr_end (int):
        The end position of the genomic region to plot.
    do_logscale (bool, optional):
        Whether to use a logarithmic scale for the y-axis. Default is False.
    genome_with_annotation (dict, optional):
        A dictionary containing genome annotation data. Default is None.

    Returns:
    None
    """
    # Plot control pileup
    # Color settings
    peak_color = "C1"
    if no_control_mask:
        ctrl_color = "dimgray"
        base_color = "C1"
    else:
        ctrl_color = "silver"
        base_color = "dimgray"
    ax.plot(
        tr_control_data[:, 0],
        tr_control_data[:, 1],
        label="Genome-seq",
        color=ctrl_color,
    )

    # Plot treat pileup
    peaks = np.ma.masked_where(
        tr_treat_data[:, 1] < tr_control_data[:, 1], tr_treat_data[:, 1]
    )
    base_line = np.ma.masked_where(
        tr_treat_data[:, 1] > tr_control_data[:, 1], tr_treat_data[:, 1]
    )
    x = tr_treat_data[:, 0]

    ax.plot(x, base_line, label="ChIP-seq", color=base_color, alpha=peak_alpha)
    ax.plot(
        x, peaks, label="ChIP-seq peaks", color=peak_color, alpha=peak_alpha
    )
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

    if ymax_given is not None:
        ymax = ymax_given
    if ymin_given is not None:
        ymin = ymin_given
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


def try_imagemagick_pdf(png_files: list[Path], pdf_path: Path) -> bool:
    """Try to create PDF using ImageMagick. Returns True if successful."""

    try:
        # Use file list to avoid ARG_MAX issues
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        ) as f:
            for png_file in png_files:
                f.write(f"{png_file}\n")
            file_list = f.name

        try:
            cmd = ["magick", f"@{file_list}", "-compress", "lzw", str(pdf_path)]
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            log.info("Successfully created PDF using ImageMagick: %s", pdf_path)
            return True
        finally:
            Path(file_list).unlink(missing_ok=True)

    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        log.warning("ImageMagick failed or not available: %s", e)
        return False


def create_pdf_matplotlib(png_files: list[Path], pdf_path: Path) -> None:
    """Create PDF from PNG files using matplotlib."""

    with PdfPages(str(pdf_path)) as pdf:
        for png_file in png_files:
            try:
                img = Image.open(png_file)
                fig = plt.figure(figsize=(img.width / 100, img.height / 100))
                ax = fig.add_axes((0, 0, 1, 1))
                ax.imshow(img)
                ax.axis("off")
                pdf.savefig(fig, bbox_inches="tight", pad_inches=0)
                plt.close(fig)
                log.info("Added %s to PDF", png_file)
            except (UnidentifiedImageError, OSError, ValueError) as e:
                log.error("Failed to add %s to PDF: %s", png_file, e)
    log.info("Successfully created PDF using matplotlib: %s", pdf_path)


def consolidate_files(
    png_files: list[Path], base_name: str, output_dir: Path
) -> None:
    """Create consolidated PDF and ZIP files, then remove intermediate PNG files."""
    if not png_files:
        log.warning("No PNG files to consolidate")
        return

    # Create PDF
    pdf_path = output_dir / f"{base_name}_plots.pdf"

    # Try ImageMagick first, fall back to matplotlib
    if not try_imagemagick_pdf(png_files, pdf_path):
        log.info("Falling back to matplotlib for PDF creation")
        create_pdf_matplotlib(png_files, pdf_path)

    # Create ZIP file
    zip_path = output_dir / f"{base_name}_plots.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
        for png_file in png_files:
            zipf.write(png_file, png_file.name)
            log.info("Added %s to ZIP", png_file)
    log.info("Successfully created ZIP file: %s", zip_path)

    # Remove intermediate PNG files
    for png_file in png_files:
        try:
            png_file.unlink()
            log.info("Removed intermediate file: %s", png_file)
        except OSError as e:
            log.error("Failed to remove %s: %s", png_file, e)

    log.info("Consolidation complete. Generated: %s and %s", pdf_path, zip_path)


def __main__():
    argparser = arg_parser()
    args = argparser.parse_args()
    if args.savefig.suffix not in [".png", ".jpg", ".jpeg", ".tiff", ".bmp"]:
        # do not use "with_suffix" as it will remove false suffix if
        # any dot in the name
        args.savefig = args.savefig.parent / (args.savefig.name + ".png")

    if not args.peak_list and not args.region and not args.gene:
        log.info("Peak list is not provided, will try to get from macs output.")
        args.peak_list = list(args.macsOutput.glob("*_peaks.xls"))[0]

    try:
        genome_with_annotation = SeqIO.read(args.genome.expanduser(), "genbank")
    except ValueError as e:
        log.warning(
            "Failed to read genome with annotation from %s: %s",
            args.genome,
            e,
        )
        SeqIO.read(args.genome.expanduser(), "fasta")
        genome_with_annotation = None

    # Track generated PNG files for consolidation
    generated_png_files = []

    if not args.peak_list:  # Plot a single region
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
            peak_alpha=args.peak_alpha,
            ymax_given=args.ymax,
            ymin_given=args.ymin,
            no_control_mask=args.no_control_mask,
            do_logscale=args.logscale,
            genome_with_annotation=genome_with_annotation,
        )
        ax.set_title(args.title)
        fig.savefig(args.savefig, dpi=args.dpi)
    else:
        peaks = read_peak_file(args.peak_list)
        for _, peak in peaks.iterrows():
            start, end = peak["start"], peak["end"]
            try:
                location = peak["summit"]
            except KeyError:
                location = (start + end) // 2
            try:
                length = peak["length"]
            except KeyError:
                length = end - start
            # Expand both side by length of "summit to edge plus peak length"
            # or 1500 bp if that value < 1500
            flanking = max(
                1500, location - start - length, end - location - length
            )
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
                ymax_given=args.ymax,
                ymin_given=args.ymin,
                peak_alpha=args.peak_alpha,
                do_logscale=args.logscale,
                no_control_mask=args.no_control_mask,
                genome_with_annotation=genome_with_annotation,
            )
            if "_temp" not in args.savefig.name:
                savefig = args.savefig.with_name(
                    args.savefig.stem + f"_{peak.name}.png"
                )
            else:
                savefig = Path(f"./{peak.name}.png")
            title = f"ChIP-seq pileup comparison - {peak.name}"
            if args.title:
                title = f"{args.title} - {peak.name}"
            log.info("Saving plot to %s", savefig)
            ax.set_title(title)
            fig.savefig(savefig, dpi=args.dpi)
            plt.close(fig)

            # Track the generated PNG file
            generated_png_files.append(savefig)

    # Handle consolidation if requested
    if args.consolidate and generated_png_files:
        base_name = (
            args.savefig.stem if args.savefig.stem != "__temp" else "peakplots"
        )
        output_dir = Path(args.savefig).parent
        consolidate_files(generated_png_files, base_name, output_dir)


if __name__ == "__main__":
    __main__()
