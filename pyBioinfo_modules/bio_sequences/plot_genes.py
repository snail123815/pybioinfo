import logging
from typing import List, Literal

import matplotlib.transforms as mtransforms
from Bio.SeqRecord import SeqRecord
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrow, Rectangle
from matplotlib.text import Text

logger = logging.getLogger(__name__)


def make_gene_appear_truncated(
    ax: Axes,
    truncate_loc: float,
    y_loc: float,
    block_direction: Literal[-1, 1],
    block_width: float,
    block_length: float,
    num_segments=2,
    transform=None,
    zorder=None,
) -> List[Rectangle]:
    """
    Create a visual representation of a truncated gene on a plot by adding
    rectangular blocks at the truncated end.

    Parameters:
    ax (Axes): The matplotlib Axes object to draw the blocks on.
    truncate_loc (float): The location on the x-axis where the truncation starts.
    y_loc (float): The location on the y-axis where the blocks will be centered.
    direction (Literal[-1, 1]): The direction in which the blocks will be drawn.
                                Use -1 for left and 1 for right.
    block_width (float): The width of each block.
    block_length (float): The length of each block.
    num_segments (int, optional): The number of blocks to draw. Default is 2.
    transform (optional): The transformation applied to the blocks. Default is None.
    zorder (optional): The z-order of the blocks. Default is None.

    Returns:
    List[Rectangle]: A list of matplotlib Rectangle objects representing
                     the blocks.
    """
    blocks = []
    for i in range(num_segments):
        segment_start = truncate_loc + block_direction * block_length * (
            2 * i + (2 if block_direction == -1 else 1)
        )
        block = Rectangle(
            (segment_start, y_loc - block_width / 2),
            block_length,
            block_width,
            linewidth=0,
            edgecolor=None,
            facecolor=ax.get_facecolor(),
            transform=transform,
            zorder=zorder,
        )
        blocks.append(block)
        ax.add_patch(block)
    return blocks


# Function to plot genes as arrows
def plot_genes(
    ax: Axes,
    genome_with_annotation: SeqRecord,
    region_start: int,
    region_end: int,
    color: List[str] = ["C0", "C2"],
    max_genes: int = 10,
) -> tuple[List[FancyArrow], dict[str:[Rectangle]], Line2D, list[Text]]:
    """
    Plots genes on a given matplotlib Axes object within a specified
    genomic region.

    Parameters:
    ax (Axes): Matplotlib Axes object where the genes will be plotted.
    genome_with_annotation (SeqRecord): A SeqRecord object containing the genome
                                        sequence and its annotations.
    region_start (int): The start position of the genomic region to be plotted.
    region_end (int): The end position of the genomic region to be plotted.
    color (List[str], optional): List of colors for the genes on forward and
                                 reverse strand, respectively.
                                 Default is ["C0", "C2"].
    max_genes (int, optional): Maximum number of gene names to display.
                               Default is 10.

    Returns:
    tuple: A tuple containing:
        - List[FancyArrow]: List of FancyArrow objects representing the genes.
        - dict[str:list[Rectangle]]: Dictionary with keys "left" and "right"
                               containing Rectangle objects for truncated genes.
        - Line2D: Line2D object representing the genome sequence blackline.
        - list[Text]: List of Text objects (axes.Text) for gene names. Including
                      the gene names that are not displayed if max_genes is
                      exceeded.
    """
    xrange = ax.get_xlim()[1] - ax.get_xlim()[0]
    arrow_y_loc = 0.1
    # Gene x-axis position will be defined by gene location itself,
    # while y-axis position will be fixed to axes coordinate
    trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    # Plot line representing the genome sequence
    genome_line = ax.plot(
        [region_start, region_end],
        [arrow_y_loc, arrow_y_loc],
        color="black",
        linewidth=1.5,
        transform=trans,
        zorder=0,
    )[0]
    gene_names = []
    arrows = []
    truncate_blocks = {"left": [], "right": []}

    gene_nr_in_region = 0
    center_gene_nr = None
    for feature in genome_with_annotation.features:
        if feature.type == "gene":
            gene_start = feature.location.start
            gene_end = feature.location.end
            gene_strand = feature.location.strand

            if gene_end < region_start:
                continue
            if gene_start > region_end:
                break

            if gene_start < (region_start / 2 + region_end / 2) < gene_end:
                center_gene_nr = gene_nr_in_region
            gene_nr_in_region += 1
            partial = None
            if gene_start < region_start:
                gene_start = region_start
                partial = "left"
            if gene_end > region_end - 1:
                gene_end = region_end - 1
                partial = "right"

            arrow_color = color[0] if gene_strand == 1 else color[1]

            # Define fixed head width and head length
            head_width = 0.05
            # Head length is 0.7% of the x-axis range, reduce if gene is shorter
            head_length = min(abs(gene_end - gene_start), 0.007 * xrange)

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
                zorder=1,
            )
            arrows.append(arrow)
            ax.add_patch(arrow)

            # Add truncated blocks
            if partial:
                trunc_loc = gene_start if partial == "left" else gene_end
                block_direction = 1 if partial == "left" else -1
                truncate_blocks[partial] = make_gene_appear_truncated(
                    ax=ax,
                    truncate_loc=trunc_loc,
                    y_loc=arrow_y_loc,
                    block_direction=block_direction,
                    block_width=head_width * 1.02,
                    # A multiplier > 1 to make the block just a bit wider than
                    # the arrow head, preventing line leakage when rendered.
                    block_length=head_length * 0.4,
                    transform=trans,
                    zorder=1,
                )
                # Reduce the genome_line from left by an arrow head length
                if gene_strand == -1:
                    genome_line.set_xdata(
                        [region_start + head_length, region_end]
                    )
                else:
                    genome_line.set_xdata(
                        [region_start, region_end - head_length]
                    )

            # Add gene name text
            gene_names.append(
                ax.text(
                    (gene_start + gene_end) / 2,
                    arrow_y_loc - 0.05,
                    feature.qualifiers.get("gene", [""])[0],
                    ha="center",
                    va="center_baseline",
                    fontsize=8,
                    transform=trans,
                )
            )

    # Reduce the length of the genome line by 3 times the pixel length to
    # prevent the line from touching the end of the plotted gene at render,
    # especially when the gene is at the edge of the plot. Leakage is bad.
    pixel_length = (region_end - region_start) / ax.get_window_extent().width
    genome_line.set_xdata(
        [
            genome_line.get_xdata()[0] + 3 * pixel_length,
            genome_line.get_xdata()[1] - 3 * pixel_length,
        ]
    )

    # Hide gene names if there are too many,
    # only show 10 genes at an average interval
    if len(gene_names) >= max_genes:
        interval = len(gene_names) // max_genes
        if center_gene_nr:
            for i, gene_name in enumerate(gene_names[center_gene_nr:]):
                if i % interval != 0:
                    gene_name.set_visible(False)
            for i, gene_name in enumerate(
                reversed(gene_names[:center_gene_nr])
            ):
                if (interval + i + 1) % interval != 0:
                    gene_name.set_visible(False)
        else:
            for i, gene_name in enumerate(gene_names):
                if (interval // 2 + i + 1) % interval != 0:
                    gene_name.set_visible(False)

    return (
        arrows,
        truncate_blocks,
        genome_line,
        gene_names,
    )
