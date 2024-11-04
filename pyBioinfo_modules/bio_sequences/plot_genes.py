import matplotlib.transforms as mtransforms
from matplotlib.patches import FancyArrow, Rectangle
from matplotlib.axes import Axes
from Bio.SeqRecord import SeqRecord
from typing import List, Literal
from matplotlib.lines import Line2D
from matplotlib.text import Text


def make_gene_appear_truncated(
    ax: Axes,
    truncate_loc: float,
    y_loc: float,
    direction: Literal[-1, 1],
    block_width: float,
    block_length: float,
    num_segments=2,
    transform=None,
    zorder=None,
) -> List[Rectangle]:
    blocks = []
    for i in range(num_segments):
        segment_start = truncate_loc + direction * block_length * (
            2 * i + (2 if direction == -1 else 1)
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
) -> tuple[List[FancyArrow], dict[str:Rectangle], Line2D, list[Text]]:
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

            if partial_left:
                truncate_blocks["left"] = make_gene_appear_truncated(
                    ax=ax,
                    truncate_loc=gene_start,
                    y_loc=arrow_y_loc,
                    direction=1,
                    block_width=head_width * 1.02,
                    # A multiplier > 1 to make the block just a bit wider than
                    # the arrow head, preventing line leakage when rendered.
                    block_length=head_length * 0.4,
                    transform=trans,
                    zorder=1,
                )
                if gene_strand == -1:
                    # Reduce the genome_line from left by an arrow head length
                    genome_line.set_xdata(
                        [region_start + head_length, region_end]
                    )
            if partial_right:
                truncate_blocks["right"] = make_gene_appear_truncated(
                    ax=ax,
                    truncate_loc=gene_end,
                    y_loc=arrow_y_loc,
                    direction=-1,
                    block_width=head_width * 1.02,
                    block_length=head_length * 0.4,
                    transform=trans,
                    zorder=1,
                )
                if gene_strand == 1:
                    # Reduce the genome_line from right by an arrow head length
                    genome_line.set_xdata(
                        [region_start, region_end - head_length]
                    )
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

    # Reduce the length of the genome line by 3 times the pixel length
    # to prevent the line from touching the end of the plotted gene at render,
    # especially when the gene is at the edge of the plot.
    pixel_length = (region_end - region_start) / ax.get_window_extent().width
    genome_line.set_xdata(
        [
            genome_line.get_xdata()[0] + 3 * pixel_length,
            genome_line.get_xdata()[1] - 3 * pixel_length,
        ]
    )

    if len(gene_names) >= max_genes:
        interval = len(gene_names) // 10
        for i, gene_name in enumerate(gene_names):
            if i % interval != 0:
                gene_name.set_visible(False)

    return (
        arrows,
        truncate_blocks,
        genome_line,
        gene_names,
    )
