import re
from copy import deepcopy

from Bio import Phylo


def change_tree_node_name(
    origional_tree: Phylo.Newick.Tree,
    relations: dict[str, str],
    name_regex: str = r"",
):
    """
    In one line
    ...,((((A0A0L0KMM2_9ACTN/6-172:0.078341483,
    ((((A0A4Q7XLU4_9ACTN/2-380:0.061149438,
    A0A1V6MSQ2_9ACTN/10-366:0.058716814)0.736:0.011935775,
    ((A0A1H9F5N7_9ACTN/13-368:0.017549879,
    A0A0Q9BBU6_9ACTN/12-365:0.057577268)0.324:0.003375517,...

    Update terminal node names in a phylogenetic tree using BioPython.

    Parameters:
        origional_tree (Bio.Phylo.Newick.Tree): Tree object.
        relations (dict[str, str]): Mapping from accession IDs to species names.
        name_regex (str): Regular expression to match the node names.
            Needs a match group for the accession ID.
            e.g. r"(.+?)\/\d{1,3}-\d{1,3}"

    Returns:
        Bio.Phylo.Newick.Tree: A new tree instance with updated node labels.
    """

    tree = deepcopy(origional_tree)
    # Iterate through all terminal nodes (leaves)
    for terminal in tree.get_terminals():
        # Extract the accession from the node name using regex
        if name_regex:
            match = re.match(name_regex, terminal.name)
            if match:
                acc = match.group(1)
                if acc in relations:
                    # Update the node name with species information
                    terminal.name = f"{acc}_{relations[acc]}"
        else:
            if terminal.name in relations:
                # Update the node name with species information
                terminal.name = f"{terminal.name}_{relations[terminal.name]}"

    return tree
