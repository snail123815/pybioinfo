import argparse
import json
from pathlib import Path

from ..pyBioinfo_modules.basic.basic import safe_name
from ..pyBioinfo_modules.basic.tree import change_tree_node_name


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            """
            Change the name of proteins in the tree to the name of
            specieces. The information is fetched from json output of phmmer
            web application.
            Will also write a tsv file with two columns:
            protein accession and species name.
            Duplicated species names will be suffixed with a number, ranked by
            the order of appearance in the json file (E values of hits, not
            domains).
            """
        )
    )
    parser.add_argument(
        "-j",
        "--json",
        type=Path,
        required=True,
        help="Input JSON file. Fetched from phmmer web application.",
    )
    parser.add_argument(
        "tree",
        type=Path,
        required=True,
        help=(
            "Input tree file. Has to be generated using the proteins from "
            "the JSON file. Newick format. No duplicated node names in the tree."
        ),
    )
    return parser.parse_args()


def parse_species_of_proteins_hmmerjson(json_path: Path):

    with open(json_path, "r") as jf:
        jinfo = json.load(jf)

    relations = {}
    unique_species = set()

    for hit in jinfo["results"]["hits"]:
        acc = hit["acc"]
        species = safe_name(hit["species"])
        if "_strain_" in species:
            species = "_".join(species.split("_strain_"))
        if species.endswith("_"):
            species = species[:-1]
        sp_split = species.split("_")
        strain_split = sp_split[2:]
        if len(strain_split) > 4:
            # If the strain is recorded like this, there mush be duplicated
            # values.
            strain_split = strain_split[:2] + strain_split[-2:]
        # Add strain to the full species name
        species = "_".join(sp_split[:2] + strain_split)
        suffixNum = 0
        while species in unique_species:
            # Add a number to the species name if duplicated, these might
            # indicate paralogues.
            if suffixNum > 0:
                species = species[: -len(str(suffixNum)) - 1]
            suffixNum += 1
            species += "_" + str(suffixNum)
        unique_species.add(species)
        relations[acc] = species
    return relations, unique_species


if __name__ == "__main__":
    args = parse_args()
    json_path = args.json
    tree_path = args.tree

    relations_table_path = json_path.with_name(
        json_path.stem + ".relations" + json_path.suffix
    )
    species_tree_path = tree_path.with_name(
        tree_path.stem + ".species" + tree_path.suffix
    )

    relations, unique_species = parse_species_of_proteins_hmmerjson(json_path)

    with open(relations_table_path, "w") as rf:
        for (
            acc,
            species,
        ) in relations.items():  # Changed to iterate over dictionary items
            rf.write(acc + "\t" + species + "\n")

    # Use the new change_tree_node_name function
    change_tree_node_name(tree_path, relations)
