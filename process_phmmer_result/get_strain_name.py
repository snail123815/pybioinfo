import json
import os
import re

import numpy as np
import termplotlib as tpl
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

jsonFile = "./C7B6EDEE-0EA7-11ED-97D6-B0F3E4368550.1.json"
relationsFile = "./C7B6EDEE-0EA7-11ED-97D6-B0F3E4368550.1.relations.tsv"
treeFile = "./tree.20220729.tre"

outputTreeFile = ".species".join(os.path.splitext(treeFile))

with open(jsonFile, "r") as jf:
    jinfo = json.load(jf)
relations: list[tuple[str, str]] = []
uniqueSpecies = set()


def safeName(name: str) -> str:
    return re.sub(r"[ _:,();{}+*'\"[\]\/\t\n]+", "_", name)


for hit in jinfo["results"]["hits"]:
    acc = hit["acc"]
    species = safeName(hit["species"])
    if "_strain_" in species:
        species = "_".join(species.split("_strain_"))
    if species.endswith("_"):
        species = species[:-1]
    spsp = species.split("_")
    spstrain = spsp[2:]
    if len(spstrain) > 4:
        spstrain = spstrain[:2] + spstrain[-2:]
    species = "_".join(spsp[:2] + spstrain)
    suffixNum = 0
    while species in uniqueSpecies:
        if suffixNum > 0:
            species = species[: -len(str(suffixNum)) - 1]
        suffixNum += 1
        species += "_" + str(suffixNum)
    uniqueSpecies.add(species)
    relations.append((acc, species))

with open(relationsFile, "w") as rf:
    for acc, species in relations:
        rf.write(acc + "\t" + species + "\n")

with open(treeFile, "r") as inTree:
    with open(outputTreeFile, "w") as outTree:
        for l in inTree:
            newl = l
            for acc, species in relations:
                nodeName = re.compile(acc + r"\/\d{1,3}-\d{1,3}:")
                newl = nodeName.sub(acc + "_" + species + ":", newl)
            outTree.write(newl)
