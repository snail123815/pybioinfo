# This file is licensed under the MIT License

############################
# by Chao DU (杜超)
# c.du@biology.leidenuniv.nl
# durand.dc@gmail.com
############################

import argparse
import os
from collections import OrderedDict

from BCBio import GFF
from Bio import SeqIO

parser = argparse.ArgumentParser(
    usage="Add info from GFF file to sequence file (.fa, .gbk)"
)
parser.add_argument(
    "target", help="Sequence file which the annotation info will be added"
)
parser.add_argument(
    "gffs", nargs="+", help="GFF format file contain annotation."
)

args = parser.parse_args()
target = args.target
gffs = args.gffs

# need the id of all sequences in the target sequence, so load them into memory
try:
    seqs = list(SeqIO.parse(target, "genbank"))
except Exception as e:
    print("Genbank read fail.")
    print(type(e), e)
if len(seqs) == 0:
    try:
        seqs = list(SeqIO.parse(target, "fasta"))
    except Exception as e:
        print("fasta read fail.")
        print(type(e), e)
assert len(seqs) > 0, "There is no sequence in target file" + target

seqDict = OrderedDict()
for seq in seqs:
    seqDict[seq.id] = seq
    if "molecule_type" not in seq.annotations:
        seqDict[seq.id].annotations["molecule_type"] = "DNA"

gffDict = OrderedDict()  # values are list
for gff in gffs:
    with open(gff, "r") as handle:
        grecs = list(GFF.parse(handle))
        for rec in grecs:
            if rec.id not in gffDict:
                gffDict[rec.id] = [rec]
            else:
                gffDict[rec.id].append(rec)


changed = 0


def addGffAnn(seqRec, gffRecs):
    global changed
    for grec in gffRecs:
        seqRec.features.extend(grec.features)
    print(f"Info added to sequence {seqRec.id}")
    changed += 1
    return seqRec


listTar = list(seqDict.keys())
listGff = list(gffDict.keys())
comKeys = sorted(list(set(listTar).intersection(set(listGff))))
difKeys = sorted(list(set(listTar).symmetric_difference(set(listGff))))
diffTar = sorted(
    list(set(listTar).difference(set(listGff)))
)  # in tar not in gff
diffGff = sorted(
    list(set(listGff).difference(set(listTar)))
)  # in gff not in tar

for k in comKeys:
    seqDict[k] = addGffAnn(seqDict[k], gffDict[k])
if len(difKeys) > 0:
    if len(seqDict) == 1:
        keyTar = listTar[0]
        for keyGff in diffGff:
            force = ""
            while force not in ["y", "n"]:
                force = input(
                    f"Source id ({keyTar}) != gff id ({keyGff}), force merge? (y/n)"
                )
            if force == "y":
                seqDict[keyTar] = addGffAnn(seqDict[keyTar], gffDict[keyGff])
            else:
                print("Ignored.")
    else:
        print(
            f"Unmatched sequence ID target ({diffTar}) and gff ({diffGff}) ignored."
        )


if changed > 0:
    out = os.path.splitext(target)[0] + ".appended.gbk"
    print(f"Write to {out}.")
    SeqIO.write(seqDict.values(), out, "genbank")
