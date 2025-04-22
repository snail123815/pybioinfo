# This tool is for single use, process json file from HMMER web
# It is not thoroughly considered for other cases
# TODO make it work with domtbl output from HMMER, wrapper.hmmer.filter...

import json
import os

import numpy as np
import termplotlib as tpl
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

jsonFile = "./C7B6EDEE-0EA7-11ED-97D6-B0F3E4368550.1.json"
alignmentFasta = "./C7B6EDEE-0EA7-11ED-97D6-B0F3E4368550.1.afa"

tStart = 110
tEnd = 310
eThresh = 1e-50

targetHitsAlignmentFasta = ".targets".join(os.path.splitext(alignmentFasta))

with open(jsonFile, "r") as jf:
    jinfo = json.load(jf)
relations: list[tuple[str, str]] = []
uniqueSpecies = set()
allAlignments = SeqIO.parse(alignmentFasta, "fasta")

hits: list[tuple[str, str]] = []
for hit in jinfo["results"]["hits"]:
    isCompleteHit = False
    for dom in hit["domains"]:
        qs = dom["alihmmfrom"]  # query start
        qe = dom["alihmmto"]  # query end
        if qs < tStart and qe > tEnd:
            isCompleteHit = True
            break
    if isCompleteHit:
        hits.append((hit["name"], hit))
print(f"Found {len(hits)} proteins cover the target region {tStart}-{tEnd}")

evalues = []
for acc, hit in hits:
    evalues.append(float(hit["evalue"]))
print("\nlog(evalue) distribution of found proteins:")
counts, bin_edges = np.histogram(np.log(evalues), bins=30)
fig = tpl.figure()
fig.hist(counts, bin_edges, orientation="horizontal", force_ascii=False)
fig.show()

hitAccs = [h[0] for h in hits]
targetSeqs: list[tuple[str, int, SeqRecord]] = []
for i, seq in enumerate(allAlignments):
    if i == 0:
        continue
    acc, location = seq.id.split("/")
    s, e = (int(i) for i in location.split("-"))
    l = e - s + 1
    if acc in hitAccs:
        targetSeqs.append((acc, l, seq))

hitSeqs = []
for acc in hitAccs:
    maxl = 0
    for seqAcc, l, seq in targetSeqs:
        if seqAcc == acc:
            if l > maxl:
                maxl = l
                hitSeq = seq
    if maxl > 0:
        hitSeqs.append(hitSeq)
assert len(hitSeqs) == len(hits), f"{len(hitSeqs)}, {len(hits)}"
SeqIO.write(hitSeqs, targetHitsAlignmentFasta, "fasta")
