# This file is licensed under the MIT License

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from pyBioinfo_modules.wrappers.antismash import find_NRPS_TE_domain

parser = argparse.ArgumentParser()
parser.add_argument("-p", type=Path, help="Root path of all AntiSMASH folders")
args = parser.parse_args()

rootAllAntismashes = args.p

allAntismashes = [d for d in rootAllAntismashes.iterdir() if d.is_dir()]
asResultJsons = []
for d in allAntismashes:
    jsons = list(d.glob("*.json"))
    assert len(jsons) == 1, (
        f'AntiSMASH result in "{d.name}" has some problem:'
        f'    Found {len(jsons)} of ".json" files, but there should be one.'
    )
    asResultJsons.append(jsons[0])

teDomains: list[SeqRecord] = []

for asResultJson in tqdm(asResultJsons):
    teDomains.extend(find_NRPS_TE_domain(asResultJson))

SeqIO.write(
    teDomains,
    (
        rootAllAntismashes.parent
        / (rootAllAntismashes.name + ".NRPS_TE_domains.fasta")
    ).open("w"),
    "fasta",
)
