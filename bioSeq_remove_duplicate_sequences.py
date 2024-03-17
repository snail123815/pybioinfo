import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


parser = argparse.ArgumentParser()
parser.add_argument('-f', type=Path, help='.fasta file to de-duplicate')

args = parser.parse_args()

uniqSeqRecs: list[SeqRecord] = []

totalRecN = 0
seqgroupN = 0
seqGroupPrefix = 'SeqGroup'
seqGroups: dict[str, tuple[Seq, list[str]]] = {}
for rec in SeqIO.parse(args.f, 'fasta'):
    totalRecN += 1
    isUnique = True
    for i, t in enumerate(uniqSeqRecs):
        if t.seq == rec.seq:
            isUnique = False
            if t.id.startswith(seqGroupPrefix):
                seqGroups[t.id][1].append(rec.description)
            else:
                seqgroupN += 1
                newId = f'{seqGroupPrefix}{seqgroupN:0>4}'
                newDesc = (
                    f'{t.description}|{rec.description}'
                )
                uniqSeqRecs[i] = SeqRecord(t.seq, id=newId, description="")
                seqGroups[newId] = (t.seq, [t.description, rec.description])
            break
    if isUnique:
        uniqSeqRecs.append(rec)

assert len(set([r.seq for r in uniqSeqRecs])) == len(uniqSeqRecs)

print(
    f'Parsed {totalRecN} sequences, found {len(uniqSeqRecs)}'
    f' unique sequences.'
)

SeqIO.write(uniqSeqRecs, args.f.with_suffix(
    f'.unique{args.f.suffix}'), 'fasta')
with open(args.f.with_suffix('.non_unique.tsv'), 'w') as handle:
    handle.write('SeqGroup\tn\tSequence\tSeq descriptions\n')
    for id, [s, descs] in seqGroups.items():
        handle.write(
            f'{id}\t{len(descs)}\t{s}\t{"|".join(descs)}\n'
        )


print("DONE")
