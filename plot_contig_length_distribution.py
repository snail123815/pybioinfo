import argparse
import os
from Bio import SeqIO
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(usage=('Generate a plot showing the distribution of contig lengths of all sequences in a fasta file'))
parser.add_argument('fasta', help='path to fasta file')

args = parser.parse_args()
fa = args.fasta

title = os.path.basename(fa)
outputSvg = os.path.splitext(fa)[0]+'.svg'

seqs = list(SeqIO.parse(fa,'fasta'))
seqs.sort(key=lambda x: len(x), reverse=True)
seqIds = [s.id for s in seqs]

fig, ax = plt.subplots(1,1)
ax.plot([len(s) for s in seqs])
ax.set_title(title)

cumax = ax.twinx()
cumList = []
cum = 0
for l in [len(s) for s in seqs]:
    cum += l
    cumList.append(cum)
cumax.plot(cumList)

n50 = [l for l in cumList if l >= cumList[-1]/2][0]
n50i = cumList.index(n50)
n50seq = seqs[n50i]

cumax.scatter(n50i, n50)
cumax.text(n50i, n50, f' N50=\n {n50seq.id}:{len(n50seq):,}', va='center', ha='left')

ax.text(0.98, 0.92, f'Total {cumList[-1]:,}',
        transform=ax.transAxes,
        va='top', ha='right')

i = 0
while i < 15:
    try:
        t = f"{seqIds[i]}: {len(seqs[i]):.2e}"
        ax.text(0.98, 0.82-0.05*i, t,
                transform=ax.transAxes,
                va='top', ha='right')
    except:
        break
    i += 1

fig.savefig(outputSvg)
plt.show()


