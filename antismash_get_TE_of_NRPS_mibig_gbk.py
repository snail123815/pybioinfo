from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from pyBioinfo_modules.basic.decompress import getSuffixIfCompressed
from pyBioinfo_modules.basic.decompress import decompFileIfCompressed
from tqdm import tqdm

parser = argparse.ArgumentParser()

parser.add_argument('-p', type=Path, help='Path to gbk(s)')
parser.add_argument('-o', type=Path, help='Path to output fasta file')
parser.add_argument('-t', '--outputType', type=str, default="nt",
                    choices=['nt', 'aa'],
                    help='Type of output (aa/nt)')

args = parser.parse_args()

if args.p.is_file():
    assert getSuffixIfCompressed(args.p) in ['.gbk', '.gbff']
    gbks = [decompFileIfCompressed(args.p)[0]]
elif args.p.is_dir():
    gbks = []
    for f in args.p.iterdir():
        if f.is_dir():
            continue
        if not getSuffixIfCompressed(f) in ['.gbk', '.gbff']:
            continue
        gbks.append(decompFileIfCompressed(f)[0])
    assert len(gbks) > 0
else:
    raise ValueError('Do not find valid files in(as) the given path.')


def TE_follow_PCP(aSDomains: list[dict]) -> list[dict]:
    teIndexs = []
    validTes = []
    try:
        for i, dom in enumerate(aSDomains):
            if dom['name'] == 'Thioesterase':
                teIndexs.append(i)
        for i in teIndexs:
            try:
                if aSDomains[i - 1]['name'] == 'PCP':
                    validTes.append(aSDomains[i])
                    break
            except IndexError:
                pass
            try:
                if aSDomains[i + 1]['name'] == 'PCP':
                    validTes.append(aSDomains[i])
            except IndexError:
                pass
    except ValueError:
        pass
    return validTes


validTes = {}
for gbk in tqdm(gbks):
    try:
        s = SeqIO.read(gbk, 'genbank')
    except ValueError as ve:
        if str(ve) == "More than one record found in handle":
            print(f'There are more than one record in gbk file {gbk}'
                  ' This is not an AntiSMASH result record.')
        else:
            raise ve
    aSDomains = []
    for feat in s.features:
        if feat.type != 'aSDomain':
            continue
        aSDomains.append({
            'name': feat.qualifiers['aSDomain'][0],
            'location': feat.location,
            'locus_tag': feat.qualifiers['locus_tag'][0],
            'protein_start': feat.qualifiers['protein_start'][0],
            'protein_end': feat.qualifiers['protein_end'][0],
            'protein_id': "",
            'translation': Seq(feat.qualifiers['translation'][0]),
            'seq': s[feat.location.start:feat.location.end + 1].seq
        })
    # fill in protein_id

    def all_indices(l, t):
        indices = []
        for i, x in enumerate(l):
            if x == t:
                indices.append(i)
        return indices
    locus_tags = [d['locus_tag'] for d in aSDomains]
    for feat in s.features:
        if feat.type != 'CDS':
            continue
        toTry = ['locus_tag', 'gene', 'product', 'protein_id']
        for t in toTry:
            if t in feat.qualifiers:
                try:
                    indices = all_indices(locus_tags, feat.qualifiers[t][0])
                except ValueError:
                    continue
                if len(indices) > 0:
                    break
        if indices == 0:
            continue
        for i in indices:
            aSDomains[i]['protein_id'] = feat.qualifiers['protein_id'][0]
    for i in range(len(aSDomains)):
        if aSDomains[i]['protein_id'] == '':
            aSDomains[i]['protein_id'] = aSDomains[i]['locus_tag']

    validTes[gbk] = TE_follow_PCP(aSDomains)

validTeSeqs = []

for gbk in validTes:
    for i, te in enumerate(validTes[gbk]):
        if len(validTes[gbk]) == 1:
            id = f'{gbk.stem}_TE'
        else:
            id = f'{gbk.stem}_TE_{i}'
        if args.outputType == 'nt':
            s = te['seq']
        else:
            s = te['translation']
        validTeSeqs.append(SeqRecord(s, id=id, description=f"{id}_{gbk.name}"))

SeqIO.write(validTeSeqs, args.o, 'fasta')

print("Done")
