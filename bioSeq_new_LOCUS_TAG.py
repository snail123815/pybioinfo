# From assembly annotation output of Quasan.py 25f21f58947c42b842bd228fd7252e49a68aa880
# https://gitlab.services.universiteitleiden.nl/ibl-bioinformatic/streptidy/-/blob/master/Quasan.py
# centre = 'MBT'
# species = 'sp.'
# ndate = '[datetime.datetime.now()]'
# version = ndate.strftime("V%d.%m.%y")
# tag = '[os.path.basename(args.indir)]'
# assembly_version = version + "_" + "hybrid_flye-pilon_" + tag
# assembly_version = version + "_" + "illumina_shovill_" + tag
# assembly_version = version + "_" + "pacbio_flye_" + tag
# prefix = assembly_version + "_prokka"
# assemblies = "[glob.glob(assembly_dir+'/*.f*a')]"
# assembly = "[assemblies[0]]"
# tag = 'strain_ID'
# cmd_prokka = f"prokka --compliant --centre {centre} --genus [args.genus] --species {species}"
# f" --strain {tag} --outdir [workdir] --prefix {prefix} --gcode 11 --cpu [args.threads]"
# f" --locustag \"{tag}_LOCUS_TAG\" --addgenes --gram [args.gram] --rfam --force {assembly}"
# ---------------------------------
# Prokka instruction on --compliant:
# Register your BioProject (e.g. PRJNA123456) and your locus_tag prefix (e.g. EHEC) first!
#% prokka --compliant --centre UoN --outdir PRJNA123456 --locustag EHEC --prefix EHEC-Chr1 contigs.fa
## Check to see if anything went really wrong
#% less PRJNA123456/EHEC-Chr1.err
## Add final details using Sequin
#% sequin PRJNA123456/EHEC-Chr1.sqn
# ---------------------------------
# output dir tree:
#├── V03.12.21_pacbio_flye_PG2_prokka.err
#├── V03.12.21_pacbio_flye_PG2_prokka.faa
#├── V03.12.21_pacbio_flye_PG2_prokka.ffn
#├── V03.12.21_pacbio_flye_PG2_prokka.fna
#├── V03.12.21_pacbio_flye_PG2_prokka.fsa
#├── V03.12.21_pacbio_flye_PG2_prokka.gbk
#├── V03.12.21_pacbio_flye_PG2_prokka.gff
#├── V03.12.21_pacbio_flye_PG2_prokka.log
#├── V03.12.21_pacbio_flye_PG2_prokka.sqn
#├── V03.12.21_pacbio_flye_PG2_prokka.tbl
#├── V03.12.21_pacbio_flye_PG2_prokka.tsv
#└── V03.12.21_pacbio_flye_PG2_prokka.txt
# There are a lot of errors show in the .err file.
# Fatal error is BAD_LOCUS_TAG_FORMAT, this is exactly what I will change.
# ---------------------------------

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('p', help='Path to annotation folder')
parser.add_argument('locustag_old', help='Locus tag to replace')
parser.add_argument('locustag_new', help='Correct locus tag')

args = parser.parse_args()
annoPath = os.path.abspath(args.p)
oldTag = args.locustag_old
newTag = args.locustag_new
newPath = os.path.join(os.path.dirname(annoPath), os.path.basename(annoPath) + '_' + newTag)

assert os.path.isdir(annoPath)
os.makedirs(newPath,exist_ok=True)

def replaceTag(sf, oldTag, nf, newTag):
    with open(sf, 'r') as s:
        with open(nf, 'w') as t:
            for sl in s:
                t.write(sl.replace(oldTag, newTag))

for sf in os.listdir(annoPath):
    fn, ext = os.path.splitext(sf)
    ext = ext.lower()
    if ext in ['.err', '.txt', '.log']:
        continue
    if ext in ['.faa', '.ffn', '.fna', '.fsa', '.gbk', '.gff', '.sqn', '.tbl', '.tsv']:
        sf = os.path.join(annoPath, sf)
        nf = fn + "_" + newTag
        nf = os.path.join(newPath, nf) + ext
        replaceTag(sf, oldTag, nf, newTag)

