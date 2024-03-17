import argparse
from pathlib import Path


from pyBioinfo_modules.bio_sequences.features_from_gbk \
    import getFaaFromGbk


argparser = argparse.ArgumentParser()
argparser.add_argument('file', help='genbank file')

args = argparser.parse_args()
gbkPath = Path(args.file)
faaPath = getFaaFromGbk(gbkPath)
print(faaPath)
