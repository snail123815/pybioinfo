# This file is licensed under the MIT License

import argparse
from pathlib import Path

from pyBioinfo_modules.bio_sequences.features_from_gbk import getCdsFromGbk

argparser = argparse.ArgumentParser()
argparser.add_argument("file", help="genbank file")

args = argparser.parse_args()
gbkPath = Path(args.file)
faaPath = getCdsFromGbk(gbkPath)
print(faaPath)
