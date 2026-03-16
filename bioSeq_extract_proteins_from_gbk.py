# This file is licensed under the MIT License

import argparse
from pathlib import Path

from pyBioinfo_modules.bio_sequences.features_from_gbk import getFaaFromGbk

argparser = argparse.ArgumentParser()
argparser.add_argument("file", help="genbank file")
argparser.add_argument(
    "-s",
    "--seq_id_from",
    help="feature name that will be used as sequence id",
    default="protein_id",
)
argparser.add_argument(
    "--prefix",
    help="prefix to add to each extracted sequence ID",
    default=None,
)
argparser.add_argument(
    "--suffix",
    help="suffix to add to each extracted sequence ID",
    default=None,
)

args = argparser.parse_args()
gbkPath = Path(args.file)
faaPath = getFaaFromGbk(
    gbkPath,
    getIdFrom=args.seq_id_from,
    prefix=args.prefix,
    suffix=args.suffix,
)
print(faaPath)
