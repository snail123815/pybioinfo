# This file is licensed under the MIT License

import argparse
import shutil
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
argparser.add_argument(
    "-o",
    "--output-dir",
    help="output directory for extracted sequences",
    default=None,
)

args = argparser.parse_args()
gbkPath = Path(args.file)
output_dir = Path(args.output_dir) if args.output_dir else None
if output_dir is not None:
    output_dir.mkdir(parents=True, exist_ok=True)

faaPath = getFaaFromGbk(
    gbkPath,
    getIdFrom=args.seq_id_from,
    prefix=args.prefix,
    suffix=args.suffix,
)
if output_dir is not None:
    target_path = output_dir / faaPath.name
    if faaPath.resolve() != target_path.resolve():
        shutil.move(str(faaPath), str(target_path))
    faaPath = target_path

print(faaPath)
