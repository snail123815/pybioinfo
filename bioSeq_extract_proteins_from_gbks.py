# This file is licensed under the MIT License

import argparse
import shutil
from pathlib import Path

from pyBioinfo_modules.bio_sequences.features_from_gbk import getFaaFromGbk

argparser = argparse.ArgumentParser()
argparser.add_argument("files", nargs="+", help="genbank file(s)")
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
    help=(
        "suffix to add to each extracted sequence ID. "
        "Defaults to the input file stem (e.g. '_genome' for genome.gbk)"
    ),
    default=None,
)
argparser.add_argument(
    "-j",
    "--join",
    action="store_true",
    help="write all extracted sequences into a single joint output file",
)
argparser.add_argument(
    "--join_name",
    help=(
        "name of the joint output file. "
        "Defaults to the first input file stem + '_etc.proteins.faa'"
    ),
    default=None,
)
argparser.add_argument(
    "-o",
    "--output-dir",
    help="output directory for extracted sequences",
    default=None,
)

args = argparser.parse_args()
output_dir = Path(args.output_dir) if args.output_dir else None
if output_dir is not None:
    output_dir.mkdir(parents=True, exist_ok=True)

if args.join:
    first_stem = Path(args.files[0]).stem
    join_path = Path(
        args.join_name if args.join_name is not None else f"{first_stem}_etc.proteins.faa"
    )
    if output_dir is not None:
        join_path = output_dir / join_path.name
    # Start fresh
    join_path.write_text("")

for file in args.files:
    gbkPath = Path(file)
    suffix = args.suffix if args.suffix is not None else f"_{gbkPath.stem}"
    faaPath = getFaaFromGbk(
        gbkPath,
        getIdFrom=args.seq_id_from,
        prefix=args.prefix,
        suffix=suffix,
    )
    if output_dir is not None:
        target_path = output_dir / faaPath.name
        if faaPath.resolve() != target_path.resolve():
            shutil.move(str(faaPath), str(target_path))
        faaPath = target_path

    print(faaPath)
    if args.join:
        with open(join_path, "a") as jf, open(faaPath) as ff:
            jf.write(ff.read())

if args.join:
    print(f"Joint file: {join_path}")
