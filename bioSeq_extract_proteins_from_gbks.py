# This file is licensed under the MIT License


import argparse
import shutil
import re
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
    "--globalprefix",
    help="universal prefix to add to each extracted sequence ID (appended after regex prefix if used)",
    default=None,
)
argparser.add_argument(
    "--globalsuffix",
    help="universal suffix to add to each extracted sequence ID (appended after regex suffix if used)",
    default=None,
)
argparser.add_argument(
    "--use-prefix-regex",
    action="store_true",
    help="Enable extracting prefix from filename using regex",
)
argparser.add_argument(
    "--use-suffix-regex",
    action="store_true",
    help="Enable extracting suffix from filename using regex",
)
argparser.add_argument(
    "--prefix-regex-from-filename",
    help="Regex to extract prefix from filename (first match group used). Default: '^(.*?)\\.|_'",
    default=r"^(.*?)\\.|_",
)
argparser.add_argument(
    "--suffix-regex-from-filename",
    help="Regex to extract suffix from filename (first match group used). Default: '^(.*?)\\.|_'",
    default=r"^(.*?)\\.|_",
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

def extract_from_filename(filename, regex, default=None):
    match = re.match(regex, filename)
    if match and match.groups():
        return match.group(1)
    return default

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
    fname = gbkPath.name
    # Determine prefix
    prefix = ""
    if args.use_prefix_regex:
        regex_prefix = extract_from_filename(fname, args.prefix_regex_from_filename)
        if regex_prefix:
            prefix += regex_prefix
    if args.globalprefix:
        prefix += args.globalprefix
    if not prefix:
        prefix = None
    # Determine suffix
    suffix = ""
    if args.use_suffix_regex:
        regex_suffix = extract_from_filename(fname, args.suffix_regex_from_filename)
        if regex_suffix:
            suffix += f"_{regex_suffix}"
    if args.globalsuffix:
        suffix += args.globalsuffix
    if not suffix:
        suffix = f"_{gbkPath.stem}"
    faaPath = getFaaFromGbk(
        gbkPath,
        getIdFrom=args.seq_id_from,
        prefix=prefix,
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
