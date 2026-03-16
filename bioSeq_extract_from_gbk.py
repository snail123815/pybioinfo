# This script is licensed under the MIT License

import argparse
import re
import shutil
from pathlib import Path
from typing import Literal

from pyBioinfo_modules.bio_sequences.features_from_gbk import getCdsFromGbk
from pyBioinfo_modules.bio_sequences.features_from_gbk import getFaaFromGbk


def _extract_from_filename(filename, regex, default=None):
    match = re.match(regex, filename)
    if match and match.groups():
        return match.group(1)
    return default


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", help="genbank file(s)")
    parser.add_argument(
        "-t",
        "--target_type",
        choices=["cds", "proteins"],
        default="cds",
        help="type of sequences to extract (default: cds)",
    )
    parser.add_argument(
        "-s",
        "--seq_id_from",
        help="feature name that will be used as sequence id",
        default="locus_tag/protein_id",
    )
    parser.add_argument(
        "--prefix",
        help="prefix to add to each extracted sequence ID",
        default=None,
    )
    parser.add_argument(
        "--suffix",
        help="suffix to add to each extracted sequence ID",
        default=None,
    )
    parser.add_argument(
        "--globalprefix",
        help=(
            "universal prefix to add to each extracted sequence ID "
            "(appended after regex prefix if used)"
        ),
        default=None,
    )
    parser.add_argument(
        "--globalsuffix",
        help=(
            "universal suffix to add to each extracted sequence ID "
            "(appended after regex suffix if used)"
        ),
        default=None,
    )
    parser.add_argument(
        "--use-prefix-regex",
        action="store_true",
        help="Enable extracting prefix from filename using regex",
    )
    parser.add_argument(
        "--use-suffix-regex",
        action="store_true",
        help="Enable extracting suffix from filename using regex",
    )
    parser.add_argument(
        "--prefix-regex-from-filename",
        help="Regex to extract prefix from filename (first match group used). Default: '^(.*?)\\.|_'",
        default=r"^(.*?)\\.|_",
    )
    parser.add_argument(
        "--suffix-regex-from-filename",
        help="Regex to extract suffix from filename (first match group used). Default: '^(.*?)\\.|_'",
        default=r"^(.*?)\\.|_",
    )
    parser.add_argument(
        "-j",
        "--join",
        action="store_true",
        help="write all extracted sequences into a single joint output file",
    )
    parser.add_argument(
        "--join_name",
        help=(
            "name of the joint output file. "
            f"Defaults to the first input file stem + '_etc.cds.fna' or '_etc.proteins.faa' based on target type"
        ),
        default=None,
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        help="output directory for extracted sequences",
        default=None,
    )
    args = parser.parse_args()
    if args.seq_id_from == "locus_tag/protein_id": # set default based on target type
        args.seq_id_from = "locus_tag" if args.target_type == "cds" else "protein_id"
    return args


def run_extract_cli():
    args = arg_parser()

    output_dir = Path(args.output_dir) if args.output_dir else None
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)

    extractor = getCdsFromGbk if args.target_type == "cds" else getFaaFromGbk
    output_suffix = "cds.fna" if args.target_type == "cds" else "proteins.faa"

    join_path = None
    if args.join:
        first_stem = Path(args.files[0]).stem
        join_path = Path(
            args.join_name
            if args.join_name is not None
            else f"{first_stem}_etc.{output_suffix}"
        )
        if output_dir is not None:
            join_path = output_dir / join_path.name
        join_path.write_text("", encoding="utf-8")

    for file in args.files:
        gbk_path = Path(file)
        fname = gbk_path.name

        prefix = ""
        if args.use_prefix_regex:
            regex_prefix = _extract_from_filename(fname, args.prefix_regex_from_filename)
            if regex_prefix:
                prefix += regex_prefix
        if args.globalprefix:
            prefix += args.globalprefix
        if args.prefix:
            prefix += args.prefix
        if not prefix:
            prefix = None

        suffix = ""
        if args.use_suffix_regex:
            regex_suffix = _extract_from_filename(fname, args.suffix_regex_from_filename)
            if regex_suffix:
                suffix += f"_{regex_suffix}"
        if args.globalsuffix:
            suffix += args.globalsuffix
        if args.suffix:
            suffix += args.suffix
        if not suffix:
            suffix = f"_{gbk_path.stem}"

        out_path = extractor(
            gbk_path,
            getIdFrom=args.seq_id_from,
            prefix=prefix,
            suffix=suffix,
        )
        if output_dir is not None:
            target_path = output_dir / out_path.name
            if out_path.resolve() != target_path.resolve():
                shutil.move(str(out_path), str(target_path))
            out_path = target_path

        print(out_path)
        if join_path is not None:
            with open(join_path, "a", encoding="utf-8") as jf, open(
                out_path, encoding="utf-8"
            ) as ff:
                jf.write(ff.read())

    if join_path is not None:
        print(f"Joint file: {join_path}")

    return 0

if __name__ == "__main__":
    raise SystemExit(run_extract_cli())