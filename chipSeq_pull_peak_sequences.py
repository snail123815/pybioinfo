from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse
from Bio import SeqIO
import os
from typing import Any
import concurrent.futures
from pathlib import Path
from argparse import ArgumentError

from pyBioinfo_modules.bio_sequences.bio_features import slice_sequence
from pyBioinfo_modules.chipseq.find_and_filter import change_location_to_summit
from pyBioinfo_modules.chipseq.find_and_filter import filter_peaks
from pyBioinfo_modules.chipseq.read_peak_file import read_peak_file


def slice_seq_concurrent(
    source_seqs: list[SeqRecord], peak_info: dict[str, list[str, list[int]]]
) -> list[SeqRecord]:
    """
    peak_info data structure:
    peak_info -> {peak_id: [chr, [start, end]]}
    """
    extracted_seqs = []
    sources = {}
    for seq in source_seqs:
        sources[seq.id] = seq
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        futures = []
        for peak, info in peak_info.items():
            chr = info[0]
            loc = info[1]
            if chr not in sources:
                raise Exception(f"Chromosome {chr} not found in genome file")
            source_seq = sources[chr]
            futures.append(
                executor.submit(slice_sequence, source_seq, loc, peak)
            )
        for future in concurrent.futures.as_completed(futures):
            sliceSeq = future.result()
            extracted_seqs.append(sliceSeq)
    return extracted_seqs


def get_title(title, filter_groups, around_summit):
    if around_summit:
        title = f"{title}_around_summit_{around_summit[0]}_{around_summit[1]}"
    for filter_group in filter_groups:
        filter_key, filter_limites = filter_group
        fls = []
        for limit in filter_limites:
            if limit is not None:
                fls.append(str(limit))
        if filter_key:
            title = f"{title}_{filter_key}"
        if fls:
            title += f"_{'_'.join(fls)}"
    return title


def parse_filter_args_single_group(filter_args: list[str]) -> tuple[str, list]:
    filter_args = [a.lower() for a in filter_args]
    if filter_args[0] == "none":
        filter_method = None
        method_args = [None]
    elif filter_args[0] == "chr":
        filter_method = "chr"
        method_args = filter_args[1:]
    else:
        filter_method = filter_args[0]
        if len(filter_args) > 1:
            try:
                method_args = [int(a) for a in filter_args[1:]]
            except ValueError:
                try:
                    method_args = [float(a) for a in filter_args[1:]]
                except ValueError:
                    raise ArgumentError(
                        None,
                        f"Filter method {filter_method} should have "
                        "integer or float arguments",
                    )
    return filter_method, method_args


def parse_filter_args(filter_args: list[list[str]]) -> list[tuple[str, list]]:
    filter_groups = []
    for filter_arg in filter_args:
        filter_groups.append(parse_filter_args_single_group(filter_arg))
    return filter_groups


def get_peak_info(peak_df: pd.DataFrame) -> dict[str, list[str, list[int]]]:
    peak_info = {}
    for peak_id, row in peak_df.iterrows():
        chr = row["chr"] if "chr" in row else None
        start = int(row["start"])
        end = int(row["end"])
        peak_info[peak_id] = [chr, [start, end]]
    return peak_info


# Main function:
def extract_seq_peaks_to_file(
    genome_file_path: Path,
    peak_file_path: Path,
    output_dir: Path,
    filter_groups: list[list[str, Any]] | None = None,
    around_summit: list[int] | None = None,
    overwrite: bool = False,
):

    title = peak_file_path.stem
    peak_df = read_peak_file(peak_file_path)
    print("*" * 100)
    print(peak_file_path)

    title = get_title(title, filter_groups, around_summit)
    output_fasta = output_dir / f"{title}_seqs.fasta"
    output_table = output_fasta.with_suffix(".xlsx")

    if output_fasta.exists() and not overwrite:
        print(f"The result file exists {output_fasta}, skip.")
        return output_fasta

    peak_df = filter_peaks(peak_df, filter_groups)
    if around_summit:
        peak_df = change_location_to_summit(peak_df, *around_summit)

    if genome_file_path.suffix in [".gb", ".gbk", ".gbff"]:
        genome_seqs = list(SeqIO.parse(genome_file_path, "genbank"))
    elif genome_file_path.suffix in [".fa", ".fasta", ".fna"]:
        genome_seqs = list(SeqIO.parse(genome_file_path, "fasta"))
    else:
        raise Exception(f"Genome file {genome_file_path} not supported")

    peak_seqs = slice_seq_concurrent(genome_seqs, get_peak_info(peak_df))

    print(f"Output file\n{output_fasta}\n{output_table}")
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    SeqIO.write(peak_seqs, output_fasta, "fasta")
    peak_df.to_excel(output_table)
    return output_fasta


def combine_filenames(file_paths):
    assert len(file_paths) >= 2
    file_names = [f.name for f in file_paths]
    p1 = file_paths[0].parent
    assert all(f.parent == p1 for f in file_paths[1:])
    f1 = file_names[0]

    # Get prefix
    prefix_len = 0
    for i, c in enumerate(f1):
        if not all(f[i] == c for f in file_names):
            break
        prefix_len = i + 1

    # Get suffix
    suffix_len = 0
    for i in range(1, len(f1) + 1):
        if not all(f[-i] == f1[-i] for f in file_names):
            break
        suffix_len = i

    # Extract unique parts between prefix and suffix
    unique_parts = []
    for fname in file_names:
        middle = fname[prefix_len : len(fname) - suffix_len]
        unique_parts.append(middle)

    # Combine parts
    combined_middle = "_".join(sorted(unique_parts))

    # Reconstruct filename
    new_name = Path(f1[:prefix_len] + combined_middle + f1[-suffix_len:])
    new_name = new_name.stem + "_combined" + new_name.suffix
    return p1 / new_name


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--genome", type=Path, help="genome file, fastta or genbank"
    )
    parser.add_argument("--files", type=Path, nargs="+")
    parser.add_argument(
        "--around_summit",
        nargs=2,
        help=(
            "Left and right "
            "around summit to extract sequence, will ignore "
            "start and end location in the peak file"
        ),
    )
    parser.add_argument(
        "--filter",
        nargs="+",
        action="append",
        default=[["none"]],
        help="""
    filter method, eg.
    --filter chr NC_0888.1 NC_0889.1 NC_0890.1 --filter length 100 300
    You can use multiple filters by passing multiple arguments
     """,
    )
    parser.add_argument("-o", "--output", type=Path, help="output dir")
    parser.add_argument(
        "-f",
        "--overwrite",
        action="store_true",
        help="overwrite existing files",
    )

    args = parser.parse_args()
    genome = args.genome
    files = args.files
    filter_args = args.filter
    output_dir = args.output
    if output_dir.exists():
        assert output_dir.is_dir()

    filters = parse_filter_args(filter_args)

    output_fastas = []
    for file in files:
        output_fastas.append(
            extract_seq_peaks_to_file(
                genome, file, output_dir, filters, overwrite=args.overwrite
            )
        )

    if len(output_fastas) > 1:
        new_filename = combine_filenames(output_fastas)
        combined_fasta = output_dir / new_filename
        print("*" * 100)
        print("Multiple files processed, generating a combined fasta file")
        with open(combined_fasta, "w") as out:
            for fasta in output_fastas:
                with open(fasta) as f:
                    out.write(f.read())
        print(f"Combined fasta file:\n{combined_fasta}")


if __name__ == "__main__":
    main()
