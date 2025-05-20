from pathlib import Path

import pandas as pd


def _read_macs_peakcalls(file: Path) -> pd.DataFrame:
    """Reads a peak file from direct xls output of MACS program
    and returns a DataFrame.
    """
    columns = [
        "chr",
        "name",
        "start",
        "end",
        "abs_summit",
        "-log10(pvalue)",
        "-log10(qvalue)",
        "fold_enrichment",
    ]
    data = pd.read_csv(
        file,
        delimiter="\t",
        comment="#",
        index_col="name",
        usecols=columns,
        dtype={
            "chr": str,
            "name": str,
            "start": int,
            "end": int,
            "abs_summit": int,
            "-log10(pvalue)": float,
            "-log10(qvalue)": float,
            "fold_enrichment": float,
        },
    )
    return data


def _read_macs_common_peaks_bed_file(file: Path) -> pd.DataFrame:
    """Reads macs common peak bed file"""
    with file.open("rt") as f:
        first_line = f.readline()
    if first_line.startswith("track"):
        skiprows = 1
    else:
        skiprows = 0
    if "common" in file.name:
        columns = ["start", "end", "likely_difference"]
    else:
        columns = ["start", "end", "value"]
    data = pd.read_csv(
        file,
        delimiter="\t",
        skiprows=skiprows,
        header=None,
        usecols=[1, 2, 3, 4],
        index_col=2,
    )
    data.columns = columns
    data.index.name = "name"
    # Change datatype after reading:
    # 'start' and 'end' as integers and the third column as float.
    data = data.astype({"start": int, "end": int, columns[2]: float})
    return data


def read_common_peaks_tsv(file):
    """
    ! not used by general script
    Read common peaks tsv file, used in Gbn project
    Three kind of dataframe generated in that project:
    name, start, end, abs_summit, fold_enrichment
    name, start, end, abs_summit, fold_enrichment_A, fold_enrichment_B
    name, start, end, log10_likely (of the peak presented in this file)
    thresh is the threshold of the -log10(pvalue)
    """
    columns = [
        "name",
        "start",
        "end",
        "abs_summit",
        "fold_enrichment_A",
        "fold_enrichment_B",
    ]
    data = pd.read_csv(file, delimiter="\t", usecols=columns, index_col="name")
    return data


def read_peak_file(file: Path) -> pd.DataFrame:
    # Read peak file, always have an index column "name", which is unique
    # for each peak
    if file.suffix == ".xls":
        data = _read_macs_peakcalls(file)
    elif file.suffix(".bed"):  # different peaks called by Macs2
        data = _read_macs_common_peaks_bed_file(file)
    else:
        raise Exception(f"File {file} not recognized")
    data.columns = data.columns.str.lower()
    return data
