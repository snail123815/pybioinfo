import bz2
import logging
import lzma
import pickle
from pathlib import Path
from typing import TextIO

import numpy as np
import pandas as pd
from pandas import DataFrame


class CoverageData:
    def __init__(
        self,
        coverageTxtFile: Path | None,
        coverageDataBz2File: Path | None = None,
    ) -> None:
        self.data: DataFrame | None
        self.txtFile: Path | None = coverageTxtFile
        self.dataBz2File: Path | None = coverageDataBz2File
        if self.dataBz2File is None:
            if self.txtFile is not None:
                logging.info(f"Reading coverage data...")
                self.data = self._readCoverageToDataframe(self.txtFile)
                logging.info(f"Writing to bz2 pickle file for later use...")
                self.dataBz2File = self._zipCoverageData(
                    self.data, self.txtFile.with_suffix(".pickle.bz2")
                )
                logging.info(f"Done writing.")
            else:
                self.data = None
        else:
            assert self.dataBz2File.exists()
            self.data = self._readCoverageFromZip(self.dataBz2File)

    def _readCoverageFromZip(self, zipFile: Path):
        with bz2.open(zipFile, "rb") as dataFile:
            return pickle.load(dataFile)

    def _readCoverageToDataframe(self, coverageFile: Path) -> DataFrame:
        """
        chr	pos	ChIP-48h	ChIP-25h	gDNA-25h	gDNA-48h
        NC_003888.3	1	2	5	7	27
        NC_003888.3	2	2	5	9	29
        NC_003888.3	3	2	5	10	34
        NC_003888.3	4	3	5	10	34"""
        with coverageFile.open("r") as cf:
            columns = cf.readline().strip().split("\t")
        columns.remove("chr")
        return pd.read_csv(
            coverageFile, delimiter="\t", usecols=columns, index_col="pos"
        )

    def _zipCoverageData(self, coverageData: DataFrame, zipFile: Path) -> Path:
        with bz2.open(zipFile, "wb") as zf:
            pickle.dump(coverageData, zf)
        return zipFile


def _read_one_macs_pileup(
    pileup: TextIO, tr_start: int, tr_end: int
) -> list[tuple[int, float]]:
    # Read the pileup file, a range per line
    tr_range_pileup = []
    for line in pileup:
        _, start, end, value = line.strip().split("\t")
        end = int(end)
        if end <= tr_start:
            continue
        start = int(start)
        if start > tr_end:
            break
        value = float(value)
        tr_range_pileup.append((start, end, value))
    # Fill in the the range per base
    tr_perbase_pileup = []
    for start, end, value in tr_range_pileup:
        effective_range = range(max(start, tr_start), min(end, tr_end + 1))
        for i in effective_range:
            tr_perbase_pileup.append((i, value))
    return tr_perbase_pileup


def read_macs_pileup(
    macsOutputPath: Path, tr_start: int, tr_end: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Reads and processes the input arguments for the chipSeq plot pileup
    comparisons.

    Args:
        macsOutputPath (Path): The path to the macs output dir.

        tr_start (int): The start position of the genomic region to plot.
        tr_end (int): The end position of the genomic region to plot.

    Returns:
        tuple: A tuple containing:
        The control data array with genomic positions and pileup values.
            - tr_control_data (np.ndarray): The control data
            - tr_treat_data (np.ndarray): The treatment data

    Raises:
        ValueError: If the region format is invalid or the gene is not found
        in the genome file.

    The function performs the following steps:
    1. Reads the genome file with annotations.
    2. If a region is specified, it parses the region and extracts the start
       and end positions.
    3. If a gene is specified, it finds the gene in the genome annotations and
       calculates the start and end positions based on the flanking region.
    """
    # Find all matching files
    pileup_files = [
        file
        for file in macsOutputPath.expanduser().glob("*.*")
        if file.name.endswith(".bdg") or file.name.endswith(".bdg.xz")
    ]

    assert len(pileup_files) == 2, (
        "There should be two and only two .bdg files in the macs output dir."
        f" Found {len(pileup_files)} files." + str(pileup_files)
    )

    control_lambda_path = None
    treat_pileup_path = None

    for f in pileup_files:
        if "_control_lambda" in f.name:
            control_lambda_path = f
        elif "_treat_pileup" in f.name:
            treat_pileup_path = f
        else:
            raise ValueError(
                "Unknown file: "
                + f
                + ". Expected either _control_lambda or _treat_pileup"
            )

    assert control_lambda_path, (
        "Control pileup file not found. "
        "Please make sure the file name contains '_control_lambda'."
    )
    assert treat_pileup_path, (
        "Experimental pileup file not found. "
        "Please make sure the file name contains '_treat_pileup'."
    )

    if control_lambda_path.suffix == ".xz":
        control_lambda = lzma.open(control_lambda_path, "rt")
    else:
        control_lambda = open(control_lambda_path, "r")
    if treat_pileup_path.suffix == ".xz":
        treat_pileup = lzma.open(treat_pileup_path, "rt")
    else:
        treat_pileup = open(treat_pileup_path, "r")

    tr_control_data = np.array(
        _read_one_macs_pileup(control_lambda, tr_start, tr_end)
    )
    tr_treat_data = np.array(
        _read_one_macs_pileup(treat_pileup, tr_start, tr_end)
    )
    control_lambda.close()
    treat_pileup.close()

    return tr_control_data, tr_treat_data
