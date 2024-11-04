from pathlib import Path
import pandas as pd
import bz2
from pandas import DataFrame
import logging
import pickle
from typing import TextIO


class CoverageData():
    def __init__(
        self,
        coverageTxtFile: Path | None,
        coverageDataBz2File: Path | None = None
    ) -> None:
        self.data: DataFrame | None
        self.txtFile: Path | None = coverageTxtFile
        self.dataBz2File: Path | None = coverageDataBz2File
        if self.dataBz2File is None:
            if self.txtFile is not None:
                logging.info(f'Reading coverage data...')
                self.data = self._readCoverageToDataframe(self.txtFile)
                logging.info(f'Writing to bz2 pickle file for later use...')
                self.dataBz2File = self._zipCoverageData(
                    self.data,
                    self.txtFile.with_suffix('.pickle.bz2')
                )
                logging.info(f'Done writing.')
            else:
                self.data = None
        else:
            assert self.dataBz2File.exists()
            self.data = self._readCoverageFromZip(self.dataBz2File)

    def _readCoverageFromZip(self, zipFile: Path):
        with bz2.open(zipFile, 'rb') as dataFile:
            return pickle.load(dataFile)

    def _readCoverageToDataframe(self, coverageFile: Path) -> DataFrame:
        '''
        chr	pos	ChIP-48h	ChIP-25h	gDNA-25h	gDNA-48h
        NC_003888.3	1	2	5	7	27
        NC_003888.3	2	2	5	9	29
        NC_003888.3	3	2	5	10	34
        NC_003888.3	4	3	5	10	34'''
        with coverageFile.open('r') as cf:
            columns = cf.readline().strip().split('\t')
        columns.remove('chr')
        return pd.read_csv(coverageFile, delimiter='\t',
                           usecols=columns, index_col='pos')

    def _zipCoverageData(self, coverageData: DataFrame, zipFile: Path) -> Path:
        with bz2.open(zipFile, 'wb') as zf:
            pickle.dump(coverageData, zf)
        return zipFile


def read_macs_pileup(
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
