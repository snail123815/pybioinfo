import gzip
import io
import logging
import re
import subprocess
from collections import OrderedDict, namedtuple
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from tempfile import NamedTemporaryFile

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from pyBioinfo_modules.basic.decompress import decompFileIfCompressed
from pyBioinfo_modules.wrappers._environment_settings import (
    HMMER_ENV,
    withActivateEnvCmd,
)
from pyBioinfo_modules.wrappers.hmmer_config import (
    BaseHmmerTblFilters,
    HmmerHomologousProtFillters,
)

logger = logging.getLogger(__name__)


def run_hmmsearch():
    raise NotImplementedError("run_hmmsearch is not implemented")


def _find_break_point(target: str, ref_proteome_p: Path):
    """
    Find the location of the target in the reference proteome file
    and the next gene location.
    This is used to resume the search if it was interrupted.
    You need to find what is the last sequence that was processed
    and start from there.

    target: Fasta header of the target protein, e.g. ">sp|P0A7D3|ARGH_ECOLI"
    ref_proteome_p: Path to the reference proteome file
        A fasta file containing the sequences to search against.
        This file is usually large, so it is compressed with gzip.
        The file is seekable, so we can use the tell() method to find
        the location of the target protein in the file.

    https://github.com/snail123815/proj-coll-Emtinan-aminoacylation-phenotype-gene-correlation/blob/main/step_1.1_get_location.py
    """

    # Check the first appearance
    target_seek_loc = 0
    with gzip.open(ref_proteome_p, "rt") as rpp:
        assert rpp.seekable()

        l = rpp.readline()
        while l:
            if l.startswith(target):
                target_seek_loc = rpp.tell()
                break
            l = rpp.readline()

    # Find the next id location
    next_id_found = False
    with gzip.open(ref_proteome_p, "rt") as rpp:
        rpp.seek(target_seek_loc)
        l = rpp.readline()
        next_seek_loc = target_seek_loc
        while l:
            next_seek_loc = rpp.tell()
            if l.startswith(">"):
                next_id_found = True
                break
            l = rpp.readline()

    if not next_id_found:
        raise ValueError(
            f"No next id found in the file. Maybe {target} is "
            "already the last one."
        )

    # Check the result
    with gzip.open(ref_proteome_p, "rt") as rpp:
        rpp.seek(next_seek_loc)
        for i in range(5):
            logger.debug(rpp.readline())

    return next_seek_loc


def run_jackhmmer(
    query_fasta: Path | StringIO,
    target_fasta_path: Path,
    domtblout_path: Path,
    jackhmmer_cfg: BaseHmmerTblFilters = BaseHmmerTblFilters(),
    cpus: int = 8,
):
    # Write sequences from query_fasta to an in-memory text stream using StringIO

    # Open the query_fasta file and parse as FASTA
    # Alternatively, if you already have SeqRecord objects, use them directly.
    if isinstance(query_fasta, Path):
        query_fasta = decompFileIfCompressed(query_fasta, to_temp=True)[0]
        with open(query_fasta, "rt") as f:
            query_content = f.read()
    else:
        query_content = query_fasta.read()

    jackhmmer_cmd = (
        "jackhmmer"
        f" -E {jackhmmer_cfg.T_E}"
        f" --incE {jackhmmer_cfg.T_INCE}"
        f" --domE {jackhmmer_cfg.T_DOME}"
        f" --incdomE {jackhmmer_cfg.T_INCDOME}"
        f" --cpu {cpus} --domtblout {domtblout_path} - {target_fasta_path}"
    )

    jackhmmer_run = subprocess.run(
        withActivateEnvCmd(jackhmmer_cmd, HMMER_ENV),
        input=query_content.encode("utf-8"),
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )

    if jackhmmer_run.returncode != 0:
        logger.error(
            f"Error running jackhmmer: {jackhmmer_run.stderr.decode('utf-8')}"
        )
        raise RuntimeError(
            f"Jackhmmer failed with return code {jackhmmer_run.returncode}"
        )
    else:
        return 1


def full_proteome_jackhmmer(
    query_fasta_path,
    target_fasta_path,
    domtblout_path,
    prots_per_group=10,
    cpus=8,
    hmmer_filters: BaseHmmerTblFilters = HmmerHomologousProtFillters,
):
    """
    Run jackhmmer on the full proteome, usually to find homologous proteins
    for a set of query proteins.
    The query proteins are split into groups, each group is searched
    separately to reduce memory usage.
    The results are concatenated into a single domtblout file.
    query_fasta_path: Path to the query fasta file
        A (gzipped) fasta file containing the sequences to search. Query.
    target_fasta_path: Path to the database file
        A (gzipped) fasta file containing the sequences to search against, make
        sure there is no duplicated headers in this file
    domtblout_path: Path to the output file
        Format is --domtblout, it is concatenated from each group searches.
        Please do exclude # lines in the middle when parsing.
    """
    query_temp = decompFileIfCompressed(query_fasta_path, to_temp=True)[0]
    with open(query_temp, "rt") as qh:
        query_seqs = list(SeqIO.parse(qh, "fasta"))
    output_handle = domtblout_path.open("wt")
    header_written = False
    last_header_skipped = False
    for i in range(0, len(query_seqs), prots_per_group):
        group_seqs = query_seqs[i : min(i + prots_per_group, len(query_seqs))]
        group_seqs_io = io.StringIO()
        SeqIO.write(group_seqs, group_seqs_io, "fasta")
        group_seqs_io.seek(0)
        group_seqs_out = NamedTemporaryFile()
        run_jackhmmer(
            group_seqs_io,
            target_fasta_path,
            Path(group_seqs_out.name),
            jackhmmer_cfg=hmmer_filters,
            cpus=cpus,
        )
        with open(group_seqs_out.name, "rt") as group_seqs_in:
            for line in group_seqs_in:
                if i == 0:
                    # Write the header line only once
                    if line.startswith("#"):
                        if not header_written:
                            output_handle.write(line)
                    else:
                        header_written = True
                elif i >= len(query_seqs) / prots_per_group:
                    if line.startswith("#"):
                        if last_header_skipped:
                            output_handle.write(line)
                    else:
                        last_header_skipped = True
                if line.startswith("#"):
                    continue
                output_handle.write(line)
        group_seqs_out.close()
        group_seqs_io.close()
    output_handle.close()


def calculate_domain_coverage(
    line: dict,
    # each row by _asdict() from domtbl DataFrame, generated by read_domtbl()
    dom_locations: list[int],
    hmmer_filters: BaseHmmerTblFilters = BaseHmmerTblFilters(),
) -> tuple[bool, list[int], float, float]:
    """
    Run within each single protein lines, will update dom_cov_regions list.
    Only the alignment of domains with E values lower than threshold
    will be included in the coverage calculation
    Each coverage will be converted to range of positions (int)
    Combine all positions -> deduplicate -> count = lenth of covered region

    Note: It does not care which protein ID it is, just process the line and
        update the dom_cov_regions.

    line: dict, a single line from domtblout as dict
    dom_cov_regions: list of int, covered **positions** so far
    hmmer_filters: BaseHmmerTblFilters, configuration for filtering

    return:
        is_end_dom: bool, whether this is the end of all domains for this protein
        dom_cov_regions: list of int, updated list of covered positions
        dom_covq: float, coverage on query, only calculate when
            is_end_dom is True and
            dom_i_E <= GATHER_T_DOME
        dom_covt: float, coverage on target, only calculate when
            is_end_dom is True and
            dom_i_E <= GATHER_T_DOME

    """
    dom_covq = dom_covt = 0
    is_end_dom = False
    if line["dom_total"] == 1:  # Single domain
        is_end_dom = True
        if line["dom_i_E"] <= hmmer_filters.GATHER_T_DOME:
            dom_locations = list(range(line["ali_from"], line["ali_to"] + 1))
            cov_len = line["ali_to"] - line["ali_from"] + 1
            dom_covq = cov_len / line["qlen"]
            dom_covt = cov_len / line["tlen"]
    else:  # multi domains
        if line["dom_n"] == 1:  # first domain
            dom_locations = []  # should be redundant, but reset to be safe
        if line["dom_i_E"] <= hmmer_filters.GATHER_T_DOME:
            dom_locations.extend(
                list(range(line["ali_from"], line["ali_to"] + 1))
            )
        if line["dom_n"] == line["dom_total"]:  # end of all domains
            is_end_dom = True
            cov_len = len(set(dom_locations))
            dom_covq = cov_len / line["qlen"]
            dom_covt = cov_len / line["tlen"]
    return is_end_dom, dom_locations, dom_covq, dom_covt


def remove_duplicates(file_p: Path) -> Path:
    no_dup_path = file_p.parent / f"{file_p.stem}_nodup{file_p.suffix}"
    logger.info(f"Removing duplicates from file {file_p} -> {no_dup_path}")
    ddup = subprocess.run(
        f"awk '!seen[$0]++' {file_p} > {no_dup_path}", shell=True
    )
    ddup.check_returncode()
    return no_dup_path


def read_domtbl(
    domtbl_p: Path,
    hmmer_filters: BaseHmmerTblFilters = BaseHmmerTblFilters(),
    domtbl_splitter=re.compile(r" +"),  # Split by spaces
) -> pd.DataFrame:

    logger.info(f"Counting lines in domtblout from jackhmmer: {domtbl_p}")
    domtbl_len = int(
        subprocess.run(["wc", "-l", domtbl_p], capture_output=True)
        .stdout.decode()
        .split(" ")[0]
    )
    logger.info(f"{domtbl_len} lines in total.")

    # Initialize domtbl_dict with header keys
    domtbl_dict = {}
    header = [
        "tp",
        "qp",
        "tlen",
        "qlen",
        "ali_from",
        "ali_to",
        "dom_i_E",
        "dom_n",
        "dom_total",
        "full_E",
        "anno",
    ]
    for h in header:
        domtbl_dict[h] = []

    t_e = hmmer_filters.GATHER_T_E
    len_diff = hmmer_filters.LEN_DIFF
    with domtbl_p.open("rt") as dt:
        for l in tqdm(dt, desc="Reading file", total=domtbl_len):
            if l.startswith("#"):
                continue
            line_list = domtbl_splitter.split(l.strip())

            full_E = float(line_list[6])
            if full_E > t_e:
                continue
            else:
                tlen = int(line_list[2])
                qlen = int(line_list[5])
                if abs(tlen - qlen) / min(tlen, qlen) > len_diff:
                    continue
                try:
                    domtbl_dict["tp"].append(line_list[0])
                    domtbl_dict["qp"].append(line_list[3])
                    domtbl_dict["tlen"].append(tlen)
                    domtbl_dict["qlen"].append(qlen)
                    domtbl_dict["ali_from"].append(int(line_list[17]))
                    domtbl_dict["ali_to"].append(int(line_list[18]))
                    domtbl_dict["dom_i_E"].append(float(line_list[12]))
                    domtbl_dict["dom_n"].append(int(line_list[9]))
                    domtbl_dict["dom_total"].append(int(line_list[10]))
                    domtbl_dict["full_E"].append(full_E)
                    domtbl_dict["anno"].append(" ".join(line_list[22:]))
                except ValueError as ve:
                    logger.error(f"Error parsing line: {l}")
                    raise ve
    logger.info("Converting data to dataframe...")
    domtbl_df = pd.DataFrame(domtbl_dict)
    logger.info("Done.")
    return domtbl_df


def parse_dom_table_mt(
    domtbl_df,  # pd.DataFrame, domtbl DataFrame from read_domtbl()
    cpus=8,
    hmmer_filters: BaseHmmerTblFilters = BaseHmmerTblFilters(),
) -> pd.DataFrame:
    """
    Multi-threaded processing of domtbl DataFrame

    For full proteome jackhmmer results, please use HmmerHomologousProtFillters
    as hmmer_filters.

    hmmer_filters: BaseHmmerTblFilters, configuration for filtering
    return: concatenated DataFrame of all results
        Each row shows the result of a single query protein, with coverage
        on both query and target, and other information.
    """

    def _process_single_qp(
        domtbl_q: list[OrderedDict],  # Gathered hits of a single query protein
    ):
        """
        Process hits of a single query protein, will calculate coverage
        and filter based on coverage thresholds.

        domtbl_q: list of OrderedDict, each dict is a line from domtblout

        Return: DataFrame of the processed results
        """
        # Get query protein id
        qp = domtbl_q[0]["qp"]
        # Initialize coverage regions and match dict
        dom_cov_regions = []
        match_dict = {
            "Query": [],
            "Target strain": [],
            "Target protein": [],
            "Coverage on Query": [],
            "Coverage on Target": [],
            "Expect protein": [],
            "Target description": [],
        }

        for row in domtbl_q:
            is_end_dom, dom_cov_regions, dom_covq, dom_covt = (
                calculate_domain_coverage(row, dom_cov_regions, hmmer_filters)
            )
            if is_end_dom:
                if (
                    dom_covq < hmmer_filters.GATHER_T_COV
                    or dom_covt < hmmer_filters.GATHER_T_COV
                ):
                    continue
                dom_cov_regions = []
            else:
                continue

            target_strain = row["tp"].split("_")[0]
            target_protein = row["tp"][len(target_strain) + 1 :]
            match_dict["Query"].append(qp)
            match_dict["Target strain"].append(target_strain)
            match_dict["Target protein"].append(target_protein)
            match_dict["Coverage on Query"].append(dom_covq)
            match_dict["Coverage on Target"].append(dom_covt)
            match_dict["Expect protein"].append(row["full_E"])
            match_dict["Target description"].append(
                row["anno"].replace(target_protein, "").strip()
            )
        return pd.DataFrame(match_dict)

    previous_qp = ""  # To group by query protein
    domtbl_q: list[OrderedDict] = []
    with ProcessPoolExecutor(cpus) as executer:
        futures = []
        for domrow in tqdm(
            domtbl_df.itertuples(),
            desc="Submitting jobs",
            total=domtbl_df.shape[0],
        ):
            dom_dict = domrow._asdict()
            if domrow.qp == previous_qp:
                domtbl_q.append(dom_dict)
            else:
                # This is next qp, now submit previous one
                if len(domtbl_q) > 0:
                    futures.append(
                        executer.submit(_process_single_qp, domtbl_q)
                    )
                # empty q
                domtbl_q = [dom_dict]
                previous_qp = domrow.qp
        # submit the last qp
        if len(domtbl_q) > 0:
            futures.append(executer.submit(_process_single_qp, domtbl_q))
        results = []
        for future in tqdm(futures, desc="Processing per protein hits"):
            results.append(future.result())
    return pd.concat(results, ignore_index=True)
