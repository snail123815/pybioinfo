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
    CONDAEXE,
    HMMER_ENV,
    SHELL,
    withActivateEnvCmd,
)
from pyBioinfo_modules.wrappers.hmmer_config import (
    BaseHmmerConfig,
    HmmerHomologousProtConfig,
)

logger = logging.getLogger(__name__)


def run_hmmsearch():
    return


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
            print(rpp.readline())

    return next_seek_loc


def run_jackhmmer(
    query_fasta: Path | StringIO,
    target_fasta_path: Path,
    domtblout_path: Path,
    jackhmmer_cfg: BaseHmmerConfig = BaseHmmerConfig(),
    cpus: int = 8,
):
    # Write sequences from query_fasta to an in-memory text stream using StringIO

    # Open the query_fasta file and parse as FASTA
    # Alternatively, if you already have SeqRecord objects, use them directly.
    if isinstance(query_fasta, Path):
        query_fasta = decompFileIfCompressed(query_fasta, to_temp=True)[0]
        query_input = open(query_fasta, "rt")
    else:
        query_input = query_fasta

    jackhmmer_cmd = (
        "jackhmmer"
        f" -E {jackhmmer_cfg.T_E}"
        f" --incE {jackhmmer_cfg.T_INCE}"
        f" --domE {jackhmmer_cfg.T_DOME}"
        f" --incdomE {jackhmmer_cfg.T_INCDOME}"
        f" --cpu {cpus} --domtblout {domtblout_path} - {target_fasta_path}"
    )

    jackhmmer_run = subprocess.run(
        withActivateEnvCmd(jackhmmer_cmd, HMMER_ENV, CONDAEXE, SHELL),
        input=query_input.read().encode("utf-8"),
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
    hhpc: BaseHmmerConfig = HmmerHomologousProtConfig(),
):
    """
    Run jackhmmer on the full proteome
    ref_proteome_path: Path to the reference proteome file
        A gzipped fasta file containing the sequences to search. Query.
    db_f: Path to the database file
        A gzipped fasta file containing the sequences to search against, make
        sure there is no duplicated headers in this file
    domtblout_path: Path to the output file
        Format is --domtblout, it is concatenated from each group searches.
        Please do exclude # lines in the middle when parsing.
    """
    query_seqs = list(
        SeqIO.parse(
            decompFileIfCompressed(query_fasta_path, to_temp=True)[0], "fasta"
        )
    )
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
            jackhmmer_cfg=hhpc,
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


def filter_domtblout(
    domtbl_df: pd.DataFrame, include_ali_coords: list[tuple[int, int]] = []
):
    """ """
    pass


def cal_cov(
    line: dict,
    dom_cov_regions: list[int],
    hhpc: HmmerHomologousProtConfig = HmmerHomologousProtConfig(),
):
    """
    Run within each single protein lines
    Only the alignment of domains with E values lower than threshold
    will be included in the coverage calculation
    Each coverage will be converted to range of positions (int)
    Combine all positions -> deduplicate -> count = lenth of covered region

    """
    dom_covq = dom_covt = 0
    is_end_dom = False
    if line["dom_total"] == 1:  # Single domain
        is_end_dom = True
        if line["dom_i_E"] <= hhpc.GATHER_T_DOME:
            cov_len = line["ali_to"] - line["ali_from"] + 1
            dom_covq = cov_len / line["qlen"]
            dom_covt = cov_len / line["tlen"]
    else:  # multi domains
        # if line['dom_n'] == 1:  # first multi domain
        #     dom_cov_regions = []
        if line["dom_i_E"] <= hhpc.GATHER_T_DOME:
            dom_cov_regions.extend(
                list(range(line["ali_from"], line["ali_to"] + 1))
            )
        if line["dom_n"] == line["dom_total"]:  # end of all domains
            is_end_dom = True
            cov_len = len(set(dom_cov_regions))
            dom_covq = cov_len / line["qlen"]
            dom_covt = cov_len / line["tlen"]
    return is_end_dom, dom_cov_regions, dom_covq, dom_covt


def remove_duplicates(file_p: Path) -> Path:
    no_dup_path = file_p.parent / f"{file_p.stem}_nodup{file_p.suffix}"
    print(f"Removing duplicates from file {file_p} -> {no_dup_path}")
    ddup = subprocess.run(
        f"awk '!seen[$0]++' {file_p} > {no_dup_path}", shell=True
    )
    ddup.check_returncode()
    return no_dup_path


def read_domtbl(
    domtbl_p: Path,
    hhpc: HmmerHomologousProtConfig = HmmerHomologousProtConfig(),
) -> pd.DataFrame:
    def parse_one_line(
        l,
        domtbl_dict,
        t_e=hhpc.GATHER_T_E,
        len_diff=hhpc.LEN_DIFF,
        splitter=re.compile(r" +"),
    ):
        line_list = splitter.split(l.strip())

        full_E = float(line_list[6])
        if full_E > t_e:
            return
        else:
            tlen = int(line_list[2])
            qlen = int(line_list[5])
            if abs(tlen - qlen) / min(tlen, qlen) > len_diff:
                return
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
                print(l)
                raise ve

    print(f"Counting lines in domtblout from jackhmmer: {domtbl_p}")
    domtbl_len = int(
        subprocess.run(["wc", "-l", domtbl_p], capture_output=True)
        .stdout.decode()
        .split(" ")[0]
    )
    print(f"{domtbl_len} lines in total.")
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
    with domtbl_p.open("rt") as dt:
        for l in tqdm(dt, desc="Reading file", total=domtbl_len):
            if l.startswith("#"):
                continue
            parse_one_line(l, domtbl_dict)
    print("Converting data to dataframe...", end="")
    domtbl_df = pd.DataFrame(domtbl_dict)
    print("Done.")
    return domtbl_df


def process_single_query(
    domtbl_q: list[OrderedDict],
    hhpc: HmmerHomologousProtConfig = HmmerHomologousProtConfig(),
):
    qp = domtbl_q[0]["qp"]
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
        is_end_dom, dom_cov_regions, dom_covq, dom_covt = cal_cov(
            row, dom_cov_regions
        )
        if is_end_dom:
            if dom_covq < hhpc.GATHER_T_COV or dom_covt < hhpc.GATHER_T_COV:
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


def parse_dom_table_mt(domtbl_df, cpus=8):
    previous_qp = ""
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
                # submit previous one
                if len(domtbl_q) > 0:
                    futures.append(
                        executer.submit(process_single_query, domtbl_q)
                    )
                # empty q
                domtbl_q = [dom_dict]
                previous_qp = domrow.qp
        # submit the last one
        if len(domtbl_q) > 0:
            futures.append(executer.submit(process_single_query, domtbl_q))
        results = []
        for future in tqdm(futures, desc="Processing per protein hits"):
            results.append(future.result())
    return pd.concat(results, ignore_index=True)
