# For RNA-Seq alignment, with each file or pair of files per sample. Out
# put aligned file per sample.


import time
import argparse
import logging
from pathlib import Path
from pyBioinfo_modules.basic.decompress import splitStemSuffixIfCompressed
from pyBioinfo_modules.basic.basic import getTimeStr, timeDiffStr
from pyBioinfo_modules.wrappers.bowtie2 import buildBowtie2idx, runBowtie2

logger = logging.getLogger(__name__)


def arg_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--raw", required=True, type=Path, nargs="+", help="path(s) to raw data"
    )
    parser.add_argument(
        "--out", required=True, type=Path, help="path to output alignment"
    )
    parser.add_argument(
        "--genome",
        nargs="+",
        required=True,
        type=Path,
        help="path to genome file(s), also supports indexed genome (*.bt2)",
    )
    parser.add_argument(
        "--isPe", action="store_true", help="set if you have pairend"
    )
    parser.add_argument(
        "--pesuffix",
        nargs=2,
        help=(
            "suffix to pairend file name. Will impute by default."
            " eg. a_1_.fq.gz and a_2_.fq.gz, --pesuffix _1_ _2_"
        ),
    )
    parser.add_argument(
        "--ncpu", default=1, type=int, help="number of cpu to use"
    )
    parser.add_argument(
        "--sampleNames",
        nargs="+",
        help=(
            "List of sample names."
            " Set when sample names are random string"
            " which will unpair the samples after removing pairend suffix"
        ),
    )
    parser.add_argument(
        "--dryRun",
        action="store_true",
        help="perform a dry run without executing Bowtie2",
    )
    return parser


def get_read_files_per_sample(
    parent_dir: list[Path],
    sample_names: list[str] | None = None,
    isPe: bool = False,
) -> tuple[dict[str, list[Path]], list[str], str]:
    # Gether all files, max depth 2
    file_paths: list[Path] = []
    for d in parent_dir:
        for n in d.iterdir():
            if n.is_file():
                file_paths.append(n)
            elif n.is_dir():
                file_paths.extend(f for f in n.iterdir() if f.is_file())
    for f in file_paths.copy():
        ext = splitStemSuffixIfCompressed(f)[1]
        if len(ext.split(".")) < 2 or ext.split(".")[1] not in ["fastq", "fq"]:
            file_paths.remove(f)
    assert len(file_paths) > 0, f"Files not found in {parent_dir}"
    file_paths = sorted(file_paths)
    file_parts = [
        splitStemSuffixIfCompressed(f, fullSuffix=True) for f in file_paths
    ]
    file_names = [fp[0] for fp in file_parts]
    file_fullexts = [fp[1] for fp in file_parts]
    assert len(set(file_fullexts)) == 1, (
        "All files should have same format, "
        f"yet multiple found {set(file_fullexts)}."
    )
    assert len(set(file_names)) == len(file_names), (
        "File names has duplicates"
        ", maybe from different dirs?\n"
        f"{set([fn for fn in file_names if file_names.count(fn) > 1])}"
    )

    if sample_names is not None:
        samples = sample_names
    else:
        samples = []

    sample_file_dict: dict[str, list[Path]] = {}
    peSfx: list[str] = []
    if isPe:
        peSfx = imputePeSuffix(file_names)
        if samples == []:
            samples = (
                sorted(list(set(f[: -len(peSfx[0])] for f in file_names)))
                if len(samples) == 0
                else samples
            )
            assert (len(file_paths) / len(samples)) % 2 == 0
        for s in samples:
            fns = [file_paths[i] for i, fn in enumerate(file_names) if s in fn]
            assert len(fns) > 0, f"Did not find files for sample {s}"
            assert len(fns) % 2 == 0, (
                f"Number of files for sample {s} is not paired.\n" f"{fns}"
            )
            sample_file_dict[s] = fns
    else:  # single end reads
        if samples == []:
            for fn, fp in zip(file_names, file_paths):
                sample_file_dict[fn] = [fp]
            samples = file_names
        elif len(samples) == 1:
            sample_file_dict[samples[0]] = file_paths
        else:
            for s in samples:
                fps = [
                    file_paths[i] for i, fn in enumerate(file_names) if s in fn
                ]
                assert len(fps) > 0, f"Did not find files for sample {s}"
                sample_file_dict[s] = fps

    return sample_file_dict, peSfx, file_fullexts[0]


def align_multiple_raw_bowtie2(
    raw: list[Path],
    sampleNames: list[str] | None,
    out: Path,
    isPe: bool,
    genomes: list[Path],
    ncpu: int,
    dryRun: bool,
):
    # Start logic
    logger.info("=" * 20 + getTimeStr() + "=" * 20)
    if dryRun:
        logger.debug("=" * 20 + "Dry run, will not execute bowtie2")

    tInit = time.time()

    sample_file_dict, peSfx, file_fullext = get_read_files_per_sample(
        raw, sampleNames, isPe
    )
    samples = sorted(list(sample_file_dict.keys()))
    logger.info(f"Samples to process: {samples}")

    if not out.is_dir():
        out.mkdir(exist_ok=True)

    genomeBowtie2Idx = buildBowtie2idx(genomes, out=out / "genomeIdx")

    for i, (s, fps) in enumerate(sample_file_dict.items()):
        logger.info(f"Processing {i+1}/{len(samples)}: {s}")

        # prepare align arguments
        if isPe:
            assert len(peSfx) == 2
            samples1: list[Path] = []
            samples2: list[Path] = []
            for fp in fps:
                suffex_len_with_pe: int = len(file_fullext) + len(peSfx[0])
                for j in range(suffex_len_with_pe, len(fp.name)):
                    if any(sfx in fp.name[-j:] for sfx in peSfx):
                        suffex_len_with_pe = j
                        break
                if peSfx[0] in fp.name[-suffex_len_with_pe:]:
                    samples1.append(fp)
                elif peSfx[1] in fp.name[-suffex_len_with_pe:]:
                    samples2.append(fp)
                else:
                    raise ValueError(f"File {fp} not bound to PE {peSfx}")
            assert len(samples1) == len(samples2) and len(samples1) > 0
            runBowtie2(
                genomeBowtie2Idx,
                outPut=out,
                peFiles1=samples1,
                peFiles2=samples2,
                sample=s,
                ncpu=ncpu,
                dryRun=dryRun,
            )
        else:
            runBowtie2(
                genomeBowtie2Idx,
                outPut=out,
                unpairedFiles=fps,
                sample=s,
                ncpu=ncpu,
                dryRun=dryRun,
            )

    logger.info(f"All done, time elapsed {timeDiffStr(tInit)}")
    logger.info("=" * 20 + getTimeStr() + "=" * 20 + "\n" * 2)


def imputePeSuffix(rawFileNames: list[str], peSfx: list[str] = []):
    if len(peSfx) != 0:
        return peSfx
    assert len(rawFileNames) % 2 == 0 and len(rawFileNames) != 0, (
        "Pair end reads should be in pairs, however"
        f" {len(rawFileNames)} found."
    )
    peFound = False
    for i in range(1, 12):
        s: set[str] = set(fn[-i:] for fn in rawFileNames)
        if len(s) == 2:
            peSfx = sorted(s)
            peFound = True
        if len(s) != 2 and peFound:
            break
    if peFound:
        logger.info(f"Found pairend suffix {peSfx}")
        return peSfx
    raise ValueError(f"pair end suffix not found in {rawFileNames}")


def main():
    parser = arg_setup()
    args = parser.parse_args()

    # Process loggers
    log_formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    log_file = args.out / "align.log"
    log_file_handler = logging.FileHandler(log_file)
    log_file_handler.setFormatter(log_formatter)
    logger.addHandler(log_file_handler)
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_formatter)
    logger.addHandler(console_handler)
    logger.setLevel(logging.DEBUG)

    logger.debug(args)
    align_multiple_raw_bowtie2(
        raw=args.raw,
        genomes=args.genome,
        out=args.out,
        isPe=args.isPe,
        sampleNames=args.sampleNames,
        ncpu=args.ncpu,
        dryRun=args.dryRun,
    )

    logger.handlers.clear()


if __name__ == "__main__":
    main()
