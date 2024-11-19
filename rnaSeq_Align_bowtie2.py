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


def main(args):
    # Gether all files, max depth 2
    filePaths: list[Path] = []
    for d in args.raw:
        for n in d.iterdir():
            if n.is_file():
                filePaths.append(n)
            elif n.is_dir():
                filePaths.extend(f for f in n.iterdir() if f.is_file())
    for f in filePaths.copy():
        ext = splitStemSuffixIfCompressed(f)[1]
        if len(ext.split(".")) < 2 or ext.split(".")[1] not in ["fastq", "fq"]:
            filePaths.remove(f)
    assert len(filePaths) > 0, f"Files not found in {args.raw}"
    FILE_PATHS = sorted(filePaths)
    FILE_PARTS = [
        splitStemSuffixIfCompressed(f, fullSuffix=True) for f in FILE_PATHS
    ]
    FILE_NAMES = [fp[0] for fp in FILE_PARTS]
    FILE_FULLEXTS = [fp[1] for fp in FILE_PARTS]
    assert len(set(FILE_FULLEXTS)) == 1, (
        "All files should have same format, "
        f"yet multiple found {set(FILE_FULLEXTS)}."
    )
    assert len(set(FILE_NAMES)) == len(FILE_NAMES), (
        "File names has duplicates"
        ", maybe from different dirs?\n"
        f"{set([fn for fn in FILE_NAMES if FILE_NAMES.count(fn) > 1])}"
    )

    if not args.out.is_dir():
        args.out.mkdir(exist_ok=True)
    
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

    # Start logic
    logger.info("=" * 20 + getTimeStr() + "=" * 20)
    logger.debug(args)
    if args.dryRun:
        logger.debug("=" * 20 + "Dry run, will not execute bowtie2")

    if args.sampleNames is not None:
        samples = args.sampleNames
    else:
        samples = []

    sampleFileDict: dict[str, list[Path]] = {}

    peSfx: list[str] = []
    if args.isPe:
        peSfx = imputePeSuffix(FILE_NAMES)
        if samples == []:
            samples = (
                sorted(list(set(f[: -len(peSfx[0])] for f in FILE_NAMES)))
                if len(samples) == 0
                else samples
            )
            assert (len(FILE_PATHS) / len(samples)) % 2 == 0
        for s in samples:
            fns = [FILE_PATHS[i] for i, fn in enumerate(FILE_NAMES) if s in fn]
            assert len(fns) > 0, f"Did not find files for sample {s}"
            assert len(fns) % 2 == 0, (
                f"Number of files for sample {s} is not paired.\n" f"{fns}"
            )
            sampleFileDict[s] = fns

    else:  # single end reads
        if samples == []:
            for fn, fp in zip(FILE_NAMES, FILE_PATHS):
                sampleFileDict[fn] = [fp]
            samples = FILE_NAMES
        elif len(samples) == 1:
            sampleFileDict[samples[0]] = FILE_PATHS
        else:
            for s in samples:
                fps = [
                    FILE_PATHS[i] for i, fn in enumerate(FILE_NAMES) if s in fn
                ]
                assert len(fps) > 0, f"Did not find files for sample {s}"
                sampleFileDict[s] = fps

    tInit = time.time()

    genomeBowtie2Idx = buildBowtie2idx(args.genome, out=args.out / "genomeIdx")

    logger.info(f"Samples to process: {samples}")
    for i, (s, fps) in enumerate(sampleFileDict.items()):
        logger.info(f"Processing {i+1}/{len(samples)}: {s}")

        # prepare align arguments
        if args.isPe:
            assert len(peSfx) == 2
            samples1: list[Path] = []
            samples2: list[Path] = []
            for fp in fps:
                suffex_len_with_pe: int = len(FILE_FULLEXTS[0]) + len(peSfx[0])
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
                outPut=args.out,
                peFiles1=samples1,
                peFiles2=samples2,
                sample=s,
                ncpu=args.ncpu,
                dryRun=args.dryRun,
            )
        else:
            runBowtie2(
                genomeBowtie2Idx,
                outPut=args.out,
                unpairedFiles=fps,
                sample=s,
                ncpu=args.ncpu,
                dryRun=args.dryRun,
            )

    logger.info(f"All done, time elapsed {timeDiffStr(tInit)}")
    logger.info("=" * 20 + getTimeStr() + "=" * 20 + "\n" * 2)
    logger.handlers.clear()



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


if __name__ == "__main__":
    parser = arg_setup()
    args = parser.parse_args()
    main(args)
