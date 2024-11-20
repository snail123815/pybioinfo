import logging
import subprocess
import time
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import IO

from Bio import SeqIO

from pyBioinfo_modules.basic.basic import getTimeStr, timeDiffStr
from pyBioinfo_modules.basic.parse_raw_read_dir import \
    get_read_files_per_sample
from pyBioinfo_modules.wrappers._environment_settings import (
    SHELL, SHORTREADS_ENV, withActivateEnvCmd)

logger = logging.getLogger(__name__)

def buildBowtie2idx(fs: list[Path], out: Path, name=None) -> Path:
    fs = [f.resolve() for f in fs]
    out = out.resolve()

    # check if idx file exists
    genomePath = fs[0].parent
    genomeName = fs[0].stem
    logger.info(f'genome name: {genomeName}')

    out.mkdir(exist_ok=True)
    bt2_base = ('_'.join(f.stem for f in fs) if name is None else name)
    outIdxForUse = out / bt2_base

    if any((Path(str(outIdxForUse) + '.1.bt2').is_file(),
            Path(str(outIdxForUse) + '.1.bt21').is_file())):
        logger.info(
            'Found index file, will not make new ones.\n'
            f'{str(list(out.glob(bt2_base + "*"))[0])}')
        return outIdxForUse

    # convert gbk to fa
    tempFiles: list[IO] = []
    if fs[0].suffix not in ['.fa', '.fasta', '.fna', '.fsa']:
        # try gbk
        newF = NamedTemporaryFile()
        newFps: list[Path] = []
        try:
            for f in fs:
                for s in SeqIO.parse(f, 'genbank'):
                    SeqIO.write(s, newF.name, 'fasta')
                    newFps.append(Path(newF.name))
                    tempFiles.append(newF)  # for closing these files later
        except Exception as err:
            logger.error(f'Unexpected error {err=}, {type(err)=}')
            raise
        fs = newFps

    logger.info('-' * 20 + 'Indexing genome ' + getTimeStr() + '-' * 20)
    cmdList = [
        'bowtie2-build',
        ','.join(str(f) for f in fs),
        str(out / bt2_base),
    ]
    logger.info(' '.join(cmdList))
    cmd = withActivateEnvCmd(' '.join(cmdList), SHORTREADS_ENV)
    result = subprocess.run(cmd, capture_output=True,
                            shell=True, executable=SHELL)
    (f.close() for f in tempFiles)
    if result.returncode != 0 or not any(
        (Path(str(outIdxForUse) + '.1.bt2').is_file(),
         Path(str(outIdxForUse) + '.1.bt21').is_file())):
        logger.info('stderr: ' + result.stderr.decode())
        logger.info('stdout: ' + result.stdout.decode())
        logger.info(
            '-' *
            20 +
            'Error Indexing genome' +
            getTimeStr() +
            '-' *
            20)
        raise Exception
    logger.info('-' * 20 + 'DONE Indexing genome' + getTimeStr() + '-' * 20)
    logger.info('\n' * 2)

    return outIdxForUse


def runBowtie2(
    genomeBowtie2Idx: Path,
    outPut: Path,
    peFiles1:  list[Path] = [],
    peFiles2:  list[Path] = [],
    unpairedFiles:  list[Path] = [],
    sample: str = '',
    ncpu: int = 2,
    dryRun: bool = False
):
    ts = time.time()
    allFiles = peFiles1 + peFiles2 + unpairedFiles
    assert len(allFiles) > 0, 'Files are needed.'
    # prepare align arguments
    cmdList = ['bowtie2', '-x', str(genomeBowtie2Idx), '-p', str(ncpu)]
    if any(len(sps) > 0 for sps in [peFiles1, peFiles2]):
        assert len(peFiles1) == len(peFiles2)
        cmdList.extend([
            '-1', ','.join(str(s) for s in peFiles1),
            '-2', ','.join(str(s) for s in peFiles2)
        ])
    if len(unpairedFiles) > 0:
        cmdList.extend([
            '-U', ','.join(str(s) for s in unpairedFiles)
        ])

    # prepare convert to bam arguments
    if sample == '':
        # inpute sample name
        sample = allFiles[0].stem[0]
        for i in range(len(allFiles[0].stem)-1):
            if len(set([f.stem[:i+2] for f in allFiles])) == 1:
                sample = allFiles[0].stem[:i+2]
            else:
                break
    target = outPut / f'{sample.strip("_")}.bam'
    targetFinishedFlag = outPut / f'{sample.strip("_")}.done'

    toBamNcpu = max(ncpu // 8, 1)
    if targetFinishedFlag.is_file():
        if target.is_file():
            logger.info(f'Found finished bam file {str(target)}')
        else:
            logger.info(f'Found finished flag but not bam file.')
            raise FileNotFoundError(str(target))
    else:
        cmdList.extend([
            '|', 'samtools', 'view', '-bS', '-@', str(toBamNcpu),
            "|", "samtools", "sort", '-@', str(
                toBamNcpu), "--write-index", '-o', str(target)
        ])
        logger.info(' '.join(cmdList))

        # Start running both
        cmd = withActivateEnvCmd(' '.join(cmdList), SHORTREADS_ENV)
        if dryRun:
            return
        result = subprocess.run(cmd, shell=True,
                                capture_output=True, executable=SHELL)
        if result.returncode != 0:
            logger.info('stderr: ' + result.stderr.decode())
            logger.info('stdout: ' + result.stdout.decode())
        # stderr has logging.info info from bowtie2
        logger.info(result.stderr.decode())
        logger.info(f'Finished in {timeDiffStr(ts)}\n')
        targetFinishedFlag.touch()


def multiple_raw_align_bowtie2(
    raw: list[Path],
    sampleNames: list[str] | None,
    out: Path,
    isPe: bool,
    genomes: list[Path],
    ncpu: int,
    dryRun: bool,
):
    # Start logic
    logger.addHandler(logging.FileHandler(out / "align.log"))
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
    # Remove logger handler with path "out / align.log"
    for handler in logger.handlers[:]:
        if isinstance(
            handler, logging.FileHandler
        ) and handler.baseFilename == str(out / "align.log"):
            logger.removeHandler(handler)
            handler.close()