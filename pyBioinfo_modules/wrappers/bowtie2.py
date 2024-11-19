from pathlib import Path
import logging
from tempfile import NamedTemporaryFile
from Bio import SeqIO
from typing import IO
from pyBioinfo_modules.basic.basic import getTimeStr
from pyBioinfo_modules.wrappers._environment_settings \
    import SHORTREADS_ENV, SHELL, withActivateEnvCmd
import subprocess
import time
from pyBioinfo_modules.basic.basic import getTimeStr, timeDiffStr

def buildBowtie2idx(fs: list[Path], out: Path, name=None) -> Path:
    fs = [f.resolve() for f in fs]
    out = out.resolve()

    # check if idx file exists
    genomePath = fs[0].parent
    genomeName = fs[0].stem
    logging.info(f'genome name: {genomeName}')

    out.mkdir(exist_ok=True)
    bt2_base = ('_'.join(f.stem for f in fs) if name is None else name)
    outIdxForUse = out / bt2_base

    if any((Path(str(outIdxForUse) + '.1.bt2').is_file(),
            Path(str(outIdxForUse) + '.1.bt21').is_file())):
        logging.info(
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
            logging.error(f'Unexpected error {err=}, {type(err)=}')
            raise
        fs = newFps

    logging.info('-' * 20 + 'Indexing genome ' + getTimeStr() + '-' * 20)
    cmdList = [
        'bowtie2-build',
        ','.join(str(f) for f in fs),
        str(out / bt2_base),
    ]
    logging.info(' '.join(cmdList))
    cmd = withActivateEnvCmd(' '.join(cmdList), SHORTREADS_ENV)
    result = subprocess.run(cmd, capture_output=True,
                            shell=True, executable=SHELL)
    (f.close() for f in tempFiles)
    if result.returncode != 0 or not any(
        (Path(str(outIdxForUse) + '.1.bt2').is_file(),
         Path(str(outIdxForUse) + '.1.bt21').is_file())):
        logging.info('stderr: ' + result.stderr.decode())
        logging.info('stdout: ' + result.stdout.decode())
        logging.info(
            '-' *
            20 +
            'Error Indexing genome' +
            getTimeStr() +
            '-' *
            20)
        raise Exception
    logging.info('-' * 20 + 'DONE Indexing genome' + getTimeStr() + '-' * 20)
    logging.info('\n' * 2)

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
            logging.info(f'Found finished bam file {str(target)}')
        else:
            logging.info(f'Found finished flag but not bam file.')
            raise FileNotFoundError(str(target))
    else:
        cmdList.extend([
            '|', 'samtools', 'view', '-bS', '-@', str(toBamNcpu),
            "|", "samtools", "sort", '-@', str(
                toBamNcpu), "--write-index", '-o', str(target)
        ])
        logging.info(' '.join(cmdList))

        # Start running both
        cmd = withActivateEnvCmd(' '.join(cmdList), SHORTREADS_ENV)
        if dryRun:
            return
        result = subprocess.run(cmd, shell=True,
                                capture_output=True, executable=SHELL)
        if result.returncode != 0:
            logging.info('stderr: ' + result.stderr.decode())
            logging.info('stdout: ' + result.stdout.decode())
        # stderr has logging.info info from bowtie2
        logging.info(result.stderr.decode())
        logging.info(f'Finished in {timeDiffStr(ts)}\n')
        targetFinishedFlag.touch()
