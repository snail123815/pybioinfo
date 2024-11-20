import subprocess
import sys
import logging
from pathlib import Path
import shutil

from pyBioinfo_modules.wrappers._environment_settings \
    import CONDAEXE, SHELL, BIGSCAPE_ENV, PFAM_DB, withActivateEnvCmd


logger = logging.getLogger(__name__)

try:
    if BIGSCAPE_ENV is None:
        bigscapeExe = Path([p for p in [
            shutil.which('bigscape'), shutil.which('bigscape.py')
        ] if p is not None][0]).resolve()
    else:
        bigscapeExe = next((BIGSCAPE_ENV / 'bin').glob('bigscape*')).resolve()
except (IndexError, StopIteration):
    raise FileNotFoundError('bigscape executable not found.')


def get_bigscape_version(
    bigscape_exe=bigscapeExe,
    condaexe=CONDAEXE,
    bigscape_env=BIGSCAPE_ENV,
    shell=SHELL,
) -> str:
    logger.info("Getting BiG-SCAPE version")
    logger.debug(f"bigscape_exe: {bigscape_exe}")
    logger.debug(f"condaexe: {condaexe}")
    logger.debug(f"bigscape_env: {bigscape_env}")
    logger.debug(f"shell: {shell}")
    cmd = withActivateEnvCmd(
        f"{bigscape_exe} --version", bigscape_env, condaexe, shell
    )
    logger.debug(f"cmd: {cmd}")
    version_run = subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        executable=shell,
    )
    version_split = version_run.stdout.decode().strip().split(" ")
    version = " ".join(version_split)
    date = None
    if version_split[0].lower() in ["BiG-SCAPE".lower(), "bigscape"]:
        version = version_split[1]
    if len(version_split) == 3:
        date = version_split[2]
    logger.info(f"BiG-SCAPE version: {version}, {date}")
    return version, date


def runBigscape(
    inputPath: Path,
    outputPath: Path,
    cpus: int = 4,
    cutoffs: list[float] = [0.2, ],
    silent: bool = True,
    bigscapeEnv=BIGSCAPE_ENV,
    pfamDb=PFAM_DB,
    condaExe=CONDAEXE,
    verbose=False,
    shell=SHELL
) -> Path:
    if not outputPath.is_dir():
        outputPath.mkdir(parents=True, exist_ok=False)
    cmd = f'{bigscapeExe} --mode auto -c {cpus}'
    cmd += f' --cutoffs {" ".join([str(c) for c in cutoffs])}'
    if pfamDb is not None:
        cmd += f' --pfam_dir {pfamDb}'
    if verbose:
        cmd += ' --verbose'
    cmd += f' -i {inputPath}'
    cmd += f' -o {outputPath}'
    # print(withActivateEnvCmd(cmd, bigscapeEnv, condaExe, shell))
    try:
        if silent:
            bigscapeRun = subprocess.run(
                withActivateEnvCmd(cmd, bigscapeEnv, condaExe, shell),
                shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, executable=shell
            )
        else:
            bigscapeRun = subprocess.run(
                withActivateEnvCmd(cmd, bigscapeEnv, condaExe, shell),
                shell=True, stdout=sys.stdout, stderr=sys.stderr, executable=shell
            )

    except subprocess.CalledProcessError:
        print(bigscapeRun.stdout.decode())
        print(bigscapeRun.stderr.decode())
    return outputPath
