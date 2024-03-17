import subprocess
import sys
from pathlib import Path
import shutil

from pyBioinfo_modules.wrappers._environment_settings \
    import CONDAEXE, SHELL, BIGSCAPE_ENV, PFAM_DB, withActivateEnvCmd

try:
    if BIGSCAPE_ENV is None:
        bigscapeExe = Path([p for p in [
            shutil.which('bigscape'), shutil.which('bigscape.py')
        ] if p is not None][0]).resolve()
    else:
        bigscapeExe = next((BIGSCAPE_ENV / 'bin').glob('bigscape*')).resolve()
except (IndexError, StopIteration):
    raise FileNotFoundError('bigscape executable not found.')


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
