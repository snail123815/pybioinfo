import subprocess
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Literal

from pyBioinfo_modules.wrappers._environment_settings import (
    CONDAEXE, MASH_ENV, SHELL, withActivateEnvCmd)


def mashSketchFiles(
    inputFiles: list[Path],
    output: Path,
    kmer: int,
    sketch: int,
    nthreads: int = 1,
    mashEnv=MASH_ENV,
    condaExe=CONDAEXE,
    shell=SHELL,
    molecule: Literal["DNA", "protein"] = "DNA",
) -> Path:
    """
    Calculates the distance between the query fasta files
    stored in the sketch file by using mash.
    """
    fileList = NamedTemporaryFile()
    with open(fileList.name, "w") as fl:
        for f in inputFiles:
            fl.write(f'{f.resolve().relative_to(Path(".").resolve())}\n')
    cmd = f"mash sketch -o {output} -k {kmer} -p {nthreads} -s {sketch}"
    cmd += " -a" if molecule == "protein" else ""
    cmd += f" -l {fileList.name}"
    mashSketchRun = subprocess.run(
        withActivateEnvCmd(cmd, mashEnv, condaExe, shell),
        shell=True,
        capture_output=True,
        executable=shell,
    )
    assert mashSketchRun.returncode == 0, (
        withActivateEnvCmd(cmd, mashEnv, condaExe, shell)
        + (mashSketchRun.stdout + mashSketchRun.stderr).decode()
    )
    fileList.close()
    outputMsh = Path(str(output) + ".msh")
    assert outputMsh.is_file()
    return outputMsh


def mashDistance(
    inputMsh: Path,
    outputFile: Path,
    nthreads: int = 1,
    mashEnv=MASH_ENV,
    condaExe=CONDAEXE,
    shell=SHELL,
) -> Path:
    """
        Calculates the distance between the query fasta files
        stored in the sketch file by using mash.
        Parameters
        ----------
        outdir
            string, the path to output directory
        returns
            Path for the output file
       The output fields are
    reference-ID  query-ID     distance  p-value  shared-hashes
    genome1.fna   genome3.fna  0         0        1000/1000
    genome2.fna   genome3.fna  0.022276  0        456/1000
        ----------
    """
    cmd = f"mash dist -p {nthreads} {inputMsh} {inputMsh} > {outputFile}"
    mashDistRun = subprocess.run(
        withActivateEnvCmd(cmd, mashEnv, condaExe, shell),
        capture_output=True,
        shell=True,
        executable=shell,
    )
    assert mashDistRun.returncode == 0, (
        withActivateEnvCmd(cmd, mashEnv, condaExe, shell)
        + (mashDistRun.stdout + mashDistRun.stderr).decode()
    )
    assert outputFile.exists
    return outputFile
