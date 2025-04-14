import os
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Literal

from pyBioinfo_modules.wrappers._environment_settings import (
    CONDAEXE, PROKKA_ENV, SHELL, withActivateEnvCmd)


def runProkka(
    fastaPath: Path,
    gcode: int = 11,
    gram: Literal["pos", "neg"] = "pos",
    center: str | None = None,
    genus: str | None = None,
    species: str | None = None,
    strain: str | None = None,
    locustag: str | None = None,
    prokkaEnv: Path | None = PROKKA_ENV,
    condaExe: str = CONDAEXE,
    cpu: int = 4,
    shell: Literal["bash", "zsh"] = SHELL,
    output: Path | None = None,
    prefix: str = "prokka",
    dry: bool = False,
    silent: bool = False,
) -> Path:

    if not silent:
        print(f"Running prokka for {fastaPath}")

    unzip = False
    if fastaPath.suffix == ".gz":
        gzip = subprocess.run(
            f"gzip -dkf {fastaPath}".split(" "), capture_output=True
        )
        assert gzip.returncode == 0
        fastaPath = fastaPath.with_suffix("")
        assert fastaPath.exists(), "\n".join(
            [
                f"Unzip file {fastaPath} failed with error message:\n",
                gzip.stderr.decode(),
                gzip.stdout.decode(),
            ]
        )
        unzip = True

    timeStr = datetime.now().strftime(r"%Y%m%d%H%M")
    prefix = "_".join(
        item
        for item in [prefix, genus, species, strain, timeStr]
        if item is not None
    )
    if output is None:
        outdir = fastaPath.parent / (fastaPath.stem + "_prokka")
    else:
        outdir = output
    cmd = (
        "prokka --compliant --addgenes --mincontiglen 200 --rfam"
        + f" --gcode {gcode}"
        + f" --gram {gram}"
        + f" --cpu {cpu}"
        + f" --outdir {outdir}"
        + f" --prefix {prefix}"
        + (f" --centre {center}" if center is not None else "")
        + (f" --genus {genus}" if genus is not None else "")
        + (f" --strain {strain}" if strain is not None else "")
        + (f" --species {species}" if species is not None else "")
        + (f" --locustag {locustag}" if locustag is not None else "")
        + f" {fastaPath}"
    )

    if not silent:
        print(cmd)

    cmd = withActivateEnvCmd(cmd, prokkaEnv, condaExe, shell)
    if dry:
        print(cmd)
    else:
        commandResult = subprocess.run(
            cmd, capture_output=True, shell=True, executable=shell
        )
        if commandResult.returncode != 0:
            print("Failed prokka:")
            print(cmd)
            print(commandResult.stdout.decode())
            print(commandResult.stderr.decode())

    if unzip:
        os.remove(str(fastaPath))

    return outdir.resolve()
