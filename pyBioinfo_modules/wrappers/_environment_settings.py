# TODO: write a script change this file
from pathlib import Path
from typing import Literal

SHELL: Literal["bash", "zsh"] = "zsh"
CONDAEXE: Literal["conda", "mamba", "micromamba"] = "micromamba"

MACS_ENV: Path | None = Path.home() / "genvs/macs3"
MACS_PROGRAM = "macs3"

ANTISMASH_ENV: Path | None = Path.home() / "genvs/quasan"
BUSCO_ENV: Path | None = Path.home() / "genvs/quasan"
PROKKA_ENV: Path | None = Path.home() / "genvs/quasan"
MASH_ENV: Path | None = Path.home() / "genvs/phylophlan"
BIGSCAPE_ENV: Path | None = Path.home() / "genvs/bigscape"
SHORTREADS_ENV: Path | None = Path.home() / "genvs/shortReads"

PFAM_DB: Path | None = None
try:
    PFAM_DB = sorted(
        [
            p
            for p in (Path.home() / "dbMisc/antismash_databases/pfam").iterdir()
            if p.is_dir()
        ]
    )[-1]
except FileNotFoundError:
    pass


def withActivateEnvCmd(
    cmd: str, condaEnv: Path | None = None, condaExe=CONDAEXE, shell=SHELL
) -> str:
    """
    Join commands that activate your environment and your command.

    Running this with popen needs shell=True.

    Parameters:
    cmd (str): The command to run after activating the environment.
    condaEnv (Path | None): The path to the conda environment to activate.
    condaExe (str): The conda executable to use (default is CONDAEXE).
    shell (str): The shell to use (default is SHELL).

    Returns:
    str: The combined command to activate the environment and run the given command.
    """
    if condaEnv is not None:
        if condaExe == "micromamba":
            activateEnvCmd = f'eval "$(micromamba shell hook --shell={shell})"'
        else:
            activateEnvCmd = f'eval "$({condaExe} shell.{shell} hook)"'
        activateEnvCmd += f" && micromamba activate {condaEnv}"
        if isinstance(cmd, list):
            cmd = " ".join(cmd)
        cmd = activateEnvCmd + " && " + cmd
    return cmd
