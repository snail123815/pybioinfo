import logging
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Tuple

from pyBioinfo_modules.wrappers._environment_settings import (
    BIGSCAPE_ENV, CONDAEXE, PFAM_DB, SHELL, withActivateEnvCmd)

logger = logging.getLogger(__name__)


class BigscapeWrapper:
    def __init__(
        self,
        conda_exe=CONDAEXE,
        bigscape_env=BIGSCAPE_ENV,
        pfam_db=PFAM_DB,
        shell=SHELL,
        bigscape_exe=None,  # Allow injecting a path for testing
    ):
        self.conda_exe = conda_exe
        self.bigscape_env = bigscape_env
        self.pfam_db = pfam_db
        self.shell = shell
        self._exe = bigscape_exe  # When None, will be resolved

        logger.info("BiG-SCAPE wrapper initialized")
        logger.debug(f"bigscape_exe: {self._exe}")
        logger.debug(f"condaexe: {self.conda_exe}")
        logger.debug(f"bigscape_env: {self.bigscape_env}")
        logger.debug(f"shell: {self.shell}")

    def get_bigscape_exe(self) -> Path:
        """Set path to BigScape executable, caching the result."""
        if self._exe is not None:
            return self._exe

        try:
            if self.bigscape_env is None:
                self._exe = Path(
                    [
                        p
                        for p in [
                            shutil.which("bigscape"),
                            shutil.which("bigscape.py"),
                        ]
                        if p is not None
                    ][0]
                ).resolve()
            else:
                self._exe = next(
                    (BIGSCAPE_ENV / "bin").glob("bigscape*")
                ).resolve()
        except (IndexError, StopIteration):
            raise FileNotFoundError("bigscape executable not found.")

        return self._exe

    def get_version(self) -> Tuple[str, Optional[str]]:
        """Get BigScape version as (version, date)."""
        logger.info("Getting BiG-SCAPE version")
        self.get_bigscape_exe()
        cmd = withActivateEnvCmd(
            f"{self._exe} --version",
            self.bigscape_env,
            self.conda_exe,
            self.shell,
        )

        version_run = subprocess.run(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            executable=self.shell,
        )

        version_split = version_run.stdout.decode().strip().split(" ")
        version = " ".join(version_split)
        date = None

        if version_split[0].lower() in ["big-scape", "bigscape"]:
            version = version_split[1]

        if len(version_split) == 3:
            date = version_split[2]

        logger.info(f"BiG-SCAPE version: {version}, {date}")
        return version, date

    def run_bigscape(
        self,
        inputPath: Path,
        outputPath: Path,
        cpus: int = 4,
        cutoffs: list[float] = [
            0.2,
        ],
        verbose=False,
    ) -> Path:
        if not outputPath.is_dir():
            outputPath.mkdir(parents=True, exist_ok=False)
        cmd = f"{self._exe} --mode auto -c {cpus}"
        cmd += f' --cutoffs {" ".join([str(c) for c in cutoffs])}'
        if self.pfam_db is not None:
            cmd += f" --pfam_dir {self.pfam_db}"
        if verbose:
            cmd += " --verbose"
        cmd += f" -i {inputPath}"
        cmd += f" -o {outputPath}"

        bigscapeRun = subprocess.run(
            withActivateEnvCmd(
                cmd, self.bigscape_env, self.conda_exe, self.shell
            ),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            executable=self.shell,
        )
        if bigscapeRun.returncode != 0:
            logger.warning(bigscapeRun.stdout.decode())
        else:
            logger.info(bigscapeRun.stdout.decode())
            logger.info(f"BiG-SCAPE finished successfully")
        return outputPath
