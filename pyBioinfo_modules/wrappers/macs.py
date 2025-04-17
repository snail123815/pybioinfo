import logging
import os
import shutil
import subprocess
from pathlib import Path

from pyBioinfo_modules.wrappers._environment_settings import (
    CONDAEXE, MACS_ENV, MACS_PROGRAM, SHELL, withActivateEnvCmd)

logger = logging.getLogger(__name__)


def macs_predictd(experimentDict: Path, outputDir: Path, gsize):
    expfiles = [experimentDict[e]["exp"] for e in experimentDict]
    # exps = [ef.stem for ef in expfiles]
    fragsizes = []
    for e in experimentDict:
        inputFile = Path(experimentDict[e]["exp"])
        sample = inputFile.stem
        rfileName = f"{sample}_preidctd.r"
        rfile = outputDir / rfileName
        # if rfile.exists():
        #     continue
        argsPredictd = [
            MACS_PROGRAM,
            "predictd",
            "-i",
            str(inputFile),
            "--gsize",
            gsize,
            "--rfile",
            rfileName,
            "--outdir",
            str(outputDir),
        ]
        print("Running:")
        print(" ".join(argsPredictd))
        p1 = subprocess.Popen(
            withActivateEnvCmd(
                argsPredictd, condaEnv=MACS_ENV, condaExe=CONDAEXE, shell=SHELL
            ),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        while True:
            output = p1.stdout.readline()
            if len(output) == 0 and p1.poll() is not None:
                break
            if output:
                print(output.decode("utf-8"), end="")
        if p1.returncode != 0:
            print(f"Error running {MACS_PROGRAM} predictd.")
            print(f"Return code: {p1.returncode}")
            raise Exception(f"{MACS_PROGRAM} predictd failed")

        # Test if "Rscript" is in PATH
        if shutil.which("Rscript") is not None:
            print("\nFinished, plotting with R script...")
            argsRscript = ["Rscript", rfileName]
            p2 = subprocess.Popen(
                argsRscript,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                cwd=outputDir,
            )
            for line in p2.stdout:
                print(line.decode("utf-8"), end="")
        else:
            print("Rscript not found in PATH, skip plotting")
        if rfile.exists():
            print(f"Fragment size prediction for {sample} is done")
            with open(rfile, "r") as rscript:
                for line in rscript.readlines():
                    if "alt lag(s) : " in line:
                        num = int(line.split(" : ")[1].split("'")[0])
                        fragsizes.append(num)
                        break
            print(f"fragement size predicition for {sample} is {num}")
        else:
            print(f"Fragment size prediction for {sample} failed")
    if len(fragsizes) > 0:
        meansize = int(sum(fragsizes) / len(fragsizes))
        print(f"Average fragment size predicted is {meansize}")
        return meansize
    else:
        return None


# predictd


def parse_comparison_file(compFile: Path, bamPath: Path) -> dict[str:Path]:
    """
    name    ctr exp
    G24 G24C_G24C.sam   G24E_G24E.sam
    G48 G48C_G48C.sam   G48E_G48E.sam
    M24 M24C_M24C.sam   M24E_M24E.sam
    M48 M48C_M48C.sam   M48E_M48E.sam

    experimentDict = {
        'exp1': {'exp':'filePathA', 'ctr':'filePathB'}
        'exp2': {'exp':'filePathC', 'ctr':'filePathD'}
        }


    """
    experimentDict = {}
    with open(compFile, "r") as f:
        for i, l in enumerate(f.readlines()):
            ls = l.strip().split("\t")
            if i == 0:
                assert all(x in ls for x in ["name", "ctr", "exp"]), (
                    "The first line of the comparisons file should have the "
                    "headers: 'name', 'ctr', 'exp' "
                    "but it has: \n"
                    f"{ls}"
                )
                ctri = ls.index("ctr")
                expi = ls.index("exp")
                continue
            experimentDict[ls[0]] = {}
            for n, i in zip(["ctr", "exp"], [ctri, expi]):
                f = bamPath / ls[i]
                if not f.exists():
                    if f.suffix == "":
                        logger.warning(f"Try to add .bam to the file name: {f}")
                        f = f.with_suffix(".bam")
                    elif f.suffix != ".bam":
                        # it is possible that the file name contains '.'
                        logger.warning(
                            f"Target file from comparison file is"
                            f"not of .bam format: {f}"
                        )
                        f = bamPath / ls[i] + ".bam"
                if not f.exists():
                    raise FileNotFoundError(
                        f"Target file from comparison file does not exist: {f}"
                    )
                experimentDict[ls[0]][n] = f
    return experimentDict


def macs_callPeak(
    experimentDict,
    outputDir,
    gsize,
    isPe: bool = False,
    fragsize: None | int = None,
    fdr="1e-20",
):
    for e in experimentDict:
        ctr = experimentDict[e]["ctr"]
        exp = experimentDict[e]["exp"]
        singleExpOutputDir = os.path.join(outputDir, e)
        if not os.path.exists(outputDir):
            os.mkdir(outputDir)

        args = [
            MACS_PROGRAM,
            "callpeak",
            "-t",
            str(exp),
            "-c",
            str(ctr),
            "-n",
            str(e),
            "--outdir",
            str(singleExpOutputDir),
            "-f",
            "BAM",
            "--gsize",
            str(gsize),
            "-q",
            str(fdr),
            "--call-summits",
            "-B",
            # -B Save extended fragment pileup, and local lambda
            # tracks (two files) at every bp into a bedGraph file
            "--keep-dup",
            "all",
        ]
        if fragsize is not None:
            args.extend(
                [
                    "--extsize",
                    str(fragsize),
                    "--nomodel",
                ]
            )
        if isPe:
            if "--nomodel" not in args:
                args.extend(["--nomodel"])
        print("Running:")
        print(" ".join(args))
        p = subprocess.Popen(
            args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        # Check return code and process output
        for line in p.stdout:
            print(line.decode("utf-8"), end="")
        p.wait()
        if p.returncode == 0:
            print("=" * 80)
            print()
            continue

        print(f"Error running {MACS_PROGRAM} callpeak with default settings.")
        print(f"Return code: {p.returncode}")
        if p.stderr:
            print("Error message:")
            for line in p.stderr:
                print(line.decode("utf-8"), end="")
        newargs = args.copy()
        assert (
            "--nomodel" not in newargs
        ), "--nomodel already in args, unknow reason for failure"
        newargs.extend(["--nomodel", "--extsize", "200"])
        print("Trying again with extsize 200")
        p = subprocess.Popen(
            newargs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        # Check return code and process output
        for line in p.stdout:
            print(line.decode("utf-8"), end="")
        p.wait()
        if p.returncode == 0:
            print("=" * 80)
            print()
        else:
            raise Exception(f"{MACS_PROGRAM} callpeak failed with extsize 200")
