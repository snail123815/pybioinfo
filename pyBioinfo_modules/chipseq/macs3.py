import os
import subprocess


def predictd(experimentDict, outputDir, gsize):
    expfiles = [experimentDict[e]["exp"] for e in experimentDict]
    exps = [os.path.splitext(ef)[0] for ef in expfiles]
    fragsizes = []
    for e in experimentDict:
        inputFile = experimentDict[e]["exp"]
        sample = os.path.splitext(os.path.split(inputFile)[1])[0]
        rfileName = f"{sample}_preidctd.r"
        rfile = os.path.join(outputDir, rfileName)
        if os.path.isfile(rfile):
            continue
        argsPredictd = [
            "macs3",
            "predictd",
            "-i",
            inputFile,
            "--gsize",
            gsize,
            "--rfile",
            rfileName,
            "--outdir",
            outputDir,
        ]
        print("Running...")
        print(" ".join(argsPredictd))
        p1 = subprocess.Popen(
            argsPredictd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        for line in p1.stdout:
            print(line.decode("utf-8"), end="")
        try:
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
            with open(rfile, "r") as rscript:
                for line in rscript.readlines():
                    if "alt lag(s) : " in line:
                        num = int(line.split(" : ")[1].split("'")[0])
                        fragsizes.append(num)
                        break
            print(f"fragement size predicition for {sample} is {num}")
        except BaseException:
            break
    try:
        meansize = sum(fragsizes) / len(fragsizes)
        print(f"Average fragment size is {meansize}")
    except BaseException:
        meansize = None
    return meansize


# predictd


def readComps(compFile, bamPath):
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
                ctri = ls.index("ctr")
                expi = ls.index("exp")
                continue
            ctr = os.path.join(bamPath, ls[ctri])
            exp = os.path.join(bamPath, ls[expi])
            experimentDict[ls[0]] = {"ctr": ctr, "exp": exp}
    return experimentDict


def callPeak(
    experimentDict, outputDir, gsize, isPe=False, fragsize=None, fdr="1e-20"
):
    for e in experimentDict:
        ctr = experimentDict[e]["ctr"]
        exp = experimentDict[e]["exp"]
        singleExpOutputDir = os.path.join(outputDir, e)
        if not os.path.exists(outputDir):
            os.mkdir(outputDir)

        args = [
            "macs3",
            "callpeak",
            "-t",
            exp,
            "-c",
            ctr,
            "-n",
            e,
            "--outdir",
            singleExpOutputDir,
            "-f",
            "BAM",
            "--gsize",
            gsize,
            "-q",
            fdr,
            "--call-summits",
            "-B",
            "--keep-dup",
            "all",
        ]
        if not isinstance(fragsize, type(None)):
            args.extend(
                [
                    "--extsize",
                    fragsize,
                ]
            )
        if isPe:
            args.extend(["--nomodel"])
        print("Running...")
        print(" ".join(args))
        p = subprocess.Popen(
            args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        for line in p.stdout:
            print(line.decode("utf-8"), end="")


# callPeak
