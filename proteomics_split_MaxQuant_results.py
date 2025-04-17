#######
# For each dir, split the files based on:
# rawD = dname+'_raw'
# mqD = dname+'_mqAnalysis'
# mqRes = dname+'_mqResult'


import argparse
import os
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("-p", help="Root path")
parser.add_argument("--dry", action="store_true", help="dry run")
parser.add_argument(
    "--rmMQinter", action="store_true", help="remove mq intermediate result"
)

args = parser.parse_args()

root = args.p
dryRun = args.dry
rmInter = args.rmMQinter
ru = not dryRun
root = os.path.realpath(root)


def fileop(r, n, t, op, ru, mkt=True):
    op = op[0].upper() + op[1:].lower()
    p = os.path.join(r, n)
    if not mkt:
        if not os.path.isdir(p):
            raise ValueError(f"{p} needs to be a dir!")
        np = os.path.join(r, t)
    else:
        if ru:
            os.makedirs(os.path.join(r, t), exist_ok=True)
        np = os.path.join(r, t, n)
    print(f"\n{op}")
    print(p)
    if not op in ["Delete", "Remove"]:
        print(np)
    if op == "Copy":
        if ru:
            if os.path.isdir(p):
                shutil.copytree(p, np)
            else:
                shutil.copy(p, np)
    elif op in ["Delete", "Remove"]:
        if ru:
            if os.path.isdir(p):
                shutil.rmtree(p)
            else:
                os.remove(p)
    elif op == "Move":
        if ru:
            shutil.move(p, np)


for r, ds, fs in os.walk(root):
    dname = os.path.split(r)[1]
    if dname.endswith("_raw") or dname == "raw":
        continue
    # Check if .raw file in current folder:
    rawFs = []
    indexFs = []
    indexDs = []
    fastaFs = []
    paraF = "mqpar.xml"
    metaFs = []
    combD = "combined"
    resD = os.path.join(combD, "txt")
    for f in fs:
        n, ext = os.path.splitext(f)
        ext = ext.lower()
        if ext == ".raw":
            rawFs.append(f)
            indexDs.append(n)
        elif ext == ".index":
            indexFs.append(f)
        elif ext in [".fa", ".fasta"]:
            fastaFs.append(f)
        elif ext in [
            ".md",
            ".txt",
            ".docx",
            "doc",
            "xls",
            ".xlsx",
            ".sld",
            ".lyt",
        ]:
            metaFs.append(f)
    rawD = dname + "_raw"
    mqD = dname + "_mqAnalysis"
    mqRes = dname + "_mqResult"

    indexDs = [d for d in indexDs if os.path.isdir(os.path.join(r, d))]

    if len(rawFs) > 0 and len(indexDs) == len(indexFs):
        if os.path.isdir(os.path.join(r, resD)):
            if rmInter:
                op = "Move"
            else:
                op = "Copy"
            fileop(r, resD, mqRes, op, ru, mkt=False)
        # Move raw files to raw folder
        for f in rawFs:
            fileop(r, f, rawD, "Move", ru)
        # Metadata
        for f in metaFs + fastaFs + [paraF]:
            # Copy metadata to rawfolder
            fileop(r, f, rawD, "Copy", ru)
            # Copy metadata to mq folder
            if not rmInter:
                fileop(r, f, mqD, "Copy", ru)
        # Move mqdata to mq folder
        for fd in indexDs + [combD] + indexFs + [paraF]:
            if rmInter:
                if not fd == paraF:
                    op = "Delete"
                else:
                    op = "Move"
            else:
                op = "Move"
            fileop(r, fd, mqD, op, ru)
