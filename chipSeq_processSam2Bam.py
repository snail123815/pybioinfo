# This file is licensed under the MIT License

import argparse
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pathSam", help="path to sam files")
parser.add_argument("-o", "--output", help="path to sorted bam files")
parser.add_argument("-t", "--threads", help="number of threads", type=int)

args = parser.parse_args()

pathSam = args.pathSam
pathOut = args.output
ncpu = str(args.threads)

if not os.path.isdir(pathOut):
    os.makedirs(pathOut)

fileList = []
for path, subdirs, files in os.walk(pathSam):
    for f in files:
        if f.endswith(".sam"):
            filePath = os.path.join(path, f)
            fileList.append(filePath)

for f in fileList:
    name = os.path.split(f)[1]
    bam = os.path.join(pathOut, f"{name[:-4]}.bam")
    bai = os.path.join(pathOut, f"{name[:-4]}.bai")

    cmd = [
        "samtools",
        "sort",
        "-m",
        "20G",
        "-o",
        bam,
        "-O",
        "bam",
        "-@",
        ncpu,
        f,
    ]
    print("Running...\n", " ".join(cmd))
    p = subprocess.run(cmd, capture_output=True)
    stdout = p.stdout.decode("utf-8")
    stderr = p.stderr.decode("utf-8")
    if len(stdout) > 0:
        print("stdout:\n", stdout)
    if len(stderr) > 0:
        print("stderr:\n", stderr)
    print("output file: ", bam, "\n")

    cmd = ["samtools", "index", bam, bai]
    print("Running...\n", " ".join(cmd))
    stdout = p.stdout.decode("utf-8")
    stderr = p.stderr.decode("utf-8")
    if len(stdout) > 0:
        print("stdout:\n", stdout)
    if len(stderr) > 0:
        print("stderr:\n", stderr)
    print("\n")
    print("output file: ", bai, "\n")
