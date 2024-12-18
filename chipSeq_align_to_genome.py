import subprocess
import os
from collections import OrderedDict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--rawReadsDir", type=str, help="path to raw reads")
parser.add_argument(
    "-g", "--genome", type=str, help="path to genome *fasta* file"
)
parser.add_argument("-o", "--output", type=str, help="path to output sam files")
parser.add_argument("-p", "--processers", type=int, help="number of cpu to use")
parser.add_argument(
    "-pe", "--pairend", action="store_true", help="set action to pair end"
)
parser.add_argument(
    "-re",
    "--rawReadExt",
    type=str,
    default=".fq.gz",
    help="extension of raw reads file",
)

args = parser.parse_args()

rawReadsDir = args.rawReadsDir
pe = args.pairend
rawReadsExt = args.rawReadExt
genome = args.genome
samOutDir = args.output
ncpu = str(args.processers)

fileList = {}

# index genome
genomeBowtie2idx = f"{genome}.bowtie2idx"
doneidx = False
try:
    idxdir = os.path.dirname(os.path.realpath(genomeBowtie2idx))
    for f in os.listdir(idxdir):
        if f.endswith("bt2"):
            doneidx = True
            break
except:
    pass
if not doneidx:
    p = subprocess.run(
        ["bowtie2-build", genome, genomeBowtie2idx], capture_output=True
    )
    if not p.returncode == 0:
        print(f"Command: {' '.join(p.args)}")
        print(f"Error message: {p.stderr.decode()}")
        raise RuntimeError


if not os.path.isdir(samOutDir):
    os.makedirs(samOutDir)

# generate file list to align
# one data file per sample
pairSep = "_"
TO_REMOVE_SUFEIXES = ["_1", "_2", "_R1", "_R2", "_R1_001", "_R2_001"]
samples = OrderedDict()
for path, subdirs, files in os.walk(rawReadsDir):
    for name in files:
        if name.endswith(rawReadsExt):
            fqFile = os.path.join(path, name)
            sname = name[: -len(rawReadsExt)]
            if pe:
                sname, pairn = sname.split(pairSep)
            else:
                pairn = "0"
                raw_data_dir_name = os.path.split(path)[-1]
                for suf in TO_REMOVE_SUFEIXES:
                    if sname.endswith(suf):
                        sname = sname[: -len(suf)]
                        break
                if sname not in raw_data_dir_name:
                    sname = os.path.split(path)[-1] + "_" + sname
            if sname not in samples:
                samples[sname] = {}
            samples[sname][pairn] = fqFile


for s in samples:
    samOut = os.path.join(samOutDir, f"{s}.sam")
    args = [
        "bowtie2",
        "-p",
        ncpu,
        "-x",
        genomeBowtie2idx,
    ]
    if pe:
        p1, p2 = sorted(samples[s].keys())
        args.extend(
            [
                "-1",
                samples[s][p1],
                "-2",
                samples[s][p2],
            ]
        )
    else:
        args.extend(
            [
                "-U",
                samples[s]["0"],
            ]
        )
    args.extend(["-S", samOut])
    print("Running...\n", " ".join(args))
    p = subprocess.run(args, capture_output=True)
    bowtie2_output = p.stdout.decode("utf-8")
    if len(bowtie2_output.strip()) > 0:
        print("Output of bowtie2:\n", bowtie2_output, "\n")
    print("Output file:\n", samOut, "\n")
