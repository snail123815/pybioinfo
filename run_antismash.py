import argparse
import logging
from multiprocessing import Pool
from pathlib import Path

from tqdm import tqdm

from pyBioinfo_modules.basic.decompress import (getRootAndFiles,
                                                getStemIfCompressed)
from pyBioinfo_modules.bio_sequences.bio_seq_file_extensions import (
    FNA_EXTENSIONS, GBK_EXTENSIONS)
from pyBioinfo_modules.wrappers.antismash import runAntismash

parser = argparse.ArgumentParser()
parser.add_argument(
    "inputFiles",
    nargs="+",
    help="Input files/folders containing gbk files (can be gz or xz)",
)
parser.add_argument(
    "--parrllel",
    type=int,
    default=8,
    help="Number of antismashes run in parallel",
)
parser.add_argument(
    "--threads",
    type=int,
    default=4,
    help="Number of threads for each antismash run",
)
parser.add_argument(
    "--taxon", type=str, choices=["bacteria", "fungi"], default="bacteria"
)
parser.add_argument(
    "--completeness",
    type=int,
    choices=[1, 2, 3, 10],
    help=f"Completeness of antismash run",
    default=2,
)
parser.add_argument("--dry", action="store_true")
parser.add_argument("--overwrite", action="store_true")
parser.add_argument(
    "--geneFinding",
    type=str,
    choices=["glimmerhmm", "prodigal", "prodigal-m", "auto", "error"],
    help='--geneFinding-tool for antismash, except "auto".'
    + 'It means "none" when the input is gbk file'
    + " which should contain annotation.\n"
    + "The antismash logic based on the fact that gbk file canbe both"
    + "with or without annotation. But sometimes it miss interprate that."
    + "Implementation of this program is that:"
    + 'If input is genbank file, use "none" if set to "auto"'
    + 'If input is DNA sequence file, use "prodigal" if set to "auto"',
    default="auto",
)
args = parser.parse_args()

ALLOW_GENE_FINDING = args.geneFinding != "error"
ALLOWED_EXTENSIONS = (
    GBK_EXTENSIONS + FNA_EXTENSIONS if ALLOW_GENE_FINDING else GBK_EXTENSIONS
)

outputRoot: Path
targetFiles: list[Path]
outputRoot, targetFiles = getRootAndFiles(args.inputFiles, ALLOWED_EXTENSIONS)
pathOut = (
    outputRoot.parent / f"{outputRoot.name}_antismash_level{args.completeness}"
)
pathOut.mkdir(exist_ok=True)
logging.basicConfig(
    filename=pathOut / "antismash.log", filemode="a", level="INFO"
)
logging.info("#" * 100)

runnerPool = Pool(args.parrllel)
results = []
for f in targetFiles:
    antismashDir = pathOut / getStemIfCompressed(f)
    results.append(
        runnerPool.apply_async(
            runAntismash,
            kwds={
                "inputFilePath": f,
                "title": f.name,
                "taxon": args.taxon,
                "completeness": args.completeness,
                "cpu": args.threads,
                "output": antismashDir,
                "silent": True,
                "dry": args.dry,
                "geneFinding": args.geneFinding,
                "overwrite": args.overwrite,
                "existsOk": True,
            },
        )
    )
runnerPool.close()

returns = []
for res in tqdm(results):
    returns.append(res.get())
logging.info("Result files in:")
returns.sort(key=lambda x: x.name)
for ret in returns:
    logging.info(ret)
logging.info("#" * 100)
