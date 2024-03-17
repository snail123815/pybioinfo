from pathlib import Path
import argparse
from pyBioinfo_modules.wrappers.prokka import runProkka
from multiprocessing import Pool
from pyBioinfo_modules.bio_sequences.bio_seq_file_extensions import FNA_EXTENSIONS
from pyBioinfo_modules.basic.decompress \
    import getStemIfCompressed, getRootAndFiles



parser = argparse.ArgumentParser(description='Run prokka for fasta files')
parser.add_argument(
    'inputFiles',
    nargs="+",
    help='Input files/folders containing fasta files (can be gz or xz)'
)
parser.add_argument(
    '--parrllel',
    type=int, default=1,
    help="Number of prokka run in parallel"
)
parser.add_argument(
    '--threads', type=int, default=4,
    help="Number of threads for each prokka run"
)
parser.add_argument(
    '--dry', action='store_true'
)

args = parser.parse_args()


outputRoot, targetFiles = getRootAndFiles(args.inputFiles, FNA_EXTENSIONS)

runnerPool = Pool(args.parrllel)
results = []
for f in targetFiles:
    prokkaDir = outputRoot / (getStemIfCompressed(f) + '_prokka')
    results.append(
        runnerPool.apply_async(
            runProkka,
            kwds={
                "fastaPath": f,
                "gcode": 11,
                "gram": 'pos',
                "center": 'MBT',
                "genus": None,
                "species": None,
                "strain": None,
                "locustag": None,
                "cpu": args.threads,
                "output": prokkaDir,
                "dry": args.dry,
                "silent": True
            }
        )
    )

returns = []
for res in results:
    returns.append(res.get())