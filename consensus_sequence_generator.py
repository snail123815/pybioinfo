import argparse
from pathlib import Path

from Bio import SeqIO

from pyBioinfo_modules.bio_sequences.vcf_parser import (
    applyVariancesOnSeqRecords, vcfParser)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path)
    parser.add_argument("--genome", type=Path, help="genbank format")
    parser.add_argument("--output", type=Path)

    args = parser.parse_args()

    varianceDatas = vcfParser(args.vcf)

    variantSeqRecordDict = applyVariancesOnSeqRecords(
        varianceDatas, SeqIO.parse(args.genome, "genbank")
    )

    SeqIO.write(variantSeqRecordDict.values(), args.output, "genbank")


if __name__ == "__main__":
    main()
