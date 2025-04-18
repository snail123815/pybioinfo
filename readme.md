# pyBioinfo

A Python toolkit for bioinformatics analysis workflows.

[![Python Tests](https://github.com/snail123815/pybioinfo/actions/workflows/Python-tests.yml/badge.svg)](https://github.com/snail123815/pybioinfo/actions/workflows/Python-tests.yml)
[![codecov](https://codecov.io/gh/snail123815/pybioinfo/branch/main/graph/badge.svg)](https://codecov.io/gh/snail123815/pybioinfo)

## Features

- RNA-Seq analysis pipeline
- ChIP-Seq processing tools
- Genome annotation utilities
- Sequence manipulation tools
- Wrappers for popular bioinformatics tools:
  - Bowtie2
  - featureCounts
  - antiSMASH
  - BigScape
  - MASH

## Installation

```bash
# Clone the repository
git clone https://github.com/snail123815/pybioinfo.git
cd pybioinfo

# Install dependencies
pip install -r requirements.txt
```

## Usage Examples

### RNA-Seq Analysis

```bash
python rnaSeq_raw_to_counts.py --raw path/to/reads/*.fastq.gz --out results/ --gbk genome.gbk --isPe --ncpu 4
```

See [RNA-Seq documentation](howto_RNA-Seq.md) for more details.

### antiSMASH Analysis

```bash
python run_antismash.py path/to/genomes/*.gbk --completeness 2 --threads 4
```

## Module Coverage

| Module | Coverage |
|--------|----------|
| bio_sequences | ![Coverage](https://img.shields.io/badge/coverage-69%25-brightgreen) |
| wrappers | ![Coverage](https://img.shields.io/badge/coverage-39%25-brightgreen) |
| basic | ![Coverage](https://img.shields.io/badge/coverage-16%25-brightgreen) |

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.