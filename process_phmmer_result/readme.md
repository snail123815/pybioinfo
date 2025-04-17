# Process phmmer result - to make a phylogentic tree

## Why

`phmmer` provide multiple options but output table is not well under control. For example:

I want to find, based on the hmm sequence similarity algorithm, all homologous proteins of ProtA on UniProt Reference Proteomes database. Do a search on https://www.ebi.ac.uk/Tools/hmmer/search/phmmer resulted in many hits, a lot of may have multiple motif or part of the query sequence matched, each matched part will have a good score (very low E-value).

However, the goal is to find good homologous protein, which represents and only represents the query protein. Threre may be a possiblity to elimiate partial hits by playing with difference "sequence" and "Hit" E-values, or bit scores. But I don't think it is valid and easy way.

## How

From search result, download the JSON result file and the alignment fasta file (.afa).

Parse JSON file, find the hit domains, if this domain covers the presumed essential region, then keep the protein, else discard it.

After we have the selected protein list, keep the longest alignment from alignment file. Use the result file to make a tree.

For each tree node, add species information for easier interpration.

## Depedencies

```yml
dependencies:
  - numpy
  - termplotlib
  - biopython
  - fasttree
```

## How to use

### 1. `filter_phmmer_alignment.py`

TODO: add argument parse

Change parameters:
- `jsonFile`, `alignmentFasta` for input files
- `tStart`, `tEnd` for start and end position of your target region on query protein sequence. One 'domain' of a 'hit' must cover full region.
- [not implemented]`eThresh` for threshold of the full 'hit', a valid 'hit' must have a 'evalue' lower than this value.

### 2. FastTree

TO BE INCORPORATED

### 3. `get_strain_name.py`

TODO: add argument parse

Change parameters:
- `jsonFile` for input phmmer result file
- `treeFile` for input tree
- `relationsFile` for output table of 'proteinID' => 'species name'

Result tree will be added a '.species' before the extension.
