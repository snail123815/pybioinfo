# Tree for all homologues protein

Make a tree for homologues proteins found by `phmmer` from a defined database. (Use sensitive if comprehensiveness is needed)

## Challanges

1. `phmmer` output domain hits only, while homologous protein needs to be much more similar. A hit should represent and only represent the query protein.
2. Output table is not clearly parseable.
3. Using "sensitive" mode is essential, may result too many hits.

## How

From search result, download
1. JSON result file
2. alignment fasta file (.afa).

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
