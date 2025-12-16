# Tree for all homologues protein

Make a tree for homologues proteins found by `phmmer` from a defined database. (Use sensitive if comprehensiveness is needed)

## Challanges

1. `phmmer` output domain hits only, while homologous protein needs to be much more similar. A hit should represent and only represent the query protein.
2. Output table is not clearly parseable.
3. Using "sensitive" mode is essential, may result too many hits, or too slow. (`diamond blastp --ultra-sensitive`, `jackhmmer --max`)

## How

### HMMER web search

(Not available now)

From search result, download
1. JSON result file
2. alignment fasta file (.afa).

Parse JSON file, find the hit domains, if this domain covers the presumed essential region, then keep the protein, else discard it.

After we have the selected protein list, keep the longest alignment from alignment file. Use the result file to make a tree.

For each tree node, add species information for easier interpration.

### HMMER search

Use `-A <f> : save multiple alignment of hits to file <f>` to get alignment file. This is a `# STOCKHOLM 1.0` file, should name `*.stockholm`. This file cannot be used by `fasttree`, need to be converted to fasta file using:

```sh
esl-reformat afa a.stockholm > a.afa
esl-reformat -o a.afa afa a.stockholm
# Usage: esl-reformat [-options] <format> <seqfile>
# afa for aligned fasta format
```

Use `--max` for max sensitivity. The cost is too much, increase one digit to `--F1`, `--F2`, `--F3` default values can be an option if missing hits.

Parse database or output header (e.g. `#` lines in stockholm file) for strain information if exists.

## How to use

### 1. `filter_phmmer_alignment.py`

TODO: add argument parse

Change parameters:
- `jsonFile`, `alignmentFasta` for input files
- `tStart`, `tEnd` for start and end position of your target region on query protein sequence. One 'domain' of a 'hit' must cover the selected region.
- [not implemented]`eThresh` for threshold of the full 'hit', a valid 'hit' must have a 'evalue' lower than this value.

Use `esl-alimask` - remove columns from a multiple sequence alignment, then to fasta, then make tree, depends on the needs.

```
esl-alimask -t [options] msafile coords
(remove a contiguous set of columns at the start and end of an alignment)
for example, 23..100, 23/100, or 23-100 all work
```

### 2. FastTree

TO BE INCORPORATED

### 3. `get_strain_name.py`

TODO: add argument parse

Change parameters:
- `jsonFile` for input phmmer result file
- `treeFile` for input tree
- `relationsFile` for output table of 'proteinID' => 'species name'

Result tree will be added a '.species' before the extension.
