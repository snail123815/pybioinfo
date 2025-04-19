<!-- This file is licensed under the MIT License -->

# RNA-Seq data processing

From raw reads to different gene expression.

## Prerequisites

### Using `Align_RNASeq_bowtie2.py`
    
`bowtie2` - v2.3.5.1 used

### Using `CountReads_RNASeq_featureCounts.py`

`featureCounts` - v2.0.1 used

### Jupyter notebooks

Require python3.9 and R kernel.

## Usage

Copy the `*.py` scripts and notebooks to the root folder of your data. Scripts needs to be edited before run. Notebooks are suitable for adjustments in data analysis.

## RUN script on ALICE from Leiden University

Check `testJob_on_ALICE.slurm` for test run setup.

### Prerequisites on ALICE

Time stamp: 20210305 (Not actual commands, so no output listed, please mind for errors)

```shell
[duc@nodelogin01 ~]$ module load Miniconda3
[duc@nodelogin01 ~]$ conda init bash
(base) [duc@nodelogin01 ~]$ source .bashrc
(base) [duc@nodelogin01 ~]$ conda config --add channels defaults # should already be there
(base) [duc@nodelogin01 ~]$ conda config --add channels bioconda
(base) [duc@nodelogin01 ~]$ conda config --add channels conda-forge
(base) [duc@nodelogin01 ~]$ conda config --show # not required
(base) [duc@nodelogin01 ~]$ conda config --set auto_activate_base False # not required
(base) [duc@nodelogin01 ~]$ conda deactivate # not required
[duc@nodelogin01 ~]$ conda create rnaseq
[duc@nodelogin01 ~]$ conda activate rnaseq
(rnaseq) [duc@nodelogin01 ~]$ conda install bowtie2=2.4.2 # Latest now, Python included
(rnaseq) [duc@nodelogin01 ~]$ conda install subread
(rnaseq) [duc@nodelogin01 ~]$ conda install tbb=2020.3 # Will override current version.
                                                       # Current version 2021.1.1 cannot be
                                                       # reconized by bowtie2=2.4.2
(rnaseq) [duc@nodelogin01 ~]$ conda install samtools
(rnaseq) [duc@nodelogin01 ~]$ python # test
(rnaseq) [duc@nodelogin01 ~]$ bowtie2 # test
(rnaseq) [duc@nodelogin01 ~]$ samtools # test
(rnaseq) [duc@nodelogin01 ~]$ featureCounts # test
```

### Time consumption

Pair end raw data ~640M × 2, align to M145 genome ~8.7M

Assigning 4 G memory,
2 core time usage >22 min
4 core time usage >10 min

Please consider time spend in copying files to and from node storage. The output below is from test job of 4G memory, 2 tasks (cpus). Final notification indicates the total run time is 00:22:32. 

The same job running on a desktop i7 6700 @3.4 GHz 64 GB ram using all resources will cost >13 min.

## TODOs

- Combine two python scripts
- Change behaviour to use arguments as meaningful logging is done by the script.
