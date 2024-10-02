"""Build correlation between two genomes, find homologue for each protein.
BLAST every protein in query file against target database,
give the best hit that score (e and coverage) higher than required threshold.

By DU Chao
c.du@biology.leidenuniv.nl
durand.dc@hotmail.com
"""

import argparse
import os
import subprocess
from tempfile import NamedTemporaryFile
from Bio.Blast import NCBIXML


Parser = argparse.ArgumentParser(
    description=(
        "Build correlation between two genomes, "
        "find homologue for each protein."
        "BLAST every protein in query file against target database,"
        "give the best hit that score (e and coverage) higher than threshold."
    ),
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
Parser.add_argument(
    "targetDatabase", type=str, help="Target protein BLAST database path"
)
Parser.add_argument(
    "queryProteinFasta", type=str, help="Query protein file in fasta format"
)
Parser.add_argument(
    "-qname", type=str, help="name of input proteome", default="query"
)
Parser.add_argument(
    "-tDbName",
    type=str,
    help="name of target proteome",
    default="target_proteome",
)
Parser.add_argument(
    "-out",
    type=str,
    help=(
        "tsv format output file, default to "
        "[queryFile]_[tDbName]_corr_e[evalue]_cover[coverage].tsv"
    ),
    default="",
)
Parser.add_argument(
    "-eThresh", type=float, help="evalue threshold", default=0.001
)
Parser.add_argument("-cThresh", type=int, help="coverage cutoff", default=60)
Parser.add_argument("-cpu", type=int, help="Threads to use.", default=8)

args = Parser.parse_args()

dbFile = args.targetDatabase.strip()
protsFasta = args.queryProteinFasta.strip()
qname = args.qname.strip()
tDbName = args.tDbName.strip()
tableOut = args.out.strip()
eThresh = args.eThresh
cThresh = args.cThresh

for n in [dbFile, protsFasta, qname, tDbName, tableOut]:
    fn = os.path.splitext(n)[0]
    for s in [".", " "]:
        assert (
            s not in fn
        ), f'"{s}" not allowed in the file names. Check: {n}\n{fn}.'


if tableOut == "":
    tableOut = os.path.splitext(protsFasta)[0] + f"_{tDbName}"
tableOut += "_corr"
tableOut += f"_e{eThresh}_cover{cThresh}.tsv"

dbPath, dbName = os.path.split(dbFile)
dbName = os.path.splitext(dbName)[0]

if dbPath != "":  # current dir
    os.chdir(dbPath)

blastResTemp = NamedTemporaryFile()

blastpCmd = [
    "blastp",
    "-query",
    protsFasta,
    "-db",
    dbName,
    "-outfmt",
    "5",
    "-out",
    blastResTemp.name,
    "-num_threads",
    str(args.cpu),
    "-evalue",
    str(eThresh),
    "-max_hsps",
    "1",
    "-max_target_seqs",
    "1",
]


print("Running BLAST")
print(" ".join(blastpCmd))

result = subprocess.run(blastpCmd, capture_output=True, text=True)

if result.returncode != 0:
    print("Error in BLAST")
    print(result.stderr)
    exit(1)

print(f"Writing result in table {tableOut}")
with open(tableOut, "w", encoding="utf-8") as outHandle:
    outHandle.write(f"{qname}\t{tDbName}\tCoverage\n")
    with open(blastResTemp.name, "r", encoding="utf-8") as resHandle:
        records = NCBIXML.parse(resHandle)
        for record in records:
            queryName = record.query
            target = "NotFound"
            coverage = 0
            if len(record.alignments) != 0:
                firstAlignment = record.alignments[0]
                hsp = firstAlignment.hsps[0]
                hspCoverage = hsp.align_length / record.query_length
                if hspCoverage >= cThresh / 100:
                    target = firstAlignment.hit_def
                    coverage = hspCoverage
            outHandle.write(f"{queryName}\t{target}\t{coverage:.0%}\n")

# TODO
# If database not provided, make one in the temperary directory.
# TODO
# phmmer output parsing
# Comparing blastp
# Produces similar results in terms of homolog detection. Searching with the
# sequence from PDB ID 2abl, chain A against PDB yields 244 matches compared
# with 214 matches for phmmer and blastp, respectively, using an E-value
# threshold of 0.01 and default search parameters. The matches were inspected
# for the presence of an SH3 (Src homology 3) and/or SH2 (Src homology 2)
# domain(s). phmmer results have the added advantage of scoring each residue
# in the alignment, giving users an indication of the parts of the alignment
# that are trustworthy. HMMER web server allows configuration of the cut-offs
# and provides access to all matches.

# --tblout
# #                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# # target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
# #------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
# tr|A0A178UTI0|A0A178UTI0_ARATH -          TMLOC_00011          -            3.9e-18   68.1   0.0   4.9e-18   67.8   0.0   1.1   1   0   0   1   1   1   1 Carbonic anhydrase OS=Arabidopsis thaliana OX=3702 GN=At5g14740 PE=3 SV=1
# tr|A0A178UT81|A0A178UT81_ARATH -          TMLOC_00011          -              4e-18   68.1   0.0   5.1e-18   67.7   0.0   1.1   1   0   0   1   1   1   1 Carbonic anhydrase OS=Arabidopsis thaliana OX=3702 GN=At5g14740 PE=3 SV=1
# tr|A0A654G125|A0A654G125_ARATH -          TMLOC_00011          -            4.4e-18   67.9   0.0   5.7e-18   67.6   0.0   1.1   1   0   0   1   1   1   1 Carbonic anhydrase OS=Arabidopsis thaliana OX=3702 GN=At5g14740 PE=3 SV=1
# tr|A0A178V525|A0A178V525_ARATH -          TMLOC_00011          -            4.6e-18   67.9   0.0   6.9e-18   67.3   0.0   1.2   1   0   0   1   1   1   1 Carbonic anhydrase OS=Arabidopsis thaliana OX=3702 GN=AXX17_At4g38400 PE=3 SV=1
# tr|Q56X90|Q56X90_ARATH         -          TMLOC_00011          -            4.6e-18   67.9   0.0   5.7e-18   67.5   0.0   1.0   1   0   0   1   1   1   1 Carbonic anhydrase OS=Arabidopsis thaliana OX=3702 GN=At3g01500 PE=2 SV=1
# sp|Q94CE3|BCA5_ARATH           -          TMLOC_00011          -            4.6e-18   67.9   0.0   6.9e-18   67.3   0.0   1.2   1   0   0   1   1   1   1 Beta carbonic anhydrase 5, chloroplastic OS=Arabidopsis thaliana OX=3702 GN=BCA5 PE=2 SV=1
# tr|A0A654FVP3|A0A654FVP3_ARATH -          TMLOC_00011          -            4.6e-18   67.9   0.0   6.9e-18   67.3   0.0   1.2   1   0   0   1   1   1   1 Carbonic anhydrase OS=Arabidopsis thaliana OX=3702 GN=AN1_LOCUS20184 PE=3 SV=1
# tr|A0A5S9YAV2|A0A5S9YAV2_ARATH -          TMLOC_00012          -           9.5e-189  631.5   1.2  1.1e-188  631.3   1.2   1.0   1   0   0   1   1   1   1 Acetyl-coenzyme A synthetase OS=Arabidopsis thaliana OX=3702 GN=C24_LOCUS23666 PE=3 SV=1
# tr|A0A654G619|A0A654G619_ARATH -          TMLOC_00012          -           4.8e-188  629.2   1.1  5.6e-188  629.0   1.1   1.0   1   0   0   1   1   1   1 Acetyl-coenzyme A synthetase OS=Arabidopsis thaliana OX=3702 GN=AN1_LOCUS23823 PE=3 SV=1
# tr|A0A7G2FDW6|A0A7G2FDW6_ARATH -          TMLOC_00012          -           1.7e-181  607.6   0.9    2e-181  607.3   0.9   1.0   1   0   0   1   1   1   1 Acetyl-coenzyme A synthetase OS=Arabidopsis thaliana OX=3702 GN=AT9943_LOCUS20608 PE=3 SV=1
# tr|A0A1P8BCV2|A0A1P8BCV2_ARATH -          TMLOC_00012          -           7.4e-170  569.1   2.7  8.7e-170  568.9   2.7   1.0   1   0   0   1   1   1   1 Acetyl-coenzyme A synthetase OS=Arabidopsis thaliana OX=3702 GN=ACS PE=1 SV=1
# tr|A0A654EMB9|A0A654EMB9_ARATH -          TMLOC_00012          -            3.7e-46  160.3   0.0   2.5e-21   78.3   0.0   2.9   1   1   1   2   2   2   2 4-coumarate--CoA ligase-like 4 OS=Arabidopsis thaliana OX=3702 GN=AN1_LOCUS2180 PE=4 SV=1
# tr|A0A7G2DS40|A0A7G2DS40_ARATH -          TMLOC_00012          -              7e-39  136.3   0.0   8.9e-21   76.5   0.0   2.9   1   1   0   2   2   2   2 (thale cress) hypothetical protein OS=Arabidopsis thaliana OX=3702 GN=AT9943_LOCUS1770 PE=4 SV=1
# sp|F4KBF3|AAE17_ARATH          -          TMLOC_00012          -            1.1e-38  135.6   0.0   1.5e-38  135.2   0.0   1.1   1   0   0   1   1   1   1 Probable acyl-activating enzyme 17, peroxisomal OS=Arabidopsis thaliana OX=3702 GN=AAE17 PE=2 SV=1
# tr|A0A654G3L4|A0A654G3L4_ARATH -          TMLOC_00012          -            1.1e-38  135.6   0.0   1.5e-38  135.2   0.0   1.1   1   0   0   1   1   1   1 AMP-dependent synthetase/ligase domain-containing protein OS=Arabidopsis thaliana OX=3702 GN=AN1_LOCUS23043 PE=4 SV=1

# --domtblout
# #                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# # target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
# #------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
# tr|A0A178UTI0|A0A178UTI0_ARATH -            256 TMLOC_00011          -            782   3.9e-18   68.1   0.0   1   1   3.5e-21   4.9e-18   67.8   0.0   576   758    71   239    48   255 0.78 Carbonic anhydrase OS=Arabidopsis thaliana OX=3702 GN=At5g14740 PE=3 SV=1
# tr|A0A178UT81|A0A178UT81_ARATH -            259 TMLOC_00011          -            782     4e-18   68.1   0.0   1   1   3.6e-21   5.1e-18   67.7   0.0   576   758    74   242    51   258 0.78 Carbonic anhydrase OS=Arabidopsis thaliana OX=3702 GN=At5g14740 PE=3 SV=1
# tr|A0A654G125|A0A654G125_ARATH -            267 TMLOC_00011          -            782   4.4e-18   67.9   0.0   1   1     4e-21   5.7e-18   67.6   0.0   576   758    82   250    60   266 0.78 Carbonic anhydrase OS=Arabidopsis thaliana OX=3702 GN=At5g14740 PE=3 SV=1
# tr|A0A178V525|A0A178V525_ARATH -            301 TMLOC_00011          -            782   4.6e-18   67.9   0.0   1   1   4.8e-21   6.9e-18   67.3   0.0   584   758   108   270    98   293 0.78 Carbonic anhydrase OS=Arabidopsis thaliana OX=3702 GN=AXX17_At4g38400 PE=3 SV=1
# tr|Q56X90|Q56X90_ARATH         -            259 TMLOC_00011          -            782   4.6e-18   67.9   0.0   1   1     4e-21   5.7e-18   67.5   0.0   577   754    75   238    55   250 0.81 Carbonic anhydrase OS=Arabidopsis thaliana OX=3702 GN=At3g01500 PE=2 SV=1
# sp|Q94CE3|BCA5_ARATH           -            301 TMLOC_00011          -            782   4.6e-18   67.9   0.0   1   1   4.8e-21   6.9e-18   67.3   0.0   584   758   108   270    98   293 0.78 Beta carbonic anhydrase 5, chloroplastic OS=Arabidopsis thaliana OX=3702 GN=BCA5 PE=2 SV=1
# tr|A0A654FVP3|A0A654FVP3_ARATH -            301 TMLOC_00011          -            782   4.6e-18   67.9   0.0   1   1   4.8e-21   6.9e-18   67.3   0.0   584   758   108   270    98   293 0.78 Carbonic anhydrase OS=Arabidopsis thaliana OX=3702 GN=AN1_LOCUS20184 PE=3 SV=1
# tr|A0A178VMM4|A0A178VMM4_ARATH -            658 TMLOC_00011          -            782    0.0004   21.8  12.2   1   1   1.1e-06    0.0015   19.8  12.2    45   435    81   488    74   501 0.75 SULTR3 OS=Arabidopsis thaliana OX=3702 GN=AXX17_At3g46280 PE=3 SV=1
# tr|A0A5S9XKR0|A0A5S9XKR0_ARATH -            658 TMLOC_00011          -            782    0.0004   21.8  12.2   1   1   1.1e-06    0.0015   19.8  12.2    45   435    81   488    74   501 0.75 (thale cress) hypothetical protein OS=Arabidopsis thaliana OX=3702 GN=AT9943_LOCUS13378 PE=3 SV=1
# tr|A0A654FEZ8|A0A654FEZ8_ARATH -            658 TMLOC_00011          -            782    0.0004   21.8  12.2   1   1   1.1e-06    0.0015   19.8  12.2    45   435    81   488    74   501 0.75 STAS domain-containing protein OS=Arabidopsis thaliana OX=3702 GN=AN1_LOCUS15552 PE=3 SV=1
# tr|A0A1I9LNF9|A0A1I9LNF9_ARATH -            703 TMLOC_00011          -            782   0.00047   21.5  12.5   1   1   1.1e-06    0.0016   19.7  12.5    45   435   126   533   119   548 0.75 Sulfate transporter 31 OS=Arabidopsis thaliana OX=3702 GN=SULTR3 PE=1 SV=1
# tr|A0A178UCJ4|A0A178UCJ4_ARATH -            693 TMLOC_00012          -            659  6.9e-189  632.0   1.2   1   1  1.3e-191    8e-189  631.8   1.2    30   649    57   688    25   692 0.91 Acetyl-coenzyme A synthetase OS=Arabidopsis thaliana OX=3702 GN=AXX17_At5g33780 PE=3 SV=1
# tr|A0A178UC31|A0A178UC31_ARATH -            743 TMLOC_00012          -            659  8.8e-189  631.6   1.2   1   1  1.7e-191    1e-188  631.4   1.2    30   649   107   738    75   742 0.91 Acetyl-coenzyme A synthetase OS=Arabidopsis thaliana OX=3702 GN=AXX17_At5g33780 PE=3 SV=1
# sp|B9DGD6|ACS_ARATH            -            743 TMLOC_00012          -            659  9.5e-189  631.5   1.2   1   1  1.9e-191  1.1e-188  631.3   1.2    30   649   107   738    75   742 0.91 Acetyl-coenzyme A synthetase, chloroplastic/glyoxysomal OS=Arabidopsis thaliana OX=3702 GN=ACS PE=1 SV=1
# tr|A0A5S9YAV2|A0A5S9YAV2_ARATH -            743 TMLOC_00012          -            659  9.5e-189  631.5   1.2   1   1  1.9e-191  1.1e-188  631.3   1.2    30   649   107   738    75   742 0.91 Acetyl-coenzyme A synthetase OS=Arabidopsis thaliana OX=3702 GN=C24_LOCUS23666 PE=3 SV=1
# tr|A0A654G619|A0A654G619_ARATH -            743 TMLOC_00012          -            659  4.8e-188  629.2   1.1   1   1  9.3e-191  5.6e-188  629.0   1.1    18   649    96   738    75   742 0.92 Acetyl-coenzyme A synthetase OS=Arabidopsis thaliana OX=3702 GN=AN1_LOCUS23823 PE=3 SV=1
# tr|A0A7G2FDW6|A0A7G2FDW6_ARATH -            683 TMLOC_00012          -            659  1.7e-181  607.6   0.9   1   1  3.3e-184    2e-181  607.3   0.9    18   649    46   678    25   682 0.91 Acetyl-coenzyme A synthetase OS=Arabidopsis thaliana OX=3702 GN=AT9943_LOCUS20608 PE=3 SV=1
# tr|A0A1P8BCV2|A0A1P8BCV2_ARATH -            610 TMLOC_00012          -            659  7.4e-170  569.1   2.7   1   1  1.4e-172  8.7e-170  568.9   2.7    29   559    56   597    25   607 0.91 Acetyl-coenzyme A synthetase OS=Arabidopsis thaliana OX=3702 GN=ACS PE=1 SV=1
# tr|A0A654EMB9|A0A654EMB9_ARATH -            824 TMLOC_00012          -            659   3.7e-46  160.3   0.0   1   2   1.1e-23   6.6e-21   76.9   0.0   112   586    62   508    55   517 0.79 4-coumarate--CoA ligase-like 4 OS=Arabidopsis thaliana OX=3702 GN=AN1_LOCUS2180 PE=4 SV=1
# tr|A0A654EMB9|A0A654EMB9_ARATH -            824 TMLOC_00012          -            659   3.7e-46  160.3   0.0   2   2   4.1e-24   2.5e-21   78.3   0.0   351   627   552   819   515   823 0.88 4-coumarate--CoA ligase-like 4 OS=Arabidopsis thaliana OX=3702 GN=AN1_LOCUS2180 PE=4 SV=1
# tr|A0A7G2DS40|A0A7G2DS40_ARATH -            761 TMLOC_00012          -            659     7e-39  136.3   0.0   1   2   1.5e-23   8.9e-21   76.5   0.0   112   582    62   504    55   513 0.78 (thale cress) hypothetical protein OS=Arabidopsis thaliana OX=3702 GN=AT9943_LOCUS1770 PE=4 SV=1
# tr|A0A7G2DS40|A0A7G2DS40_ARATH -            761 TMLOC_00012          -            659     7e-39  136.3   0.0   2   2   6.2e-17   3.7e-14   54.6   0.0   351   562   552   758   514   761 0.86 (thale cress) hypothetical protein OS=Arabidopsis thaliana OX=3702 GN=AT9943_LOCUS1770 PE=4 SV=1
# sp|F4KBF3|AAE17_ARATH          -            721 TMLOC_00012          -            659   1.1e-38  135.6   0.0   1   1   2.5e-41   1.5e-38  135.2   0.0    69   627   156   711   128   719 0.81 Probable acyl-activating enzyme 17, peroxisomal OS=Arabidopsis thaliana OX=3702 GN=AAE17 PE=2 SV=1
# tr|A0A654G3L4|A0A654G3L4_ARATH -            721 TMLOC_00012          -            659   1.1e-38  135.6   0.0   1   1   2.5e-41   1.5e-38  135.2   0.0    69   627   156   711   128   719 0.81 AMP-dependent synthetase/ligase domain-containing protein OS=Arabidopsis thaliana OX=3702 GN=AN1_LOCUS23043 PE=4 SV=1

# hmm coord - hmm profile coordinates - query sequence
# env coord - envelope coordinates - target sequence

# Doubts:
# Domain table, always 1 of 1?
#     What happens if 1 of 2 or more?
# Target table, some have two domains, what and why?

# Parsing tblout could be achieved by simple IO. Making use of know column
# numbers (or get it from the third line), split by " ", join the exceeding
# elements by " ".

# TODO
# PSI-BLAST vs. jackhmmer
# These are iterative search methods.
# When comparing PSI-BLAST and JACKHMMER:
# - Algorithm: PSI-BLAST uses a position-specific scoring matrix (PSSM) approach,
#   while JACKHMMER uses profile hidden Markov models (HMMs).
# - Iterative Process: Both tools use an iterative approach to refine the search
#   results and improve sensitivity.
# - Applications: While both tools are used for sequence similarity searches,
#   PSI-BLAST is often used for comparing individual protein sequences, while
#   JACKHMMER is particularly useful for analyzing entire protein families.
# - Sensitivity: JACKHMMER is often considered more sensitive than PSI-BLAST in
#   detecting remote homologs.
