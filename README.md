# Tools for Processing ARIBA's Outputs

*Yu Wan*



This repository consists of helper scripts for processing outputs of [ARIBA](https://github.com/sanger-pathogens/ariba).



## compileMLST.py

Compile outputs (`mlst_report.tsv`, `mlst_report.details.tsv`, and `report.tsv`) from MLST analysis by ARIBA, assuming the names of ARIBA's output directories are genome names. Namely, the MLST reports should be generated with the command 

```
ariba run db genome1_1.fastq.gz genome1_2.fastq.gz ./genome1
```

where an output directory shares the name of a genome.

**Example command line**

```bash
find ./* -maxdepth 1 -type d | python ~/ARIBA_toolkit/compileMLST.py
```

**Outputs**

- `STs.tsv`: compiled from `mlst_report.tsv` files;
- `scores.tsv`: compiled from `report.tsv` files, which include values of nucleotide identities;
- `hits.tsv`: compiled from `mlst_report.details.tsv`, which contain values of reference coverage (%) and read depths of allele hits.



## concMLSTalleles.py

Concatenate allele sequences of a given list of STs. It is useful for comparing STs at the nucleotide
level and MLST-based phylogenetic construction. Currently this script does not support allele variants,
namely, requires exact hits of query genomes to reference allele sequences.

Example command lines:
```bash
python concMLSTalleles.py -p st_profiles.tsv -s ./sequence_dir -e tfa > mlstSeqs.fna

# Append sequences to an extant FASTA file
python concMLSTalleles.py -p st_profiles_new.tsv -s ./sequence_dir -e tfa >> mlstSeqs.fna
```

This script assumes input ST profiles are offered in a tab-delimited file with columns of STs and
names of MLST loci. (No other columns should be allowed) For example:

| ST   | _locus\_1_ | *locus\_2* | *locus\_3* | *locus\_4* | *locus\_5* | *locus\_6* | *locus\_7* |
| ---- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- |
| 1    | 11         | 21         | 31         | 41         | 51         | 61         | 71         |
| 2    | 12         | 21         | 31         | 41         | 51         | 61         | 71         |
| ...  | ...        | ...        | ...        | ...        | ...        | ...        | ...        |

The directory of MLST allele sequences is expected to be the same as that downloaded from the PubMLST
database using ARIBA. (`pubmlst_download`). The '-s' parameter does not require an end forward slash
for the path. By default, the delimiter for locus name (for instance, gene _adk_) and allele number in
an allele name (for example, _adk.1_) is a period, although some databases use an underscore as the
delimiter.

Dependencies: Python v3, packages [pandas](https://pandas.pydata.org/) and [biopython](https://biopython.org/).

