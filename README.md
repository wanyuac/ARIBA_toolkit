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

