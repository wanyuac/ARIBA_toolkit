# Scripts for ARIBA-based gene detection

*Yu Wan*



This repository consists of helper scripts for submitting [ARIBA](https://github.com/sanger-pathogens/ariba) jobs or processing their outputs.

<br />

## 1. A Nextflow pipeline for ARIBA

This pipeline submits PBS (Portable Batch System) jobs of ARIBA to a high-performance computing (HPC) cluster with a control for queue sizes. Users may also opt to run ARIBA locally (without submitting jobs to the HPC system) using this pipeline. Two Nextflow scripts implement this pipeline:

- `ariba.nf`: Main Nextflow script for the pipeline;
- `ariba.config`: Configuration file of `ariba.nf`. By default, each PBS job requests four CPUs and 16 GB RAM from the system. The same request applies to local jobs as well.



**Dependencies**

- Singularity-compatible Docker image of ARIBA
- Nextflow v20.0 and above



**Parameters**

- `fastq`: A glob for read files. It must be braced by a pair of quotation marks and no square brackets should be used. Default: `"*_{1,2}.fastq.gz"`.
- `output_dir`: Name and path for the parental directory of output files. Default: `output`.
- `db_dir`: Name and path to access an existing ARIBA-compatible ResFinder database. Default: `resfinder`.
- `queue_size`: Maximum number of concurrent PBS jobs submitted to the HPC system. Some computation facilities have restrictions for the queue size for fair use. Default: 10.



**Example commands**

```bash
# Install ARIBA's Singularity image through Docker
singularity pull docker://staphb/ariba  # Rename the Image file to ariba.sif

# Run the ARIBA Nextflow pipeline for gene detection
nextflow -Djava.io.tmpdir=$PWD run ariba.nf --fastq "$PWD/reads/*_{1,2}.fastq.gz" --db_dir $HOME/db/resfinder --output_dir output -c ariba.config -profile pbs --queue_size 15 -with-singularity $HOME/software/docker/ariba.sif
```



**Outputs**

Five subdirectories are created under the parental output directory by the pipeline:

- `/report`: `[sample]_report.tsv`
- `/gene`: `[sample]_genes.fna`, decompressed and renamed from ARIBA's output file `assembled_genes.fa.gz`.
- `/stats`: `[sample]_stats.tsv`, renamed from`debug.report.tsv`.
- `/contig`: `[sample]_assemblies.fna.gz` and `[sample]_seqs.fna.gz`, renamed from output files `assemblies.fa.gz` and `assembled_seqs.fa.gz` of each sample, respectively.
- `/log`: `[sample]_log_clusters.gz` and `[sample]_log_version.txt`, renamed from output files `log.clusters.gz` and `version_info.txt` of each sample, respectively.



Users may refer to [instructions for PAMmaker](https://github.com/wanyuac/PAMmaker#ariba) for a method converting this pipeline's outputs to an allelic presence-absence matrix.

<br />

## 2. Processing results of MLST analysis

### compileMLST.py

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



### concMLSTalleles.py

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

