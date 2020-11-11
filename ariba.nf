#!/usr/bin/env nextflow

/*
Submit ARIBA jobs for gene detection. Detection of non-coding sequences is not applicable using this pipeline.

[Use guide]
Dependency: Singularity

To run this pipeline in a screen session:
    nextflow -Djava.io.tmpdir=$PWD run ariba.nf --fastq "$PWD/*_{1,2}.fastq.gz" --our_dir sample1 \
    -c ariba.config -profile pbs -with-singularity $HOME/software/docker/ariba.sif

[Declaration]
Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public License v3.0
Publication: 1 Nov 2020; latest update: 11 Nov 2020.

Reference:
    Martin Hunt (@martibartfast) and Anthony Underwood (@bioinformant) github.com/aunderwo/nextflow_ariba
*/

/********** Function definition **********/
def mkdir(dir_path) {
    /*
    Creates a directory and returns a File object.
    */
    def dir_obj = new File(dir_path)
    
    if ( !dir_obj.exists() ) {
        result = dir_obj.mkdir()
        println result ? "Successfully created directory ${dir_obj}" : "Cannot create directory ${dir_obj}"
    } else {
        println "Directory ${dir_obj} exists."
    }
    
    return dir_obj
}

/********** Pipeline that connects processes through channels **********/
/* 
Creating output directories
An explanation of outputs: github.com/sanger-pathogens/ariba/wiki/Task:-run
In summary:
  - assembled_genes.fa.gz: assembled protein-coding genes, which are actually alleles of reference sequences.
                           Sequence headers use gene names (namely, cluster names) rather than allele names from
                           the reference database. Such gene names are used in ARIBA's compiled report as well.
                           This notation avoids misleading allele names when reference sequences are highly similar.
                           Multiple alleles of the same cluster may present.
  - assemblies.fa.gz: complete, unedited contigs matching to reference sequences. Flanking regions of these
                      sequences are also included.
  - assembled_seqs.fa.gz: contigs that fully or partially match to reference sequences. Allele names and gene names
                          from the reference database are both used in sequence headers.
*/
main_output_dir = params.output_dir  // Parental directory of all outputs
report_dir =  main_output_dir + "/report"  // Key output 1: report.tsv
gene_dir = main_output_dir + "/gene"  // Key output 2: assembled_genes.fa.gz
stats_dir = main_output_dir + "/stats"  // debug.report.tsv
log_dir = main_output_dir + "/log"  // log.clusters.gz and version_info.txt
contig_dir = main_output_dir + "/contig"  // assemblies.fa.gz and assembled_seqs.fa.gz

mkdir(params.output_dir)  // Returns true or false
mkdir(report_dir)
mkdir(gene_dir)
mkdir(stats_dir)
mkdir(log_dir)
mkdir(contig_dir)

/* Get inputs */
read_sets = Channel.fromFilePairs(params.fastq)
ariba_db = file(params.db_dir)

/* Processes */
process gene_detection {
    /* Grep and copy outputs of this process (independent of the output declaration blok) */
    publishDir report_dir, mode: "copy", pattern: "*_report.tsv"
    publishDir gene_dir, mode: "copy", pattern: "*_genes.fna"
    publishDir contig_dir, mode: "copy", pattern: "*.fna.gz"
    publishDir stats_dir, mode: "copy", pattern: "*_stats.tsv"
    publishDir log_dir, mode: "copy", pattern: "*_log_*.*"

    input:
    set genome, file(paired_fastq) from read_sets  // Genome ID and its paired read files

    /*
    Output files must be registered in this block so publishDir can take effect, or no file
    will be copied to output directories. These output files do not have to be sent to any
    channel.
    */
    output:
    file "${genome}_report.tsv" into report_channel  // Grep the file of interest using this glob string and feed it into the next process
    file "${genome}_genes.fna"
    file "${genome}_*.fna.gz"
    file "${genome}_stats.tsv"
    file "${genome}_log_*.*"

    /*
    Runs under the working directory of the current job. This directory differs from $PWD, the directory where
    this Nextflow job is launched. (Try adding "echo $PWD" in the block of bash commands)
    */
    script:
    """
    ariba run --assembler spades --spades_mode wgs --spades_options \"--only-assembler --careful\" --assembly_cov 70 --threads 4 --assembled_threshold 0.95 --min_scaff_depth 10 --force --nucmer_min_id 90 --tmp_dir $PWD ${ariba_db} ${paired_fastq[0]} ${paired_fastq[1]} $genome
    mv ${genome}/report.tsv ${genome}_report.tsv
    gunzip -c ${genome}/assembled_genes.fa.gz > ${genome}_genes.fna
    mv ${genome}/assembled_seqs.fa.gz ${genome}_seqs.fna.gz
    mv ${genome}/assemblies.fa.gz ${genome}_assemblies.fna.gz
    mv ${genome}/debug.report.tsv ${genome}_stats.tsv
    mv ${genome}/log.clusters.gz ${genome}_log_clusters.gz
    mv ${genome}/version_info.txt ${genome}_log_version.txt
    """
}

process summarise_reports {
    publishDir main_output_dir, mode: "copy", pattern: "*.csv"

    input:
    file ariba_report from report_channel.collect()

    output:
    file "*.csv"

    script:
    """
    ariba summary --row_filter n --no_tree --min_id 90.0 ariba_summary ${ariba_report}
    """
}