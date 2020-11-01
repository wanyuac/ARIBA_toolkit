#!/usr/bin/env nextflow

/*
Submit ARIBA jobs for gene detection.

[Use guide]
Dependency: Singularity

To run this pipeline in a screen session:
    nextflow -Djava.io.tmpdir=$PWD run ariba.nf --fastq "$PWD/*_{1,2}.fastq.gz" --our_dir sample1 \
    -c ariba.config -profile pbs -with-singularity $HOME/software/docker/ariba.sif

[Declaration]
Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public License v3.0
Publication: 1 Nov 2020

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
ariba_out_dir = mkdir(params.out_dir)

read_sets = Channel.fromFilePairs(params.fastq)

process Ariba {
    input:
    set genome, file(paired_fastq) from read_sets  // Genome ID and paired read files
    
    script:    
    """
    ariba run --assembler spades --spades_mode wgs --assembly_cov 50 --spades_options '--only-assembler --careful' --threads 2 --assembled_threshold 0.95 --min_scaff_depth 10 --force --nucmer_min_id 90 --tmp_dir $PWD ${params.db} ${paired_fastq[0]} ${paired_fastq[1]} $genome
    mv ./$genome ${ariba_out_dir}
    """
}
