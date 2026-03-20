process DADA2 {
    publishDir "${params.outdir}/03_dada2", mode: 'copy'
    
    input:
    val sample_ids
    path read1_files
    path read2_files
    
    output:
    path "asv_table.tsv", emit: asv_table
    path "track_control.tsv", emit: track_control
    path "ASV_sequences.fasta", emit: fasta
    
    script:
    // IMPORTANT: pairing relies on positional alignment
    // ids[i], r1s[i], r2s[i] must correspond to the same sample
    def ids = sample_ids.join(' ')
    def r1s = read1_files.join(' ')
    def r2s = read2_files.join(' ')
    
    """
    Rscript ${projectDir}/bin/dada2.R \\
        --sample_ids ${ids} \\
        --input_R1 ${r1s} \\
        --input_R2 ${r2s} \\
        --nproc ${task.cpus} \\
        --truncQ ${params.truncQ} \\
        --truncLen_R1 ${params.truncLen_R1} \\
        --truncLen_R2 ${params.truncLen_R2} \\
        --maxEE_R1 ${params.maxEE_R1} \\
        --maxEE_R2 ${params.maxEE_R2} \\
        --maxMismatch ${params.maxMismatch} \\
        --minOverlap ${params.minOverlap}

    """
}