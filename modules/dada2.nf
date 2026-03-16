process DADA2 {
    publishDir "${params.outdir}/03_dada2", mode: 'copy'
    
    input:
    path('*')
    
    output:
    path "asv_table.tsv", emit: asv_table
    path "track_control.tsv", emit: track_control
    path "ASV_sequences.fasta", emit: fasta
    
    script:
    """
    R1_files=\$(ls *R1.trimmed.fastq.gz | sort | tr '\\n' ' ')
    R2_files=\$(ls *R2.trimmed.fastq.gz | sort | tr '\\n' ' ')

    Rscript ${projectDir}/bin/dada2.R \\
        --input_R1 \${R1_files} \\
        --input_R2 \${R2_files} \\
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