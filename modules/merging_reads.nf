process MERGING_READS {
    publishDir "${params.outdir}/04_merged", mode: 'copy'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("${sample_id}-mergedpairs.fastq.gz"), emit: reads
    
    script:
    """
    vsearch -fastq_mergepairs ${read1} -reverse ${read2} \\
    -fastq_maxdiffs ${params.maxdiff} \\
    -fastq_maxdiffpct ${params.maxdiffpct} \\
    -fastqout ${sample_id}-mergedpairs.fastq \\
    -fastq_truncqual ${params.trunctail} \\
    -fastq_minmergelen ${params.minmergelen} \\
    -fastq_maxmergelen ${params.maxmergelen}

    gzip -c ${sample_id}-mergedpairs.fastq > ${sample_id}-mergedpairs.fastq.gz  
    rm ${sample_id}-mergedpairs.fastq
    """
}