process MERGING_READS {
    publishDir "${params.outdir}/04_bbmap", mode: 'copy'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("${sample_id}-mergedpairs.fastq.gz"), emit: reads
    
    script:
    """
    bbmerge.sh \\
    in1=${read1} \\
    in2=${read2} \\
    qtrim=r \\
    trimq=${params.bbmap_trimq} \
    maxlength=${params.bbmap_maxlength} \
    mininsert=${params.bbmap_mininsert} \
    threads=${task.cpus} \
    out=${sample_id}-mergedpairs.fastq.gz \
    outu=/dev/null
    """
}