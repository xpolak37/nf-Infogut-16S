process VSEARCH_UNOISE3 {
    publishDir "${params.outdir}/06_unoise", mode: 'copy'

    input:
    path(reads)

    output:
    path("asv_table.tsv"), emit: asv_table
    path("ASV_sequences.fasta"), emit: fasta

    script:
    """
    # 1. Concatenate and decompress all oriented reads, labeling each read
    #     with its sample name so vsearch --otutabout produces per-sample
    #     columns.  Uses ;sample=NAME; annotation (usearch/vsearch standard)
    #     to avoid dot-splitting ambiguity.  Strip -oriented suffix.
    for f in *.fastq.gz; do
        sample=\$(basename "\$f" -mergedpairs.fastq.gz)
        gunzip -c "\$f" | awk -v s="\$sample" '
            NR%4==1 { print "@" s ";sample=" s ";" substr(\$0,2); next }
            { print }
        '
    done > all_reads.fastq

    # 2. Quality filtering ──
    vsearch --fastq_filter all_reads.fastq \
        --fastq_trunclen ${params.trunclen} \
        --fastaout filtered.fasta \
        --threads ${task.cpus}

    # 3. Dereplicate ──
    vsearch --derep_fulllength filtered.fasta \
        --output derep.fasta \
        --sizeout \
        --minuniquesize 1 \
        --threads ${task.cpus}

    # 4. Denoise with UNOISE3 ──
    vsearch --cluster_unoise derep.fasta \
        --unoise_alpha ${params.unoise_alpha} \
        --minsize ${params.unoise_minsize} \
        --centroids zotus_raw.fasta \
        --threads ${task.cpus}

    # 5. De novo chimera removal (uchime3) ──
    vsearch --uchime3_denovo zotus_raw.fasta \
        --nonchimeras zotus.fasta

    # 6. Assign IDs
    awk '
/^>/ {
    if (seq != "") {
        print ">" seq
        print seq
    }
    seq=""
    next
}
{
    seq = seq \$0
}
END {
    if (seq != "") {
        print ">" seq
        print seq
    }
}
' zotus.fasta > ASV_sequences_tmp.fasta

tr '[:lower:]' '[:upper:]' < ASV_sequences_tmp.fasta > ASV_sequences.fasta

    # 8. Map all reads back to ZOTUs to build an abundance table ──
    vsearch --usearch_global filtered.fasta \
        --db ASV_sequences.fasta \
        --id ${params.unoise_id} \
        --otutabout otu_table_raw.tsv \
        --threads ${task.cpus}

    # 9. Rename header column; ASV IDs are already the row keys ──
    awk 'BEGIN{FS=OFS="\\t"}
        NR==1 { \$1="SeqID"; print; next }
        { print }
    ' otu_table_raw.tsv > asv_table.tsv

    
    """
}
