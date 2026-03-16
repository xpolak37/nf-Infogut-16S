process QIIME_NAIVE_BAYES {
    publishDir "${params.outdir}/07a_qiime_naive_bayes", mode: 'copy'

    input:
    path(asv_fasta)
    val denoising_tool  // 'dada2' or 'deblur' or 'unoise'

    output:
    path("${denoising_tool}_taxa_table.tsv")

    script:
    """

    export NUMBA_CACHE_DIR=\${PWD}/numba_cache
    mkdir -p \${PWD}/numba_cache

    export TMPDIR=\${PWD}/tmp
    mkdir -p \${PWD}/tmp

    # IMPORT
    qiime tools import \
    --input-path ${asv_fasta} \
    --output-path rep-seqs.qza \
    --type 'FeatureData[Sequence]'

    # TAXONOMY
    qiime feature-classifier classify-sklearn \\
        --i-reads rep-seqs.qza \\
        --i-classifier ${params.classifiers_dir}/silva-138.2-ssu-nr99-341F-805R-classifier.qza \\
        --p-n-jobs ${task.cpus} \\
        --p-confidence ${params.qiime_naive_bayes_confidence} \\
        --o-classification taxonomy.qza 

    qiime tools export --input-path taxonomy.qza --output-path .
    
    python3 - <<'EOF'
import csv

input_file = "taxonomy.tsv"
output_file = "${denoising_tool}_taxa_table.tsv"

with open(input_file, newline='') as fin, \
     open(output_file, "w", newline='') as fout:

    reader = csv.reader(fin, delimiter="\\t")
    writer = csv.writer(fout, delimiter="\\t")

    header = next(reader)
    header[0] = "SeqID"
    header[1] = "Taxonomy"
    header[2] = "Confidence"

    writer.writerow(header)

    for row in reader:
        writer.writerow(row)
EOF
    """
} 

process QIIME_BLAST {
    publishDir "${params.outdir}/07b_qiime_blast", mode: 'copy'

    input:
    path(asv_fasta)
    val denoising_tool  // 'dada2' or 'deblur' or 'unoise'

    output:
    path("${denoising_tool}_taxa_table.tsv")

    script:
    """
    export NUMBA_CACHE_DIR=\${PWD}/numba_cache
    mkdir -p \${PWD}/numba_cache

    export TMPDIR=\${PWD}/tmp
    mkdir -p \${PWD}/tmp
    
    # IMPORT
    qiime tools import \
    --input-path ${asv_fasta} \
    --output-path rep-seqs.qza \
    --type 'FeatureData[Sequence]'

    # TAXONOMY
    qiime feature-classifier classify-consensus-blast  \\
        --i-query rep-seqs.qza \\
        --i-reference-reads ${params.classifiers_dir}/silva-138.2-ssu-nr99-seqs-filt.qza \\
        --i-reference-taxonomy ${params.classifiers_dir}/silva-138.2-ssu-nr99-tax.qza \\
        --p-num-threads ${task.cpus} \\
        --p-perc-identity ${params.blast_percidentity} \\
        --p-strand both \\
        --p-min-consensus ${params.blast_minconsensus} \\
        --o-classification  taxonomy.qza \\
        --o-search-results blast_results.qza

    qiime tools export --input-path taxonomy.qza --output-path .

        python3 - <<'EOF'
import csv

input_file = "taxonomy.tsv"
output_file = "${denoising_tool}_taxa_table.tsv"

with open(input_file, newline='') as fin, \
     open(output_file, "w", newline='') as fout:

    reader = csv.reader(fin, delimiter="\\t")
    writer = csv.writer(fout, delimiter="\\t")

    header = next(reader)
    header[0] = "SeqID"
    header[1] = "Taxonomy"
    header[2] = "Confidence"
    
    writer.writerow(header)

    for row in reader:
        writer.writerow(row)
EOF
    """
} 

process ASSIGNTAXONOMY {
    publishDir "${params.outdir}/07c_assigntaxonomy", mode: 'copy'

    input:
    path(asv_fasta)
    val denoising_tool  // 'dada2' or 'deblur' or 'unoise'

    output:
    path("${denoising_tool}_taxa_table.tsv"), emit: taxa_table

    script:
    """
    Rscript ${projectDir}/bin/assigntaxonomy.R \
        ${asv_fasta} \
       ${params.classifiers_dir}/dada2_trainset_v3v4_silva_nr99_v138_2.fa.gz \
        ${task.cpus} \
        ${denoising_tool}
    """
}

process IDTAXA {
    publishDir "${params.outdir}/07d_idtaxa", mode: 'copy'

    input:
    path(asv_fasta)
    val denoising_tool  // 'dada2' or 'deblur' or 'unoise'

    output:
    path("${denoising_tool}_taxa_table.tsv")

    script:
    """
    Rscript ${projectDir}/bin/idtaxa.R \
        ${asv_fasta} \
        ${params.classifiers_dir}/idtaxa_trainingSet_V3V4_silva_138_2.RData \
        ${task.cpus} \
        ${denoising_tool}
    """
}

process SINTAX {
    publishDir "${params.outdir}/07e_sintax", mode: 'copy'

    input:
    path(asv_table)
    val denoising_tool  // 'dada2' or 'deblur' or 'unoise'

    output:
    path("${denoising_tool}_taxa_table.tsv")

    script:
    """
    # TO BE IMPLEMENTED
    """
} 