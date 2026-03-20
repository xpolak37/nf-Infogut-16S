#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    16S PROFILING PIPELINE
========================================================================================
    FastQC -> Cutadapt -> DADA2 + QIIME- DEBLUR + USEARCH -> Naive bayes + BLAST + IDTAXA + assignTaxonomy + SINTAX -> MultiQC
----------------------------------------------------------------------------------------
*/

// Print pipeline header
log.info """\
    ===================================
    16S PROFILING PIPELINE
    ===================================
    Input samplesheet : ${params.input}
    Output directory  : ${params.outdir}
    ===================================
    """
    .stripIndent()

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { FASTQC as FASTQC_RAW }        from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED }    from './modules/fastqc'
include { MULTIQC }                     from './modules/multiqc'
include { CUTADAPT }                    from './modules/cutadapt'
include { DADA2 }                       from './modules/dada2'
include { MERGING_READS }               from './modules/merging_reads.nf'
include { QIIME_IMPORT; QIIME_DEBLUR }  from './modules/deblur.nf'
include { QIIME_NAIVE_BAYES as QIIME_NAIVE_BAYES_DADA2 }  from './modules/tax_classifiers.nf'
include { QIIME_NAIVE_BAYES as QIIME_NAIVE_BAYES_DEBLUR }  from './modules/tax_classifiers.nf'
include { QIIME_BLAST as QIIME_BLAST_DADA2 }  from './modules/tax_classifiers.nf'
include { QIIME_BLAST as QIIME_BLAST_DEBLUR }  from './modules/tax_classifiers.nf'
include { IDTAXA as IDTAXA_DADA2 }  from './modules/tax_classifiers.nf'
include { IDTAXA as IDTAXA_DEBLUR }  from './modules/tax_classifiers.nf'
include { ASSIGNTAXONOMY as ASSIGNTAXONOMY_DADA2 }  from './modules/tax_classifiers.nf'
include { ASSIGNTAXONOMY as ASSIGNTAXONOMY_DEBLUR }  from './modules/tax_classifiers.nf'
include { METASTANDARD as NB_DADA2_METASTANDARD }        from './modules/MetaStandard16S'
include { METASTANDARD as NB_DEBLUR_METASTANDARD }       from './modules/MetaStandard16S'
include { METASTANDARD as BLAST_DADA2_METASTANDARD }        from './modules/MetaStandard16S'
include { METASTANDARD as BLAST_DEBLUR_METASTANDARD }       from './modules/MetaStandard16S'
include { METASTANDARD as IDTAXA_DADA2_METASTANDARD }        from './modules/MetaStandard16S'
include { METASTANDARD as IDTAXA_DEBLUR_METASTANDARD }       from './modules/MetaStandard16S'
include { METASTANDARD as AT_DADA2_METASTANDARD }        from './modules/MetaStandard16S'
include { METASTANDARD as AT_DEBLUR_METASTANDARD }       from './modules/MetaStandard16S'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
     
    // Read and parse samplesheet
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            def sample_id = row.sample
            def read1 = file(row.read1)
            def read2 = file(row.read2)
            
            // Validate files exist
            if (!read1.exists()) exit 1, "ERROR: Read1 file does not exist: ${read1}"
            if (!read2.exists()) exit 1, "ERROR: Read2 file does not exist: ${read2}"
            
            return tuple(sample_id, read1, read2)
        }
        .set { ch_input_reads }
    
    // FastQC on raw reads
    FASTQC_RAW(ch_input_reads, "raw")
    
    // cutadapt for primer and adapter trimming
    CUTADAPT(ch_input_reads)

    // fastqc TRIMMED
    FASTQC_TRIMMED(CUTADAPT.out.reads, "trimmed")
    
    // DADA2
    dada2_input = CUTADAPT.out.reads
        .collect(flat: false)
        // IMPORTANT: pairing relies on positional alignment
        // ids[i], r1s[i], r2s[i] must correspond to the same sample
        .map { samples -> 
            def sorted = samples.toSorted { sample -> sample[0]}
            [
                ids: sorted.collect { sample -> sample[0] },
                r1s: sorted.collect { sample -> sample[1] },
                r2s: sorted.collect { sample -> sample[2] }
            ]
        }

    DADA2(dada2_input.ids, dada2_input.r1s, dada2_input.r2s)

    // PREPROCESSING FOR DEBLUR AND USEARCH
    MERGING_READS(CUTADAPT.out.reads)

    // DEBLUR
    qiime_input = MERGING_READS.out.reads
    .map { sample_id, reads -> reads }
    .collect()

    QIIME_IMPORT(qiime_input)
    QIIME_DEBLUR(QIIME_IMPORT.out)


    // USEARCH - TO DO


    // TAXONOMIC ASSIGNMENT
    ////////////////////////////////////////////
    // 1 QIIME NAIVE BAYES
    qiime_NB_classifier = file(params.classifiers_dir + "/silva-138.2-ssu-nr99-341F-805R-classifier.qza")

    // 1A DADA2
    QIIME_NAIVE_BAYES_DADA2(DADA2.out.fasta,"dada2",
                            qiime_NB_classifier)

    // 1B DEBLUR
    QIIME_NAIVE_BAYES_DEBLUR(QIIME_DEBLUR.out.fasta,"deblur",
                            qiime_NB_classifier)

    // 1C UNOISE - TO DO 

    ////////////////////////////////////////////
    // 2 QIIME BLAST
    qiime_BLAST_classifier_reads = file(params.classifiers_dir + "/silva-138.2-ssu-nr99-seqs-filt.qza")
    qiime_BLAST_classifier_tax = file(params.classifiers_dir + "/silva-138.2-ssu-nr99-tax.qza")

    // 2A DADA2
    QIIME_BLAST_DADA2(DADA2.out.fasta,"dada2",
                    qiime_BLAST_classifier_reads,
                    qiime_BLAST_classifier_tax)
   
    // 2B DEBLUR
    QIIME_BLAST_DEBLUR(QIIME_DEBLUR.out.fasta,"deblur",
                    qiime_BLAST_classifier_reads,
                    qiime_BLAST_classifier_tax)

    // 2C UNOISE - TO DO 


    ///////////////////////////////////////////////
    // 3 IDTAXA
    idtaxa_classifier = file(params.classifiers_dir + "/idtaxa_trainingSet_V3V4_silva_138_2.RData")

    // 3A DADA2
    IDTAXA_DADA2(DADA2.out.fasta,"dada2",
                idtaxa_classifier)

    // 3B DEBLUR
    IDTAXA_DEBLUR(QIIME_DEBLUR.out.fasta,"deblur",
                idtaxa_classifier)

    // 3C UNOISE - TO DO 

    ///////////////////////////////////////////////
    // 4 ASSIGNTAXONOMY
    assigntaxonomy_classifier = file(params.classifiers_dir + "/dada2_trainset_v3v4_silva_nr99_v138_2.fa.gz")

    // 4A DADA2
    ASSIGNTAXONOMY_DADA2(DADA2.out.fasta,"dada2",
                        assigntaxonomy_classifier)

    // 4B DEBLUR
    ASSIGNTAXONOMY_DEBLUR(QIIME_DEBLUR.out.fasta,"deblur",
                        assigntaxonomy_classifier)

    // 4C UNOISE - TO DO 

    ///////////////////////////////////////////////
    // 5 SINTAX 

    // 5A DADA2

    // 5B DEBLUR

    // 5C SINTAX


    ///////////////////////////////////////////////
    // 6 METASTANDARD 16S
    // 6A DADA2
    AT_DADA2_METASTANDARD(DADA2.out.asv_table,ASSIGNTAXONOMY_DADA2.out,"AssignTaxonomy")
    NB_DADA2_METASTANDARD(DADA2.out.asv_table,QIIME_NAIVE_BAYES_DADA2.out,"NaiveBayes")
    BLAST_DADA2_METASTANDARD(DADA2.out.asv_table,QIIME_BLAST_DADA2.out,"BLAST")
    IDTAXA_DADA2_METASTANDARD(DADA2.out.asv_table,IDTAXA_DADA2.out,"IDTAXA")

    // 6B DEBLUR
    AT_DEBLUR_METASTANDARD(QIIME_DEBLUR.out.asv_table,ASSIGNTAXONOMY_DEBLUR.out,"AssignTaxonomy")
    NB_DEBLUR_METASTANDARD(QIIME_DEBLUR.out.asv_table,QIIME_NAIVE_BAYES_DEBLUR.out,"NaiveBayes")
    BLAST_DEBLUR_METASTANDARD(QIIME_DEBLUR.out.asv_table,QIIME_BLAST_DEBLUR.out,"BLAST")
    IDTAXA_DEBLUR_METASTANDARD(QIIME_DEBLUR.out.asv_table,IDTAXA_DEBLUR.out,"IDTAXA")

    // 6C UNOISE

    // Collect all QC files for MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect().ifEmpty([]))
    
    // MultiQC aggregation
    MULTIQC(ch_multiqc_files.collect())
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """\
        Pipeline completed at: ${workflow.complete}
        Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Duration: ${workflow.duration}
        Output directory: ${params.outdir}
        """
        .stripIndent()
}