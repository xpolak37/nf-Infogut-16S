# WMGS Benchmarking Framework Tutorial

**Authors:**  
Petra Polakovicova¹, Alise Ponsero²  
¹ Institute for Clinical and Experimental Medicine, Prague, Czech Republic  
² Core Bioinformatics, Quadram Institute of Biosciences, Norwich, UK

---

## Table of Contents
1. [Installing Nextflow & Singularity](#1-installing-nextflow--singularity)
2. [Setting up the Pipeline](#2-setting-up-the-pipeline)
3. [Editing your Params](#3-editing-your-params)
4. [Running the Pipeline with Test Data](#4-running-the-pipeline-with-test-data)
5. [Looking at Results](#5-looking-at-results)
6. [Validating the Results](#6-validating-the-results)
7. [Possible Problems](#7-possible-problems)

---

## 1. Installing Nextflow & Singularity

The easiest way to get Nextflow and Singularity is to set up a dedicated Conda environment:

```bash
conda create --prefix </path/to/your/new/nf-env/> bioconda::nextflow
conda activate </path/to/your/new/nf-env/>
conda install conda-forge::singularity
```

---

## 2. Setting up the Pipeline

Clone the pipeline repository from GitHub and run the provided setup script. The script will interactively ask for an installation directory, where it will create two subdirectories: `classifiers/` for reference databases and `singularity_cache/` for container images:

```bash
git clone https://github.com/xpolak37/nf-Infogut-16S.git
cd nf-Infogut-16S
bash setup_pipeline.sh
```

You should then see a screen similar to this. Provide your custom installation path when prompted:

```
======================================================================
  16S PROFILING PIPELINE - SETUP
======================================================================

No installation directory provided.

Please enter the full path where you want to install pipeline resources:
(This will create subdirectories: databases/ and singularity_cache/)

Installation directory: </path/to/your/installation/directory/>

```

Press Enter. The installation will take a couple of minutes. At the end, you should see something like this:

```
======================================================================
  SETUP COMPLETE!
======================================================================

[SUCCESS] 2026-03-17 10:52:47 - All resources downloaded and installed successfully

Installation Summary:
  - Installation directory: </path/to/your/installation/directory>
  - Total time: 2 minutes 20 seconds

Resource Locations:
  - Classifiers:      </path/to/your/installation/directory>/classifiers/
  - Containers:       </path/to/your/installation/directory>/singularity_cache/
  - Configuration:    </path/to/your/installation/directory>/pipeline_paths.config
  - Log file:         </path/to/your/installation/directory>/logs/setup_20260317_105233.log
```

---

## 3. Editing your Params

If the installation was successful, open `nextflow.config` and update `</path/to/your/installation/directory/>` in the following parameters:

```bash
params {
    singularity_cache_dir     = '</path/to/your/installation/directory>/singularity_cache'  // Container images cache
    classifiers_dir = '</path/to/your/installation/directory>/classifiers' // path to your classifiers
}
```

---

## 4. Running the Pipeline with Test Data

Now, you are ready to run the pipeline on the prepared testing data, which is located in the "test" directory. It consists of two samples that have been subsampled to run the testing relatively quickly and easily (it should take approximately ~10 minutes to generate the results).

If you followed the steps above and updated `nextflow.config`, run:

```bash
nextflow run main.nf \
        --input test/test_samplesheet.csv \
        --outdir results
```

If everything runs fine, your screen will show something like this:

```
N E X T F L O W   ~  version 25.10.4
===================================
16S PROFILING PIPELINE
===================================
Input samplesheet : test/test_samplesheet.csv
Output directory  : results
===================================
[35/51bd12] FASTQC_RAW (SRR13005986-test)     | 4 of 4 ✔
[44/2eeefd] CUTADAPT (SRR13005987-test)       | 4 of 4 ✔
[b9/76b852] FASTQC_TRIMMED (SRR13005968-test) | 4 of 4 ✔
[d3/03c2a5] DADA2                             | 1 of 1 ✔
[d7/1a70d7] MERGING_READS (1)                 | 4 of 4 ✔
[0b/aa6f7c] QIIME_IMPORT                      | 1 of 1 ✔
[29/9f888b] QIIME_DEBLUR                      | 1 of 1 ✔
[b7/e2ec5e] QIIME_NAIVE_BAYES_DADA2           | 1 of 1 ✔
[40/282495] QIIME_NAIVE_BAYES_DEBLUR          | 1 of 1 ✔
[ca/ca4ecb] QIIME_BLAST_DADA2                 | 1 of 1 ✔
[25/25121e] QIIME_BLAST_DEBLUR                | 1 of 1 ✔
[10/4205d3] IDTAXA_DADA2                      | 1 of 1 ✔
[53/7a34dd] IDTAXA_DEBLUR                     | 1 of 1 ✔
[11/9f632b] ASSIGNTAXONOMY_DADA2              | 1 of 1 ✔
[07/898d86] ASSIGNTAXONOMY_DEBLUR             | 1 of 1 ✔
[2b/b66308] AT_DADA2_METASTANDARD             | 1 of 1 ✔
[2a/026c7e] NB_DADA2_METASTANDARD             | 1 of 1 ✔
[4e/549f43] BLAST_DADA2_METASTANDARD          | 1 of 1 ✔
[57/b705d0] IDTAXA_DADA2_METASTANDARD         | 1 of 1 ✔
[8e/c27e2e] AT_DEBLUR_METASTANDARD            | 1 of 1 ✔
[11/2f31b9] NB_DEBLUR_METASTANDARD            | 1 of 1 ✔
[65/f5f767] BLAST_DEBLUR_METASTANDARD         | 1 of 1 ✔
[03/f26d06] IDTAXA_DEBLUR_METASTANDARD        | 1 of 1 ✔
[1d/138279] MULTIQC                           | 1 of 1 ✔
Pipeline completed at: 2026-03-17T11:07:16.373240639+01:00
Execution status: SUCCESS
Duration: 11m 53s
Output directory: results

Completed at: 17-Mar-2026 11:07:16
Duration    : 11m 53s
CPU hours   : 9.6
Succeeded   : 36
```

---

## 5. Looking at Results

By following the steps in this tutorial, you should see a new directory called results in your working directory `nf-Infogut-16S`. All results are saved there. The most important output is located in the `08\_metastandard` folder, where standardized TSV files from all ASV inference and taxonomic tools created. These files collapse the results to the genus level using a unified format and merge samples that were processed with this pipeline.

```bash
cd results/08_metastandard
ls
```

You will see multiple files:

```
dada2_AssignTaxonomy_run01_genus.tsv  
dada2_NaiveBayes_run01_genus.tsv
dada2_BLAST_run01_genus.tsv
dada2_IDTAXA_run01_genus.tsv
deblur_IDTAXA_run01_genus.tsv
deblur_AssignTaxonomy_run01_genus.tsv  
deblur_NaiveBayes_run01_genus.tsv
deblur_BLAST_run01_genus.tsv
```

---

## 6. Validating the Results

TO BE ADDED

---

## 7. Possible Problems

**No space left on device error:**  
If the pipeline fails with this error, set a custom temp directory before running:

```bash
export SINGULARITY_TMPDIR=</path/to/your/desired/tmp>
```

**Other issues:**  
If you encounter a problem not covered here, please contact us:  
- petra.polakovicova@ikem.cz  
- alise.ponsero@quadram.ac.uk
