
# *** Nextflow Bioinfomatics Pipelines ***

## **before using**

These are some of my rules before automating pipelines,

1 - Only automate repetitive tasks.
2 - Automation should not interfere with robust Quality control of the analysis or the results.
3 - its automation, dont ask too much input in commandline

## **RNAseq1.nf**

Nextflow pipeline for RNAseq data analysis which uses salmon for transcriptome based quantification. Output is salmon output with multiqc report.
make sure all the programs are found in docker container

## **RNAseq2.nf**

Nextflow pipeline for RNAseq data analysis which uses kallisto for transcriptome based quantification. Output is kallisto output with multiqc report.

## Read **Quality Control**

1 - qc_reads.nf    - pre-qc reports
2 - qc_processing.nf - trimming, adapter removal and filtering of reads.

### **qc_reads.nf**

automating generation of fastqc and multiqc reports.
makesure reads are in Project Directory in a folder named rawReads with naming format of "samplename_[1/2].fq.gz"
programs should be in seperate conda environments

_before executing_

if miniconda3 folder is not in home folder make sure to update the location of tools.
update the conda environment names
results will be in qcReport folder in project directory.

### **qc_processing.nf**

