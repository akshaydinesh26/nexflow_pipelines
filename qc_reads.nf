#!/usr/local/bin/env nextflow

// quality control of rawreads with reporting in multiqc format and fastqc format html files.
// reads should be in rawReads folder in basedir(where thi script is present) and in format "name_[1/2].fq.gz"
//  fastqc and multiqc should be in conda environment.

params.reads="$projectDir/rawReads/*_{1,2}.fq.gz"
params.outdir="$projectDir/qcReport"

// #location of programs please update
params.fqc="/home/akshay/miniconda3/envs/fastqc"
params.mqc="/home/akshay/miniconda3/envs/multiqc"

log.info """\

            Quality Control Report Workflow
            Raw Reads          : "${params.reads}"
            fastqc location    : "${params.fqc}"
            multiqc location   : "${params.mqc}"
            output directory   : "${params.outdir}"

"""

process fastqc_run {
  
  conda "${params.fqc}"
  publishDir "${params.outdir}/fastqc", mode: 'copy'

  input:
  tuple val(sample_id), path(reads)

  output:
  path "fastqc_${sample_id}_logs"

  script:
  """
  mkdir fastqc_${sample_id}_logs
  fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
  """
}


process multiqc_run {
  
  conda "${params.mqc}"
  publishDir "${params.outdir}/multiqc", mode:'copy'

  input:
  path '*'

  output:
  path 'multiqc_report.html'

  script:
  """
  multiqc .
  """
}

workflow {
  channel.fromFilePairs(params.reads, checkIfExists: true).set{ reads_ch }
  fastqc_ch = fastqc_run(reads_ch)
  multiqc_run(fastqc_ch.collect())
}