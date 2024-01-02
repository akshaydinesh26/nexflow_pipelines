#!/usr/local/bin/env nextflow

// quality control of rawreads with reporting in multiqc format and fastqc format html files.
// reads should be in rawReads folder in basedir(where thi script is present) and in format "name_[1/2].fq.gz"
// processed reads should be in "hqReads" folder
//  fastqc and multiqc should be in conda environment.

params.reads="$projectDir/rawReads/*_{1,2}.fastq.gz"
params.outdir="$projectDir/hqReads"
params.qtrim=20
params.len=75

// #location of programs please update
params.fastp="/home/akshay/miniconda3/envs/fastp"

log.info """\

            QC Reads Processing Workflow
            Raw Reads          : "${params.reads}"
            fastp location     : "${params.fastp}"
            output directory   : "${params.outdir}"
	    quality filter     : "${params.qtrim}"
            minimum length     : "${params.len}"

"""

process processRun {
  cpus 8
  conda "${params.fastp}"
  publishDir "${params.outdir}", mode: 'copy', pattern: "*.fastq.gz"
  publishDir "${params.outdir}/report", mode: 'copy', pattern: "*.{html,json}"
  tag "fastp on $sample_id"

  input:
  tuple val(sample_id), path(reads)
  val qtrim
  val lenl

  
  output:
   stdout
  path "${sample_id}.html"
  path "${sample_id}.json"
  path "${sample_id}_1.fastq.gz"
  path "${sample_id}_2.fastq.gz"

  script:
  """
  fastp -i ${reads[0]} -I ${reads[1]} -o "${sample_id}_1.fastq.gz" -O "${sample_id}_2.fastq.gz" -q ${qtrim} -l ${lenl} -R "${sample_id}" -j "${sample_id}.json" -h "${sample_id}.html"
  """
}

workflow {
  channel.fromFilePairs(params.reads, checkIfExists: true).set{ reads_ch }
  channel.value(params.qtrim).set{ qtrim_ch }
  channel.value(params.len).set{ lenl_ch }
  fastp_ch = processRun(reads_ch,qtrim_ch,lenl_ch)
}