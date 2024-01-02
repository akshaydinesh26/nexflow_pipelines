// make sure required packages are installed in system conda environment with environment name as tool name.
// reads should be present in hqreads folder as output of np_processing.nf workflow.

// the folders associated are given below
// bwa index should be present in the index folder
params.reads="$projectDir/hqreads"
params.index="$projectDir/index"
params.outdir="$projectDir/result"

// #location of the programs programs
params.prgms="/home/akshay/miniconda3/envs/wgs_bwa"

log.info """\

            Quality Control Report Workflow
            Raw Reads          : "${params.reads}"
            index location     : "${params.index}"
            aligner            : "bwa mem"
            output directory   : "${params.outdir}"

"""

process alignment {

  cpus 4
  conda "${params.prgms}"

  input:
  tuple val(sample_id), path(reads)

  output:
  path "bamfile"

  script:
  """
  bwa mem 
  """
}
