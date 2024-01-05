// make sure required packages are installed in system conda environment with environment name as tool name.
// reads should be present in hqreads folder as output of np_processing.nf workflow.

// the folders associated are given below
// bwa index should be present in the index folder
// name the index prefix as genome
params.reads="$projectDir/hqReads/*_{1,2}.fastq.gz"
params.index="$projectDir/index/*.{amb,ann,pac,bwt,sa}"
params.outdir="$projectDir/result"
params.tempdir="$projectDir/tempdir"

// inputs values required
params.seed=19

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
  val seedlen
  path bwa_index
  
  output:
  path "${sample_id}.bam"

  script:
  def index_name="${bwa_index[0]}";
  def new_name=index_name.replaceAll(".amb","");
  """
  bwa mem -t ${task.cpus} "${new_name}" "${reads[0]}" "${reads[1]}" | samtools view -b -hS -@ ${task.cpus} -o "${sample_id}.bam" -
  """
}


process bam_nsort {
  conda "${params.prgms}"
  cpus 4

  input:
  path bamfile

  output:
  path "${bamfile.baseName}_nsort.bam"
  
  script:
  """
  samtools sort -@ ${task.cpus} -n -o "${bamfile.baseName}_nsort.bam" "${bamfile}"
  """
}


process bam_fixmate {
  conda "${params.prgms}"
  cpus 4

  input:
  path nsortbam
  
  output:
  path "${nsortbam.baseName}_fixmate.bam"
  
  script:
  """
  samtools fixmate "${nsortbam}" "${nsortbam.baseName}_fixmate.bam" 
  """
}


process bam_psort {
  conda "${params.prgms}"
  cpus 4

  input:
  path fixmatebam
  
  output:
  path "${fixmatebam.baseName}_psort.bam"
  
  script:
  """
  samtools sort -@ ${task.cpus} -o "${fixmatebam.baseName}_psort.bam" "${fixmatebam}"
  """
}


process markdupbam {
  publishDir "${params.outdir}/bamfile", mode: 'copy'
  conda "${params.prgms}"
  cpus 4

  input:
  path psortbam

  output:
  path "${psortbam.baseName}_markdup.bam"

  script:
  """
  samtools markdup -@ ${task.cpus} "${psortbam}" "${psortbam.baseName}_markdup.bam"
  """
 }


// process bamsummary {

// }

// process varantCalling {

// }

// process variantSummary {

// }

workflow {
  channel.fromFilePairs(params.reads, checkIfExists: true).set{ read_ch }
  channel.value(params.seed).set{ seed_ch }
  bwa_index = channel.fromPath(params.index, checkIfExists: true).collect()
  align_ch = alignment(read_ch,seed_ch,bwa_index)
  nsort_ch = bam_nsort(align_ch)
  fixmate_ch = bam_fixmate(nsort_ch)
  psort_ch = bam_psort(fixmate_ch)
  markdup_ch = markdupbam(psort_ch)
}