// runs as nextflow run my_rnaseq.nf -with-docker
// or add docker.enabled = true to nextflow.config file
params.reads = "$projectDir/data/ggal/*_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir  = "$projectDir/my_rnaseq_data/result"

log.info """\

            Simple RNAseq Workflow
            Raw Reads        : "${params.reads}"
            Transcriptome    : "${params.transcriptome_file}"
            multiqc report   : "${params.multiqc}"
            output directory : "${params.outdir}"

"""

/*
 * define the INDEX process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
    cpus 8
    memory '4 GB'
    // https://www.nextflow.io/docs/latest/process.html#directives
    // above is called directives sets resources used by process should be above all

    input:
    path transcriptome

    output:
    path 'salmon_index'
   
    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

// salmon quantification
process QUANTIFICATION {
    cpus 8
    tag "Salmon on $sample_id"
    // publish into a directory
    publishDir "my_rnaseq_data/rnaseq"
    input:
    // you can use any name "sam" also just specify its a file(path)
    path salmon_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

// can quote same terminology as bash if specifiled in input
    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """
}

process FASTQC {
    tag "FASTQC on $sample_id"

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

process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}
/* create channel using set operator and check if file exist using method from "channel factory"
The declaration of a channel can be before the workflow scope or within it. 
As long as it is upstream of the process that requires the specific channel.*/
channel.fromFilePairs(params.reads, checkIfExists: true).set{ read_pairs_ch }
workflow{
    index_ch=INDEX(params.transcriptome_file)
    index_ch.view()
    // prints the a list of prefix and a list of reads files for each prefix
    read_pairs_ch.view()
    quant_ch=QUANTIFICATION(index_ch,read_pairs_ch)
    fastqc_ch = FASTQC(read_pairs_ch)
    /*We only want one task of MultiQC to be executed to produce one report. Therefore, 
    we use the mix channel operator to combine the two channels followed by the collect operator, 
    to return the complete channel contents as a single element.*/
    MULTIQC(quant_ch.mix(fastqc_ch).collect())
}