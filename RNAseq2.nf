// runs as nextflow run my_rnaseq.nf -with-docker
// or add docker.enabled = true to nextflow.config file
params.reads = "$projectDir/test_data/*_{1,2}.fq"
params.transcriptome_file = "$projectDir/test_data/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir  = "$projectDir/rnaseq1_data/results"

// # location of programs
params.klsto="~/miniconda3/envs/kallisto"
params.fqc="~/miniconda3/envs/fastqc"
params.mqc="~/miniconda3/envs/multiqc"


log.info """\

            Simple RNAseq Workflow
            Raw Reads         : "${params.reads}"
            Transcriptome     : "${params.transcriptome_file}"
            multiqc report    : "${params.multiqc}"
            output directory  : "${params.outdir}"
            kallisto locattion: "${params.klsto}"
            fastqc location   : "${params.fqc}"
            multiqc location  : "${params.mqc}"

"""

/*
 * define the INDEX process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
    conda "${params.klsto}"
    // https://www.nextflow.io/docs/latest/process.html#directives
    // above is called directives sets resources used by process should be above all

    input:
    path transcriptome

    output:
    path 'kallisto_index'
   
    script:
    """
    kallisto index -i kallisto_index $transcriptome
    """
}

// salmon quantification
process QUANTIFICATION {
    conda "${params.fqc}"
    tag "kallisto on $sample_id"
    // publish into a directory
    publishDir "my_rnaseq_data/rnaseq"
    input:
    // you can use any name "sam" also just specify its a file(path)
    path kallisto_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

// can quote same terminology as bash if specifiled in input
    script:
    """
    kallisto quant --threads $task.cpus -i $kallisto_index -o $sample_id ${reads[0]} ${reads[1]}
    """
}

process FASTQC {
    conda "${params.mqc}"
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
    conda '/home/akshay/miniconda3/envs/multiqc'
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
