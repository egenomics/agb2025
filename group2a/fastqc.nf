#!usr/bin/env nextflow

params.raw_data = "../raw_data/*.fastq.gz"

process fastQC {
    tag "FastQC on ${reads.getBaseName()}"

    conda 'bioconda::fastqc=0.12.1'

    input:
    path reads

    output:
    path "*.zip"
    path "*.html"

    script:
    """
    fastqc ${reads}
    """
}

workflow {
    raw_data_ch = Channel.fromPath(params.raw_data)

    fastqc_out = fastQC(raw_data_ch)
    multiQC(fastqc_out)
}
