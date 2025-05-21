#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/egenomics/agb2025.git
    Slack  : https://agb2025.slack.com/archives/C08QVHBFE6A
----------------------------------------------------------------------------------------
    QC AND PREPROCESSING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.raw_data = "../raw_data"
params.out_dir  = "multiqc_report"

process multiQC {
    tag "MultiQC on ${raw_data}"

    input:
    path raw_data

    output:
    path("${params.out_dir}")

    script:
    """
    mkdir -p ${params.out_dir}
    multiqc ${raw_data} -o ${params.out_dir}
    """
}

workflow {
    raw_data = file(params.raw_data)

    multiQC(raw_data)
}
