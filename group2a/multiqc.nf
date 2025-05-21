#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/egenomics/agb2025.git
    Slack  : https://agb2025.slack.com/archives/C08QVHBFE6A
----------------------------------------------------------------------------------------
    QC AND PREPROCESSING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.raw_data = "../raw_data/*.fastq.gz"

process fastQC {
    tag "FastQC"

    conda 'bioconda::fastqc=0.12.1'

    input:
    path reads

    output:
    tuple path("*.zip"), path("*.html")


    script:
    """
    fastqc ${reads}
    """
}

params.out_dir  = "multiqc_report"

process multiQC {
    tag "MultiQC"
    conda 'bioconda::multiqc=1.14'

    input:
    path raw_data

    output:
    path("multiqc_report", type: 'dir')

    script:
    """
    mkdir -p ${params.out_dir}
    multiqc ${raw_data} -o ${params.out_dir}
    """
}

workflow {
    raw_data_ch = Channel.fromPath(params.raw_data)

    fastqc_out = fastQC(raw_data_ch)

    // Flatten the tuples to get all files into a single channel
    fastqc_files = fastqc_out.flatMap { zip, html -> [zip, html] }

    // Collect all FastQC outputs into one list
    fastqc_files.collect().set { all_fastqc_outputs }

    multiQC(all_fastqc_outputs)
}
// End of the script
// This script processes raw sequencing data using FastQC for quality control and MultiQC for report generation.
// It takes raw data files from the specified directory, runs FastQC on each file, and then generates a MultiQC report