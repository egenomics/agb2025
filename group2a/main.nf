#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/egenomics/agb2025.git
    Slack  : https://agb2025.slack.com/archives/C08QVHBFE6A
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QC AND PREPROCESSING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process multiQC {

    input:
    path raw_data

    output:
    path outputs

    script:
    """
    multiqc ${raw_data} -o ${outputs}
    """
}

workflow {
    // Define the input data
    raw_data = file("/path/to/raw/data")

    // Define the output directory
    outputs = file("/")

    // Run the multiQC process
    multiQC(raw_data: raw_data, outputs: outputs)
}