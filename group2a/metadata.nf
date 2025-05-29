#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process merge_metadata {
    publishDir "multiqc_data", mode: 'copy'
    
    input:
    path metadata_sample
    path multiqc_fastqc

    output:
    path "metadata_sample_merged.csv"

    script:
    """
    # Get number of columns in metadata header
    num_cols=\$(head -n1 ${metadata_sample} | sed 's/[^,]//g' | wc -c)

    # Create header in output
    cp ${metadata_sample} metadata_sample_merged.csv

    # Extract first 12 columns from multiqc and convert to CSV
    awk -F'\\t' 'NR > 1 {
        printf "%s", \$1
        for (i=2; i<=12; i++) printf ",%s", \$i
        for (j=13; j<=25; j++) printf ","
        printf "\\n"
    }' ${multiqc_fastqc} >> metadata_sample_merged.csv
"""
}

workflow {
    metadata_sample_ch   = Channel.fromPath("../metadata/metadata_sample.csv")
    multiqc_fastqc_ch    = Channel.fromPath("multiqc_data/multiqc_fastqc.txt")

    merge_metadata(metadata_sample_ch, multiqc_fastqc_ch)
}
