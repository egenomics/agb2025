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
    # extract + merge
    awk -F'\\t' 'NR==1 {
        for (i = 1; i <= 12; i++) {
            printf "%s%s", \$i, (i==12 ? "\\n" : "\\t")
        }
    }
    NR > 1 {
        for (i = 1; i <= 12; i++) {
            printf "%s%s", \$i, (i==12 ? "\\n" : "\\t")
        }
    }' ${multiqc_fastqc} > multiqc_qc_clean.tsv

    sed 's/\\t/,/g' multiqc_qc_clean.tsv > multiqc_qc_clean.csv

    csvjoin -c Sample_ID,Sample ${metadata_sample} multiqc_qc_clean.csv > metadata_sample_merged.csv
    """
}

workflow {
    metadata_sample_ch   = Channel.fromPath("../metadata/metadata_sample.csv")
    multiqc_fastqc_ch    = Channel.fromPath("multiqc_data/multiqc_fastqc.txt")

    merge_metadata(metadata_sample_ch, multiqc_fastqc_ch)
}
