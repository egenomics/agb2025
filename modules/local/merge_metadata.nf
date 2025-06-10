process MERGE_METADATA_MULTIQC_PROCESS {
    publishDir "runs/${params.run_id}/metadata/", mode: 'copy'

    input:
    path metadata_sample
    path multiqc_fastqc

    output:
    path "metadata.tsv"

    script:
       """
    mkdir -p temp

    tr '\\t' ',' < ${metadata_sample} > temp/metadata.csv

    awk -F'\\t' 'NR==1 {
        for (i = 1; i <= 12; i++) printf "%s%s", \$i, (i==12?"\\n":"\\t")
    }
    NR > 1 {
        for (i = 1; i <= 12; i++) printf "%s%s", \$i, (i==12?"\\n":"\\t")
    }' ${multiqc_fastqc} > temp/multiqc_qc_clean.tsv

    sed 's/\\t/,/g' temp/multiqc_qc_clean.tsv > temp/multiqc_qc_clean.csv

    csvjoin -c Sample_ID,Sample temp/metadata.csv temp/multiqc_qc_clean.csv > metadata.csv

    # Convert joined CSV to TSV
    tr ',' '\\t' < temp/metadata.csv > metadata.tsv

    echo "Merged metadata saved as metadata.tsv"
    """
}

workflow MERGE_METADATA_MULTIQC {
    take:
    metadata_sample
    multiqc_fastqc

    main:
    MERGE_METADATA_MULTIQC_PROCESS(metadata_sample, multiqc_fastqc)

    emit:
    MERGE_METADATA_MULTIQC_PROCESS.out
}
