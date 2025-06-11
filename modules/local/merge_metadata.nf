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
    ls -lh ${metadata_sample} ${multiqc_fastqc}
    head -n 3 ${metadata_sample}

    tr '\\t' ',' < ${metadata_sample} > temp/metadata.csv

    awk -F'\\t' 'NR==1 {
        for (i = 1; i <= 15; i++) printf "%s%s", \$i, (i==15?"\\n":"\\t")
    }
    NR > 1 {
        for (i = 1; i <= 15; i++) printf "%s%s", \$i, (i==15?"\\n":"\\t")
    }' ${multiqc_fastqc} > temp/multiqc_qc_clean.tsv

    sed 's/\\t/,/g' temp/multiqc_qc_clean.tsv > temp/multiqc_qc_clean.csv

    csvjoin -c "Sample ID",Sample temp/metadata.csv temp/multiqc_qc_clean.csv > metadata.csv

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
