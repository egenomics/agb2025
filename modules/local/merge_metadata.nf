process MERGE_METADATA_MULTIQC_PROCESS {
    publishDir "runs/${params.run_id}/metadata/", mode: 'copy'

    input:
    path metadata_sample
    path multiqc_fastqc

    output:
    path "sample_metadata.csv"

    script:
       """
    pip install --quiet pandas
    mkdir -p temp
    tr '\\t' ',' < ${metadata_sample} > temp/metadata.csv

    awk -F'\\t' 'NR==1 {
        for (i = 1; i <= 12; i++) printf "%s%s", \$i, (i==12?"\\n":"\\t")
    }
    NR > 1 {
        for (i = 1; i <= 12; i++) printf "%s%s", \$i, (i==12?"\\n":"\\t")
    }' ${multiqc_fastqc} > temp/multiqc_qc_clean.tsv

    sed 's/\\t/,/g' temp/multiqc_qc_clean.tsv > temp/multiqc_qc_clean.csv

    python3 <<EOF
    import pandas as pd
    meta = pd.read_csv('temp/metadata.csv')
    qc = pd.read_csv('temp/multiqc_qc_clean.csv')

    # If the right columns have different names, adjust them below
    merged = pd.merge(meta, qc, left_on='Sample ID', right_on='Sample')

    merged.to_csv('sample_metadata.csv', index=False)
    EOF
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
