process CLASSIFY_QUALITY_PROCESS {
    publishDir "runs/${params.run_id}/metadata/", mode: 'copy'

    input:
    path sample_metadata_csv

    output:
    path "metadata.tsv"

    script:
    """
    export RUN_ID=${params.run_id}

    python3 << EOF
    import pandas as pd
    df = pd.read_csv('${sample_metadata_csv}')
    df['quality_flag'] = df['%GC'].apply(lambda x: 'FAIL' if x < 40 else 'PASS')
    df.to_csv('metadata.tsv', index=False)
    EOF

    echo "Sample metadata classified and saved as metadata.tsv in runs/${RUN_ID}/metadata/"
    """

}


workflow CLASSIFY_QUALITY_SAMPLES {
    take:
    sample_metadata_csv

    main:
    CLASSIFY_QUALITY_PROCESS(sample_metadata_csv)

    emit:
    CLASSIFY_QUALITY_PROCESS.out
}
