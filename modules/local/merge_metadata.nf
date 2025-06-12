process MERGE_METADATA_MULTIQC_PROCESS {
    publishDir "runs/${params.run_id}/metadata/", mode: 'copy'

    input:
    path metadata_sample
    path multiqc_fastqc

    output:
    path "metadata.tsv"

    script:
    """
    python3 <<EOF
    import pandas as pd

    metadata = pd.read_csv('${metadata_sample}', sep='\t')
    qc_data = pd.read_csv('${multiqc_fastqc}', sep='\t')

    qc_data['Sample_base'] = qc_data['Sample'].str.extract(r'^(ERR\\d+)')

    # Compute average %GC per sample (based on all 4 lines per sample)
    gc_avg = qc_data.groupby('Sample_base')['%GC'].mean().reset_index()

    # Perform merge on the "Sample ID" and "Sample" columns
    metadata['Sample_base'] = metadata['Sample ID'].str.extract(r'^(ERR\\d+)')

    merged = metadata.merge(gc_avg, on='Sample_base', how='left')

    merged.drop(columns=['Sample_base'], inplace=True)

    # Save the merged data to a new TSV file
    merged.to_csv('metadata.tsv', sep='\t', index=False)

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
