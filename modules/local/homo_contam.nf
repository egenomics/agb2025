process HOMO_CONTAM_PROCESS {
    // Use $RUN_ID inside publishDir
    publishDir "runs/${params.run_id}/metadata/", mode: 'copy'

    input:
    path merged_metadata_csv
    path kraken_reports

    output:
    path "metadata.tsv"


    script:
    """

    mkdir -p runs/${params.run_id}/metadata/reports
    cp ${kraken_reports} runs/${params.run_id}/metadata/reports/

    python3 <<EOF
    import pandas as pd
    import os

    df = pd.read_csv('${merged_metadata_csv}', sep='\\t')
    df['Homo_Sapiens_%'] = 0.0

    for fname in os.listdir('runs/${params.run_id}/metadata/reports'):
        if not fname.endswith('.report.txt'):
            continue
        sample_id = fname.split('.kraken2.report.txt')[0]
        percent = 0.0
        found = False
        with open(os.path.join('runs/${params.run_id}/metadata/reports', fname)) as f:
            for line in f:
                if 'Homo sapiens' in line:
                    try:
                        percent = float(line.strip().split('\\t')[0])
                        found = True
                        break
                    except:
                        pass
        print(df.columns)
        df.loc[df["Sample ID"] == sample_id, 'Homo_Sapiens_%'] = percent if found else 0.0
    print(df.head())
    df.to_csv('metadata.tsv', index=False)
    EOF
    """
}


workflow HOMOSAPINENS_CONTAMINATION {
    take:
    merged_metadata_csv
    kraken_reports

    main:
    HOMO_CONTAM_PROCESS(merged_metadata_csv, kraken_reports)

    emit:
    HOMO_CONTAM_PROCESS.out
}
