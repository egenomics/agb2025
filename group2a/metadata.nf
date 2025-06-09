#!/usr/bin/env nextflow

// csvkit installed (pip install csvkit)
// pandas installed (pip install pandas)

nextflow.enable.dsl=2

params.run_id = null

process MERGE_METADATA_MULTIQC {
    publishDir "runs/${params.run_id}/metadata/", mode: 'copy'

    input:
    path metadata_sample
    path multiqc_fastqc

    output:
    path "sample_metadata.csv"

    script:
    """
    tr '\\t' ',' < ${metadata_sample} > metadata.csv

    awk -F'\\t' 'NR==1 {
        for (i = 1; i <= 12; i++) printf "%s%s", \$i, (i==12?"\\n":"\\t")
    }
    NR > 1 {
        for (i = 1; i <= 12; i++) printf "%s%s", \$i, (i==12?"\\n":"\\t")
    }' ${multiqc_fastqc} > multiqc_qc_clean.tsv

    sed 's/\\t/,/g' multiqc_qc_clean.tsv > multiqc_qc_clean.csv

    if [[ ! -s metadata.csv || ! -s multiqc_qc_clean.csv ]]; then
        echo "One or both CSV files are missing or empty!"
        exit 1
    fi

    csvjoin -c "Sample ID",Sample metadata.csv multiqc_qc_clean.csv > sample_metadata.csv
    """
}
process HOMOSAPINENS_CONTAMINATION {
    // This process adds a column to the metadata CSV with the percentage of Homo sapiens contamination based on Kraken2 reports.
    publishDir "runs/${params.run_id}/metadata/", mode: 'copy'

    input:
    path merged_metadata_csv
    path kraken_reports

    output:
    path "sample_metadata.csv"

    script:
    """
    mkdir reports
    cp ${kraken_reports} reports/

    python3 -c "
    import pandas as pd
    import os

    df = pd.read_csv('${merged_metadata_csv}')
    df['Homo_Sapiens_%'] = 0.0

    for fname in os.listdir('reports'):
        if not fname.endswith('.report.txt'):
            continue
        sample_id = fname.split('.kraken2.report.txt')[0]
        percent = 0.0
        found = False
        with open(os.path.join('reports', fname)) as f:
            for line in f:
                if 'Homo sapiens' in line:
                    try:
                        percent = float(line.strip().split('\\t')[0])
                        found = True
                        break
                    except:
                        pass
        # If not found, optionally set to NaN or leave as 0.0
        df.loc[df['Sample ID'] == sample_id, 'Homo_Sapiens_%'] = percent if found else 0.0

    df.to_csv('sample_metadata.csv', index=False)
"
    """
}
process CLASSIFY_QUALITY_SAMPLES {
    publishDir "runs/${params.run_id}/metadata/", mode: 'copy'

    input:
    path sample_metadata_csv

    output:
    path "sample_metadata.csv"

    script:
    """
    python3 -c "
    import pandas as pd
    df = pd.read_csv('${sample_metadata_csv}')
    df['quality_flag'] = df['%GC'].apply(lambda x: 'FAIL' if x < 40 else 'PASS')
    df.to_csv('sample_metadata.csv', index=False)
"
    """
}

workflow {
    def metadata_tsv_path = "runs/${params.run_id}/metadata/sample_metadata.tsv"
    def multiqc_path      = "multiqc_data/multiqc_fastqc.txt"

    // Use `file()` to test existence
    if (file(multiqc_path).exists()) {

        metadata_sample_ch = Channel.fromPath(metadata_tsv_path)
        multiqc_fastqc_ch  = Channel.fromPath(multiqc_path)
        kraken_reports_ch  = Channel.fromPath("runs/${params.run_id}/taxonomy/kraken2/*.kraken2.report.txt", checkIfExists: true)

        merged_csv_ch = MERGE_METADATA_MULTIQC(metadata_sample_ch, multiqc_fastqc_ch)
        human_augmented_ch = HOMOSAPINENS_CONTAMINATION(merged_csv_ch, kraken_reports_ch)
        CLASSIFY_QUALITY_SAMPLES(human_augmented_ch)

    } else {
        log.warn "Skipping metadata merge: multiqc_fastqc.txt not found at ${multiqc_path}"
    }

}


