#!/usr/bin/env nextflow

// csvkit installed (pip install csvkit)
// pandas installed (pip install pandas)

nextflow.enable.dsl=2

params.run_id = null

process merge_metadata {
    publishDir "multiqc_data", mode: 'copy'

    input:
    path metadata_sample
    path multiqc_fastqc

    output:
    path "sample_metadata.csv"

    script:
    """
    # Convert TSV to CSV
    tr '\\t' ',' < ${metadata_sample} > metadata.csv

    # Extract first 12 columns from FastQC summary
    awk -F'\\t' 'NR==1 {
        for (i = 1; i <= 12; i++) printf "%s%s", \$i, (i==12?"\\n":"\\t")
    }
    NR > 1 {
        for (i = 1; i <= 12; i++) printf "%s%s", \$i, (i==12?"\\n":"\\t")
    }' ${multiqc_fastqc} > multiqc_qc_clean.tsv

    sed 's/\\t/,/g' multiqc_qc_clean.tsv > multiqc_qc_clean.csv

    # Validate inputs
    if [[ ! -s metadata.csv || ! -s multiqc_qc_clean.csv ]]; then
        echo "One or both CSV files are missing or empty!"
        exit 1
    fi

    # Merge on common columns (adjust as needed)
    csvjoin -c Sample_ID,Sample metadata.csv multiqc_qc_clean.csv > sample_metadata.csv
    """
}

process classify_quality {
    publishDir "multiqc_data", mode: 'copy'

    input:
    path sample_metadata_csv

    output:
    path "sample_metadata_classified.csv"

    script:
    """
    python3 -c "
import pandas as pd
df = pd.read_csv('${sample_metadata_csv}')
df['quality_flag'] = df['%GC'].apply(lambda x: 'FAIL' if x < 40 else 'PASS')  # adjust logic
df.to_csv('sample_metadata_classified.csv', index=False)
"
    """
}


workflow {

    if (params.run_id == null) {
        error " Missing required parameter: --run_id"
    }

    // Define paths
    def metadata_tsv_path = "runs/${params.run_id}/metadata/sample_metadata.tsv"
    def multiqc_path      = "multiqc_data/multiqc_fastqc.txt"

    // Channels
    metadata_sample_ch = Channel.fromPath(metadata_tsv_path)
    multiqc_fastqc_ch  = Channel.fromPath(multiqc_path)

    // Processes
    merge_metadata(metadata_sample_ch, multiqc_fastqc_ch)

    classify_quality("sample_metadata.csv")
}

