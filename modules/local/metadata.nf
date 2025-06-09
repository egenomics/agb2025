// modules/local/metadata.nf

// ========= MERGE METADATA =========
// File: modules/local/metadata.nf

// ========== PROCESS 1 ==========
process MERGE_METADATA_MULTIQC_PROCESS {
    input:
    path metadata_sample
    path multiqc_fastqc

    output:
    path "sample_metadata.csv"

    script:
    """
    mkdir -p runs/temp
    tr '\\t' ',' < ${metadata_sample} > runs/temp/metadata.csv

    awk -F'\\t' 'NR==1 {
        for (i = 1; i <= 12; i++) printf "%s%s", \$i, (i==12?"\\n":"\\t")
    }
    NR > 1 {
        for (i = 1; i <= 12; i++) printf "%s%s", \$i, (i==12?"\\n":"\\t")
    }' ${multiqc_fastqc} > runs/temp/multiqc_qc_clean.tsv

    sed 's/\\t/,/g' runs/temp/multiqc_qc_clean.tsv > runs/temp/multiqc_qc_clean.csv

    csvjoin -c "Sample ID",Sample runs/temp/metadata.csv runs/temp/multiqc_qc_clean.csv > sample_metadata.csv
    """
}

// ========== WORKFLOW WRAPPER ==========
workflow MERGE_METADATA_MULTIQC {
    take:
    metadata_sample
    multiqc_fastqc

    main:
    MERGE_METADATA_MULTIQC_PROCESS(metadata_sample, multiqc_fastqc)

    emit:
    merged_csv_ch = MERGE_METADATA_MULTIQC_PROCESS.out
}

//
// ========= ADD HOMO SAPIENS CONTAMINATION COLUMN =========
process HOMO_CONTAM_PROCESS {
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
    df.loc[df['Sample ID'] == sample_id, 'Homo_Sapiens_%'] = percent if found else 0.0

df.to_csv('sample_metadata.csv', index=False)
"
    """
}

workflow HOMOSAPINENS_CONTAMINATION {
    take:
        merged_metadata_csv
        kraken_reports

    main:
        HOMO_CONTAM_PROCESS(merged_metadata_csv, kraken_reports)

    emit:
        human_augmented_csv -> HOMO_CONTAM_PROCESS.out
}

//
// ========= CLASSIFY SAMPLE QUALITY =========
process CLASSIFY_QUALITY_PROCESS {
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

workflow CLASSIFY_QUALITY_SAMPLES {
    take:
        sample_metadata_csv

    main:
        CLASSIFY_QUALITY_PROCESS(sample_metadata_csv)

    emit:
        classified_csv -> CLASSIFY_QUALITY_PROCESS.out
}
