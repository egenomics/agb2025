#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/egenomics/agb2025.git
    Slack  : https://agb2025.slack.com/archives/C08QVHBFE6A
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QC AND PREPROCESSING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process multiQC {
    inuput:
    path raw_data

    output:
    path fastq_trimmed

    script:
    " FASTQ="SRR11576980"

OUT_DIR="fastq_trimmed"

TRIMMOMATIC_PATH="/Trimmomatic-0.39/trimmomatic-0.39.jar"

mkdir -p ${OUT_DIR}

java -jar ${TRIMMOMATIC_PATH} PE \
    -threads 1 \
    fastq/${FASTQ}_1.fastq.gz \
    fastq/${FASTQ}_2.fastq.gz \
    ${OUT_DIR}/${FASTQ}tr_1P.fastq.gz \
    ${OUT_DIR}/${FASTQ}tr_1U.fastq.gz \
    ${OUT_DIR}/${FASTQ}tr_2P.fastq.gz \
    ${OUT_DIR}/${FASTQ}tr_2U.fastq.gz \
    ILLUMINACLIP:help/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
}
