#!/bin/bash

ADAPTERS="help/adapters/TruSeq3-PE-2.fa"
IN_DIR="raw_data"
OUT_DIR="fastq_trimmed"
mkdir -p ${OUT_DIR}

# process each sample
for SAMPLE in ERR1328363 ERR1328366 ERR1328415; do
    cutadapt \
        -j 1 \
        -a file:${ADAPTERS} \
        -A file:${ADAPTERS} \
        -q 3,3 \
        -m 36 \
        --pair-filter=any \
        -o ${OUT_DIR}/${SAMPLE}_1P.fastq.gz \
        -p ${OUT_DIR}/${SAMPLE}_2P.fastq.gz \
        ${IN_DIR}/${SAMPLE}_1.fastq.gz \
        ${IN_DIR}/${SAMPLE}_2.fastq.gz
done
