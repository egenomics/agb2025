#!/bin/bash -ue
mkdir -p multiqc_report
multiqc ERR1328415_2_fastqc.zip ERR1328415_2_fastqc.html ERR1328415_1_fastqc.zip ERR1328415_1_fastqc.html ERR1328363_1_fastqc.zip ERR1328363_1_fastqc.html ERR1328366_2_fastqc.zip ERR1328366_2_fastqc.html ERR1328366_1_fastqc.zip ERR1328366_1_fastqc.html ERR1328363_2_fastqc.zip ERR1328363_2_fastqc.html -o multiqc_report
