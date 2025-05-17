# agb2025
Repository for the AGB 2025 common class project

[![Documentation  Wiki](https://img.shields.io/static/v1?label=Documentation&message=Wiki&labelColor=black&color=blue&logo=github&logoColor=white)](https://github.com/egenomics/agb2025/wiki)


Wiki you can modify and use to document the process!

## Context (Here we will change it with the name of the pipeline)

This is a pipeline created in the Hospital del Mar bioinformaticians group (HdMBioinfo) to process microbiota data coming from stool samples processed in the hospital's analytical facilities by petitions of the Digestology Unit. The main objective is to classify the samples between healthy and non-healthy individuals, so the clinicians afterwards can decide useful interventions (ex: if that individual could be a potential stool donor, or if the individual needs from a stool transplant). The sample metadata comes from the hospital database, we retrieve it using the sampleID provided by the barcodes in the stool samples. Once the sample has been processed by a technician (mainly doing the DNA extraction), the sample goes into an Hi-Seq Illumina Sequencer. As we are performing advanced analysis of the samples, a small team made by bioinformatician technicians do the demultiplexing and send to us the FastQ files from each sample processed in that sequencing run.

## Project Directory Structure

```text
├── 📁 raw_data/            # Original FASTQ files (paired-end reads)
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   └── ...
├── 📁 metadata/            # Sample metadata (CSV or TSV)
│   └── run_metadata.csv
├── 📁 scripts/             # Custom or pipeline scripts
│   └── …
├── 📁 outputs/              # Outputs per run
│   └── 📁 run_YYYYmmdd_hhmmss
│       ├── 📁 qc_reports/          # FastQC or MultiQC reports
│       │   ├── fastqc/
│       │   └── multiqc_report.html
│       ├── 📁 processed_data/      # Trimmed, filtered reads
│       │   ├── data_filtered/
│       │   └── qiime2_demux/
│       ├── 📁 feature_tables/      # OTU or ASV tables
│       │   ├── table.qza
│       │   └── table.tsv
│       ├── 📁 taxonomy/            # Taxonomy assignment results
│       │   ├── classifier.qza
│       │   ├── taxonomy.qza
│       │   └── taxonomy.tsv
│       ├── 📁 phylogeny/           # Phylogenetic tree files
│       │   ├── rooted-tree.qza
│       │   └── unrooted-tree.qza
│       ├── 📁 diversity_analysis/  # Alpha & beta diversity
│       │   ├── core-metrics-results/
│       │   └── emperor_plots/
│       └── 📁 visualizations/      # QIIME 2 visualizations
│           ├── taxonomy_barplots.qzv
│           ├── rarefaction_curves.png
│           └── …
├── 📁 results/             # Final summary tables and plots per run (or sample?)
│   └── summary_report.html
│ 
└── 📁 logs/                # Pipeline logs
    └── run_YYYYmmdd_hhmmss_pipeline.log

