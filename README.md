# agb2025
Repository for the AGB 2025 common class project

[![Documentation  Wiki](https://img.shields.io/static/v1?label=Documentation&message=Wiki&labelColor=black&color=blue&logo=github&logoColor=white)](https://github.com/egenomics/agb2025/wiki)


Wiki you can modify and use to document the process!

## Documentation

In this wiki page you will find the information about the pipeline context, the sample processing and the decisions made through each of the modules (https://github.com/egenomics/agb2025/wiki/Pipeline-Context).

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
```

## Quick start (Docker edition)

### 1 · Prerequisites

| Tool               | macOS (Homebrew)                                                            | Ubuntu / Debian                                                           | Notes                                                         |
|--------------------|-----------------------------------------------------------------------------|---------------------------------------------------------------------------|---------------------------------------------------------------|
| **Docker ≥ 24**    | `brew install --cask docker`<br/>Launch *Docker Desktop*                     | `sudo apt install docker.io`                                              | Ensure the Docker daemon is running and your user has access. |
| **Nextflow ≥ 23.10** | `brew install nextflow`                                                    | `curl -s https://get.nextflow.io \| bash && sudo mv nextflow /usr/local/bin/` | No extra software needed – the pipeline pulls everything in containers. |


### 2 · Run the pipeline
```bash
# minimal
nextflow run main.nf -profile docker

# resume an interrupted run (skips completed tasks)
nextflow run main.nf -profile docker -resume