# agb2025
Repository for the AGB 2025 common class project

[![Documentation  Wiki](https://img.shields.io/static/v1?label=Documentation&message=Wiki&labelColor=black&color=blue&logo=github&logoColor=white)](https://github.com/egenomics/agb2025/wiki)


Wiki you can modify and use to document the process!

## Documentation

In this wiki page you will find the information about the pipeline context, the sample processing and the decisions made through each of the modules (https://github.com/egenomics/agb2025/wiki/Pipeline-Context).

## Project Directory Structure

```text
HdMBioinfo-MicrobiotaPipeline/
│
├── bin/                         # Helper scripts and small executables
│   ├── trim_reads.sh
│   ├── assign_taxonomy.py
│   └── generate_qc_report.R
│
├── conf/                        # Configuration files
│   ├── base.config              # Default config
│   ├── hospital.config          # Custom config for Hospital del Mar environment
│   └── docker.config            # Containerized setup (if needed)
│
├── docs/                        # Documentation and metadata templates
│   ├── metadata_template.yaml
│   ├── versioning_policy.md
│   └── pipeline_overview.png
│
├── workflows/                   # Main Nextflow scripts
│   └── main.nf                  # Entry point of the pipeline
│
├── modules/                     # DSL2-style modules for each pipeline step
│   ├── preprocessing/
│   │   └── trim_reads.nf
│   ├── taxonomy/
│   │   └── assign_taxonomy.nf
│   ├── classification/
│   │   └── classify_health_status.nf
│   └── qc/
│       └── fastqc.nf
│
├── raw_data
│   ├── run1
│   │    ├── sample1.fastq.gz
│   │    ├── sample2.fastq.gz
│   │    ├── sample3.fastq.gz
│   │    ├── metadata_technical.csv
│   │    └── metadata_samples.csv
│   ├── run2
│   └── run3
│
├── metadata                      # Merged technical and sample metadata
│   ├── run1
│   │    └── metadata_run1.csv
│   ├── run2
│   └── run3
│
├── outputs/
│   └── run_<run_id>/
│       ├── raw_data/                     # Raw fastq input files (paired-end)
│       ├── processed_data/
│       │   └── data_filtered/            # Filtered and trimmed fastq's
│       ├── qc_reports/
│       │   ├── fastqc/                   # Individual fastqc reports
│       │   └── multiqc_report.html       # Combined multiqc report
│       ├── metadata_sample_merged.csv    # QC-integrated metadata
│       └── logs/                         # Execution logs (optional)
├── results/
│   └── Run2025_01/
│       ├── README.md
│       ├── pipeline.log
│       └── final_report.html
│
├── .gitignore
├── nextflow.config             # Main config file
├── LICENSE
├── README.md
├── CHANGELOG.md
└── environment.yml             # Conda env (if not using containers)
```

## Quick start (Docker edition)

### 1 · Prerequisites

| Tool               | macOS (Homebrew)                                                            | Ubuntu / Debian                                                           | Notes                                                         |
|--------------------|-----------------------------------------------------------------------------|---------------------------------------------------------------------------|---------------------------------------------------------------|
| **Docker ≥ 24**    | `brew install --cask docker`<br/>Launch *Docker Desktop*                     | `sudo apt install docker.io`                                              | Ensure the Docker daemon is running and your user has access. |
| **Nextflow ≥ 23.10** | `brew install nextflow`                                                    | `curl -s https://get.nextflow.io \| bash && sudo mv nextflow /usr/local/bin/` | No extra software needed – the pipeline pulls everything in containers. |

It is also necessary to create locally a folder called databases/ and add the SILVA classifier used by group 2b.

### 2 · Download the Samples

Before running the pipeline (main.nf), you need to download sample data for development. The provided script will create a local folder called `runs/run_id/` following the run naming convention. This folder will already contain 15 paired fastqs in raw_data/ and a metadata.tsv in metadata/. If you don't remove the folder, you will only need to do that once.

To download the samples:

```bash
# make the download script executable
chmod +x create_run.sh

# run the script to download the sample files
./create_run.sh
```

This script is only needed to simulate a real sequencing run during development.
It creates a folder like outputs/run_<run_id>/ and fills it with test FASTQ files in raw_data/ and metadata in metadata/.
You don’t need to run it again unless you're creating a new dev run folder.

```bash
### 3 · Run the pipeline

# run the pipeline using the previously created run folder
nextflow run main.nf --run_id <run_id> -profile docker
# i.e. nextflow run main.nf --run_id R01310525 -profile docker

# resume an interrupted run (skips completed tasks)
nextflow run main.nf --run_id <run_id> -profile docker -resume

## scripts overview
- `create_run.sh` – prepares a run folder with raw fastq files and metadata.
- `download_samples.sh` – used to pull dev/test sample files in a raw_data/ folder.
- `merge_multiqc_metadata.sh` – merges multiqc summary with `metadata_sample.csv`.
```
## Documentation

This repository is complemented by a shared github wiki with module-specific documentation:
https://github.com/egenomics/agb2025/wiki/Group-2-A

---

### Notes

- Filtered fastq's and `outputs/` directories are not pushed to github due to size limits.
- kraken2 is under development and lives in a separate branch (`kraken_integration`).
- multiqc summary merging into metadata is done via csv tools and logged in `outputs/run_<run_id>/`.
