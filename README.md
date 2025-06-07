# agb2025
Repository for the AGB 2025 common class project

[![Documentation  Wiki](https://img.shields.io/static/v1?label=Documentation&message=Wiki&labelColor=black&color=blue&logo=github&logoColor=white)](https://github.com/egenomics/agb2025/wiki)


Wiki you can modify and use to document the process!

## Documentation

In this wiki page you will find the information about the pipeline context, the sample processing and the decisions made through each of the modules (https://github.com/egenomics/agb2025/wiki/Pipeline-Context).

## Project Directory Structure

```text
HdMBioinfo-MicrobiotaPipeline/
├── main.nf                        # Entry point of the pipeline     
├── nextflow.config 
│
├── scripts/                       # Helper scripts and small executables
│   ├── api_csv.py
│   ├── metadata_parsing.R
│   ├── Reorder.R
│   └── process_metadata.py
│
├── modules/                       # DSL2-style modules for each pipeline step
│   └── nf-core/
│       ├── fastqc/
│       ├── trimmomatic/
│       └── kraken2/kraken2
│
├── databases/
│   ├── silva-138-99-nb-classifier.qza
│   └── k2_Human_20230629.tar.gz                
│
├── conf/                          # Configuration files (??)
│
├── controls/
│   ├── metadata_controls.csv
│   ├── taxonomy_control
│   └── ...
│
├── runs/
│   └── R[01-99][DDMMYY] (e.g. R01030525)/
│       ├── raw_data/
│       │   ├── S[A-Z]{2}[1-9][DDMMYY][R1-R2].fastq.gz 
│       │   └── e.g. SAF1030525R1.fastq.gz
│       ├── metadata/
│       │   ├── metadata_cleaned.csv
│       │   ├── metadata_sample.csv
│       │   └── metadata_run.csv
│       ├── qc_reports/
│       │   ├── raw/
│       │   │   ├── *_fastqc.html
│       │   │   └── *_fastqc.zip
│       │   └── trimmed/
│       │       ├── *_trimmed_fastqc.html
│       │       └── *_trimmed_fastqc.zip
│       ├── preprocessed_data/
│       │   ├── *_trimmed.fastq.gz
│       │   └── e.g. SAF1030525R1_trimmed.fastq.gz
│       ├── qiime_output/
│       │   ├── artifacts/
│       │   │   ├── 01_imported_reads
│       │   │   ├── 02_denoised_dada2
│       │   │   ├── 03_summaries
│       │   │   ├── 04_taxonomy
│       │   │   └── 05_phylogeny
│       │   └── relevant_results/
│       └── results/
│           ├── pipeline.log ?
│           ├── multiqc_report.html
│           └── shiny_dashboard_results_app.R
│
├── .gitignore
├── LICENSE
├── README.md
└── CHANGELOG.md
```

## Quick start (Docker edition)

### 1 · Prerequisites

### 1 · Prerequisites

| Tool                   | macOS (Homebrew)                                                            | Ubuntu / Debian                                                           | Notes                                                         |
|------------------------|-----------------------------------------------------------------------------|---------------------------------------------------------------------------|---------------------------------------------------------------|
| **Docker ≥ 24**        | `brew install --cask docker`<br/>Launch *Docker Desktop*                     | `sudo apt install docker.io`                                              | Ensure the Docker daemon is running and your user has access. |
| **Nextflow ≥ 23.10**   | `brew install nextflow`                                                     | `curl -s https://get.nextflow.io \| bash && sudo mv nextflow /usr/local/bin/` | The pipeline pulls everything in containers.                |
| **Memory Requirement** |                                                                             |                                                                           | At least 8 GB of available system memory is required for Kraken2. |

### 2 · Create run and download the samples

Before running the pipeline (main.nf), you need to create the run and download sample data for development. The provided script will create a local folder called `runs/<run_id>/` following the run naming convention. This folder will contain 15 paired fastqs in raw_data/ and a metadata.tsv in metadata/. If you don't remove the folder and you reuse the same run_id, you will only need to do that once.

To create the run and download the samples:

```bash
# make the download script executable
chmod +x create_run.sh
```
```bash
# run the script to download the sample files
./create_run.sh
```

This script is only needed to simulate a real sequencing run during development.

### 3 · Download the Kraken2 db

To run Kraken2, we need to download a database built from the human library. The script will save the databse in the main project folder with the name k2_Human_20230629.tar.gz. If the file already exists, it will not download anything.

```bash
# make the download script executable
chmod +x download_krakendb.sh
```

```bash
# run the script to download the kraken2 database
./download_krakendb.sh
```

### 4 · Download the classifier

First, create the folder databases:
```bash
mkdir databases
```
For running qiime2 and obtaining the feature table, the classifier has to be downloaded, for that you can use the following link:

https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza

Place the downloaded database in the new folder in the main project folder.

### 5 · Run the pipeline
In this step, make sure you have run the following script found in the previous steps:
```bash
./create_run.sh
```
#### Run the pipeline using the previously created run folder
After executing `./create_run.sh`, a `runs/<run_id>/` folder will be created. Copy the `run_id` and use it to run the pipeline.

```bash
nextflow run main.nf --run_id <run_id> -profile docker
# i.e. nextflow run main.nf --run_id R01310525 -profile docker
```

#### Resume an interrupted run (skips completed tasks)
```bash
nextflow run main.nf --run_id <run_id> -profile docker -resume
```

## Scripts Overview
- `create_run.sh` – prepares a run folder with raw fastq files and metadata.
- `download_samples.sh` – used to download fastq files for development, testing or healthy controls. Saves them in raw_data/, test_data/ or healthy_controls/
- `download_krakendb.sh` – Downloads a kraken2 database built from just the human library. Saves it in the main project folder with the name k2_Human_20230629.tar.gz.
- `merge_multiqc_metadata.sh` – merges multiqc summary with `metadata_sample.csv`.

## Documentation

This repository is complemented by a shared github wiki with module-specific documentation:
https://github.com/egenomics/agb2025/wiki/Group-2-A

---

### Notes

- Filtered fastq's and `outputs/` directories are not pushed to github due to size limits.
- kraken2 is under development and lives in a separate branch (`kraken_integration`).
- multiqc summary merging into metadata is done via csv tools and logged in `outputs/run_<run_id>/`.
