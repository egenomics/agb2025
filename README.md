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
├── output
│   ├── run1
│   │    ├── taxonomy
│   │    ├── phylogeny
│   │    ├── imported_reads
│   │    ├── denoised_dada2
│   │    └── ...
│   ├── run2
│   └── run3
│
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

### 2 · Download the Samples

Before running the pipeline (main.nf), you need to download sample data for development. The provided script will create a local folder called `raw_data` and download 15 sample datasets into it.

To download the samples:

```bash
# make the download script executable
chmod +x download_samples.sh

# run the script to download the sample files
./download_samples.sh
```

After running this, the sample FASTQ files will be available directly inside the `raw_data` folder. This folder will not pushed to GitHub.

### 3 · Run the pipeline
```bash
# minimal
nextflow run main.nf -profile docker

# resume an interrupted run (skips completed tasks)
nextflow run main.nf -profile docker -resume
```
