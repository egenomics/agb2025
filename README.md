# agb2025
Repository for the AGB 2025 common class project

[![Documentation  Wiki](https://img.shields.io/static/v1?label=Documentation&message=Wiki&labelColor=black&color=blue&logo=github&logoColor=white)](https://github.com/egenomics/agb2025/wiki)


Wiki you can modify and use to document the process!

## Documentation

In this wiki page you will find the information about the pipeline context, the sample processing and the decisions made through each of the modules (https://github.com/egenomics/agb2025/wiki/Pipeline-Context).


## Quick start (Docker edition)

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

### 3 · Setup Installing

To run the pipeline (specially kraken2), we have to download a database built from the human library.

On the other hand, for running qiime2 and obtaining the feature table, a classifier has to be downloaded. By default, the classifier downloaded is: https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza. If you want to change the classifier, further information is in the documentation.

For install this requirement, run this command line:

```bash
chmod +x INSTALLME.sh
./INSTALLME.sh
```

The script will save both databases in databases/.

### 5 · Run the pipeline
In this step, make sure you have run the following script found in the previous steps:
```bash
./create_run.sh
```
#### Run the pipeline using the previously created run folder
After executing `./create_run.sh`, a `runs/<run_id>/` folder will be created. Copy the `run_id` and use it to run the pipeline.

```bash
nextflow run main.nf --run_id <run_id> -profile docker --auto_rarefaction TRUE
# i.e. nextflow run main.nf --run_id R01310525 -profile docker
```

#### Resume an interrupted run (skips completed tasks)
```bash
nextflow run main.nf --run_id <run_id> -profile docker -resume
```

## Scripts Overview
- `create_run.sh` – prepares a run folder with raw fastq files and metadata.
- `INSTALLME.sh` – automatically downloads and extracts kraken2 and qiime2 databases.

## Custom Python Image for Metadata Merging
This pipeline uses a custom Docker image (agb2025-python) to merge metadata with MultiQC output using pandas and csvkit.

### Build the Docker Image (once, before running the pipeline)
```bash
docker build -t agb2025-python -f Dockerfile .
```
### Usage
The image is used in the process that merges metadata.tsv with multiqc_fastqc.txt. This step will fail unless agb2025-python is built locally beforehand.


## Documentation

This repository is complemented by a shared github wiki with module-specific documentation:
https://github.com/egenomics/agb2025/wiki/Group-2-A

---

### Notes

- Filtered fastq's and `outputs/` directories are not pushed to github due to size limits.
- multiqc summary merging into metadata is done via csv tools and logged in `outputs/run_<run_id>/`.
