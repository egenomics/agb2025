# agb2025

Repository for the AGB 2025 common class project

[![Documentation  Wiki](https://img.shields.io/static/v1?label=Documentation&message=Wiki&labelColor=black&color=blue&logo=github&logoColor=white)](https://github.com/egenomics/agb2025/wiki)

In this wiki page you will find the information about the pipeline context, the sample processing and the decisions made through each of the modules.

![IMG1.png]([https://github.com/egenomics/agb2025/img/IMG1.png](https://github.com/egenomics/agb2025/blob/main/img/IMG1.png))

## Quick start (Docker edition)

### 1 · Prerequisites

| Tool                   | macOS (Homebrew)                                                            | Ubuntu / Debian                                                           | Notes                                                         |
|------------------------|-----------------------------------------------------------------------------|---------------------------------------------------------------------------|---------------------------------------------------------------|
| **Docker ≥ 24**        | `brew install --cask docker`<br/>Launch *Docker Desktop*                     | `sudo apt install docker.io`                                              | Ensure the Docker daemon is running and your user has access. |
| **Nextflow ≥ 23.10**   | `brew install nextflow`                                                     | `curl -s https://get.nextflow.io \| bash && sudo mv nextflow /usr/local/bin/` | The pipeline pulls everything in containers.                |
| **Memory Requirement** |                                                                             |                                                                           | At least 8 GB of available system memory is required for Kraken2. |

## 2 · Setup Installing

This pipeline uses a custom Docker image (agb2025-python) to merge metadata with MultiQC output using pandas and csvkit. To build the Docker Image, run this in the terminal

```bash
docker build -t agb2025-python -f Dockerfile .
```

The image is used in the process that merges metadata.tsv with multiqc_fastqc.txt. This step will fail unless agb2025-python is built locally beforehand.

In addition, run the pipeline, it is mandatory to install two databases: Kraken2 and the classifier SILVA. To install them, run this command line:

```bash
chmod +x INSTALLME.sh
./INSTALLME.sh
```

The `INSTALLME.sh` script will save both databases in databases/.


### 3 · Usage

Use this command line:

```bash
nextflow run main.nf --run_id <run_id> --sampling_depth <number> --auto_rarefaction TRUE -profile docker
# i.e. nextflow run main.nf --run_id R01310525 --sampling_depth 10000 --auto_rarefaction TRUE -profile docker
```

To resume in case of an interrupted run, to skip completed tasks, we highly recommend to use the argument -resume:

```bash
nextflow run main.nf --run_id <run_id> --sampling_depth <number> --auto_rarefaction TRUE -profile docker -resume
```

### 4 ·  Tutorial

To test the developing version of the pipeline (main.nf), you need to create the run and download sample data. Executing the `create._run.sh` script will create a local folder called `runs/<run_id>/` following the run naming convention. This folder will contain 15 paired fastqs in raw_data/ and a metadata.tsv in metadata/. If you don't remove the folder and you reuse the same run_id, you will only need to do that once.

```bash
chmod +x create_run.sh
./create_run.sh
```

This script is only needed to simulate a real sequencing run during development.

After executing `./create_run.sh`, a `runs/<run_id>/` folder will be created. Copy the `run_id` and use it to **run the pipeline**.

```bash
nextflow run main.nf --run_id R01110625 --auto_rarefaction TRUE --sampling_depth 10000 -profile docker -resume
```


## Scripts Overview
- `create_run.sh` – prepares a run folder with raw fastq files and metadata.
- `INSTALLME.sh` – automatically downloads and extracts kraken2 and qiime2 databases.


---

### Notes

- Filtered fastq's and `outputs/` directories are not pushed to github due to size limits.
- multiqc summary merging into metadata is done via csv tools and logged in `outputs/run_<run_id>/`.
