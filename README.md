# Welcome to **AGB 2025 Microbiota Classification Pipeline**

### Repository for the AGB 2025 common class project.
AGB 2025 is a **dockerised Nextflow workflow** that turns raw stool sequencing reads into clinically actionable labels (*healthy* vs *non-healthy*) and rich microbiome analytics. Built by the students of the subject AGB_2025.

---

[![Documentation  Wiki](https://img.shields.io/static/v1?label=Documentation&message=Wiki&labelColor=black&color=blue&logo=github&logoColor=white)](https://github.com/egenomics/agb2025/wiki)

In this wiki page you will find the information about the pipeline context, the sample processing and the decisions made through each of the modules.

![IMG1.png](https://github.com/egenomics/agb2025/blob/main/img/IMG1.png)
## Quick start (Docker edition)

### Prerequisites to start the pipeline

| Tool                   | macOS (Homebrew)                                                            | Ubuntu / Debian                                                           | Notes                                                         |
|------------------------|-----------------------------------------------------------------------------|---------------------------------------------------------------------------|---------------------------------------------------------------|
| **Docker ≥ 24**        | `brew install --cask docker`<br/>Launch *Docker Desktop*                     | `sudo apt install docker.io`                                              | Ensure the Docker daemon is running and your user has access. |
| **Nextflow ≥ 23.10**   | `brew install nextflow`                                                     | `curl -s https://get.nextflow.io \| bash && sudo`<br>`chmod +x nextflow`<br>` mv nextflow /usr/local/bin/` | The pipeline pulls everything in containers.                |
| **Memory Requirement** |                                                                             |                                                                           | At least 4 GB of available system memory is required for Kraken2. |

## 1 · Setup Installing

Before running the pipeline you must **(i)** clone the repository, **(ii)** make the custom Docker image available, **(iii)** download the reference databases and **(iv)** create the run.

### 1.1 Clone the repository
```bash
git clone https://github.com/egenomics/agb2025.git
cd agb2025
```
### 1.2 Docker image
This pipeline uses a **custom Docker image (agb2025-python)** to merge metadata with MultiQC output using pandas and csvkit. **To build the Docker Image, run this in the terminal**

```bash
docker build -t agb2025-python -f Dockerfile .
```

The image is used in the process that merges metadata.tsv with multiqc_fastqc.txt. This step will fail unless agb2025-python is built locally beforehand.

### 1.3 Download the reference databases
In addition, **it is mandatory to install two databases: Kraken2 and the classifier SILVA**. To install them, run this command line:

```bash
chmod +x INSTALLME.sh
./INSTALLME.sh
```

The `INSTALLME.sh` script will save both databases in databases/.

### 1.4 Create the run

**To test the developing version of the pipeline (main.nf), a run needs to be created, sample raw data downloaded and the corresponding metadata created**. Executing the `create._run.sh` script will create a local folder called `runs/<run_id>/` following the run naming convention. This folder will contain 15 paired fastqs in raw_data/ and a metadata.tsv in metadata/. If you don't remove the folder and you reuse the same run_id, you will only need to do that once.

```bash
chmod +x create_run.sh
./create_run.sh
```

*This script is only needed to simulate a real sequencing run during development.*

After executing `./create_run.sh`, a `runs/<run_id>/` folder will be created. Copy the `run_id`. It will be used to **run the pipeline**. You can alternatively copy the run_id from the last line of the output of the `./create_run.sh` script.

## 2 · Running the pipeline

After copying the `run_id`, use this command line to **run the pipeline**:

```bash
nextflow run main.nf --run_id <run_id> --sampling_depth <number> --auto_rarefaction TRUE -profile docker
# i.e. nextflow run main.nf --run_id R01310525 --sampling_depth 10000 --auto_rarefaction TRUE -profile docker
```

**To resume in case of an interrupted run, to skip completed tasks, we highly recommend to use the argument -resume:**

```bash
nextflow run main.nf --run_id <run_id> --sampling_depth <number> --auto_rarefaction TRUE -profile docker -resume
```
---

## 3 · Data Visualization

To explore the pipeline outputs interactively, you can launch the integrated Shiny app and visualize taxonomy profiles and sample metadata.

### 3.1 Launching the Shiny App

To start the app, give execution permissions and run the launch script:

```bash
chmod +x shiny_app.sh
./shiny_app.sh
```

Once the app is running, a browser window will open where you can manually upload your processed files for visualization.

### 3.2 Using Pipeline Output

To directly visualize the output generated by the pipeline, you need to convert QIIME-formatted files to the long format expected by the app. To do it, you can use the following script:

```bash
python visualization/convert_qiime_to_long.py [taxonomy.tsv path] [feature_table.tsv path] [output_file path]
```

This will create a .tsv file suitable for upload into the Shiny app.

### Scripts Overview
- `create_run.sh` – prepares a run folder with raw fastq files and metadata.
- `INSTALLME.sh` – automatically downloads and extracts kraken2 and qiime2 databases.

---

### Notes

- Filtered fastq's and `outputs/` directories are not pushed to github due to size limits.
- multiqc summary merging into metadata is done via csv tools and logged in `outputs/run_<run_id>/`.
