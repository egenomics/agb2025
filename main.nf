nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW } from './modules/nf-core/fastqc/main.nf'
include { FASTQC as FASTQC_TRIM } from './modules/nf-core/fastqc/main.nf'
include { TRIMMOMATIC } from './modules/nf-core/trimmomatic/main.nf'
include { KRAKEN2_KRAKEN2 as KRAKEN } from './modules/nf-core/kraken2/kraken2/main.nf'


process CHECK_OR_DOWNLOAD_DB {

    tag 'kraken_db'

    // a very small image with Python – lets us `pip install gdown`
    container 'docker.io/library/python:3.11-slim'

    /* where the files will be kept on the host                *
     * (the same place the rest of your pipeline already uses) */
    publishDir "refs/kraken2", mode: 'copy'

    cache 'lenient'

    output:
    path "k2_minusb_20250402", emit: db_dir

    script:
    """
        set -euo pipefail

        DB_DIR='k2_minusb_20250402'

        if [ ! -d "\$DB_DIR" ]; then
            echo "[INFO] Kraken2 DB not found – downloading …"

            # give the unprivileged user a writable HOME and PATH
            export HOME=/tmp
            mkdir -p \$HOME/.local/bin
            export PATH=\$HOME/.local/bin:\$PATH

            # install gdown only for this user
            pip install --quiet --no-cache-dir --user gdown

            # fetch the archive
            gdown --id 1C4aisqMEmUiIv-jNBzxTkKX1AqFKgYen -O kraken_db.tar.gz

            # create the DB folder and unpack everything into it
            mkdir -p "\$DB_DIR"
            tar -xzf kraken_db.tar.gz -C "\$DB_DIR"

            rm kraken_db.tar.gz
            echo "[INFO] Kraken2 DB download complete."
        else
            echo "[INFO] Kraken2 DB already present – nothing to do."
        fi
    """
}

workflow {
    def kraken2_db = file('refs/kraken2/k2_minusb_20250402')

    println("Timestamp: ${params.timestamp}")
    println("Output directory: ${params.outdir}")

    // ── only download if we don’t already have the DB on the host ──
    if (!file('refs/kraken2/k2_minusb_20250402').exists()) {
        println("[INFO] Kraken2 DB not found in refs/kraken2 → downloading")
        CHECK_OR_DOWNLOAD_DB()
    }
    else {
        println("[INFO] Found Kraken2 DB in refs/kraken2 → skipping download")
    }

    // Channel with paired-end fastq files
    Channel.fromFilePairs("raw_data/*_{1,2}.fastq.gz", size: 2)
        .ifEmpty { error("No paired FASTQ files found in raw_data/") }
        .map { sample_id, reads ->
            def meta = [
                id: sample_id,
                single_end: false,
            ]
            return tuple(meta, reads)
        }
        .set { raw_reads }

    // Debug: view the channel contents
    raw_reads.view { meta, reads -> "Sample: ${meta.id}, Files: ${reads}" }

    // Run FastQC on raw reads
    FASTQC_RAW(raw_reads)

    // Run Trimmomatic
    TRIMMOMATIC(raw_reads)

    // Run FastQC on trimmed reads
    FASTQC_TRIM(
        TRIMMOMATIC.out.trimmed_reads.map { meta, reads ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}_trimmed"
            return tuple(new_meta, reads)
        }
    )

    // Run Kraken2
    KRAKEN(
        TRIMMOMATIC.out.trimmed_reads.map { meta, reads ->
            tuple(meta, reads)
        },
        kraken2_db,
        false,
        true,
    )

    println("Workflow completed at: ${new Date()}")
    println("Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}")
    println("Output directory: ${params.outdir}")
}
