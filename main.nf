nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW } from './modules/nf-core/fastqc/main.nf'
include { FASTQC as FASTQC_TRIM } from './modules/nf-core/fastqc/main.nf'
include { TRIMMOMATIC } from './modules/nf-core/trimmomatic/main.nf'
include { KRAKEN2_KRAKEN2 as KRAKEN } from './modules/nf-core/kraken2/kraken2/main.nf'


process CHECK_OR_DOWNLOAD_DB {

    tag 'kraken_db'

    container 'docker.io/library/bash:latest' // no Python needed

    publishDir "refs/kraken2", mode: 'copy'

    cache 'lenient'

    output:
    path "k2_minusb_20250402", emit: db_dir

    script:
    """
        set -euo pipefail

        if [ ! -d "k2_minusb_20250402" ]; then
            echo "[INFO] Kraken2 DB not found – downloading …"

            curl -L 'https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20250402.tar.gz' -o kraken_db.tar.gz

            mkdir -p "k2_minusb_20250402"
            tar -xzf kraken_db.tar.gz -C "k2_minusb_20250402"
            rm kraken_db.tar.gz

            echo "[INFO] Kraken2 DB download complete."
        else
            echo "[INFO] Kraken2 DB already present – Using this DB."
        fi
    """
}


workflow {
    println("Timestamp: ${params.timestamp}")
    println("Output directory: ${params.outdir}")

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

    // ── only download if we don’t already have the DB on the host ──
    println("[INFO] Checking Kraken2 DB in refs/kraken2...")
    CHECK_OR_DOWNLOAD_DB()

    // Run Kraken2
    KRAKEN(
        TRIMMOMATIC.out.trimmed_reads.map { meta, reads ->
            tuple(meta, reads)
        },
        CHECK_OR_DOWNLOAD_DB.out.db_dir.map { db_dir ->
            return db_dir
        },
        false,
        true,
    )

    println("Workflow completed at: ${new Date()}")
    //println("Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}")
    println("Output directory: ${params.outdir}")
}
