nextflow.enable.dsl = 2

//
// ── 1. Include all nf‐core modules ─────────────────────────────────────────────
//
include { FASTQC  as FASTQC_RAW  }   from './modules/nf-core/fastqc/main.nf'
include { FASTQC  as FASTQC_TRIM  }   from './modules/nf-core/fastqc/main.nf'
include { TRIMMOMATIC                } from './modules/nf-core/trimmomatic/main.nf'
include { KRAKEN2_KRAKEN2  as KRAKEN  } from './modules/nf-core/kraken2/kraken2/main.nf'
include { MULTIQC                   } from './modules/nf-core/multiqc/main.nf'


//
// ── 2. Create “dummy” files for any MultiQC inputs we don’t actually generate ───
//
file("dummy_multiqc_config.yaml").text   = ""
file("dummy_extra_config.yaml").text     = ""
file("dummy_logo.png").withWriter { out -> /* leave empty */ }
file("dummy_summary.txt").text           = ""
// We’ll also collect real version‐ymls from each module, so no need for a dummy_versions.yml


//
// ── 3. A process to check/download the Kraken2 DB if it’s not already present ──
//
process CHECK_OR_DOWNLOAD_DB {
    tag 'kraken_db'
    container 'docker.io/library/python:3.11-slim'
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
        export HOME=/tmp
        mkdir -p \$HOME/.local/bin
        export PATH=\$HOME/.local/bin:\$PATH
        pip install --quiet --no-cache-dir --user gdown
        gdown --id 1C4aisqMEmUiIv-jNBzxTkKX1AqFKgYen -O kraken_db.tar.gz
        mkdir -p "\$DB_DIR"
        tar -xzf kraken_db.tar.gz -C "\$DB_DIR"
        rm kraken_db.tar.gz
        echo "[INFO] Kraken2 DB download complete."
    else
        echo "[INFO] Kraken2 DB already present – nothing to do."
    fi
    """
}


//
// ── 4. The main workflow ───────────────────────────────────────────────────────
//
workflow {
    //
    // a) Check/download Kraken2 DB
    //
    def kraken2_db = file('refs/kraken2/k2_minusb_20250402')

    println("Timestamp: ${params.timestamp ?: 'N/A'}")
    println("Output directory: ${params.outdir ?: 'N/A'}")

    if ( ! file('refs/kraken2/k2_minusb_20250402').exists() ) {
        println("[INFO] Kraken2 DB not found in refs/kraken2 → downloading")
        CHECK_OR_DOWNLOAD_DB()
    }
    else {
        println("[INFO] Found Kraken2 DB in refs/kraken2 → skipping download")
    }

    //
    // b) Load paired‐end FASTQ files into a channel called `raw_reads`
    //
    Channel
        .fromFilePairs("raw_data/*_{1,2}.fastq.gz", size: 2)
        .ifEmpty { error("No paired FASTQ files found in raw_data/") }
        .map { sample_id, reads ->
            def meta = [ id: sample_id, single_end: false ]
            return tuple(meta, reads)
        }
        .set { raw_reads }

    raw_reads.view { meta, reads -> "Sample: ${meta.id}, Files: ${reads}" }


    //
    // c) RUN FASTQC on raw reads → returns two channels:
    //      1) a channel of (meta, report_file) under “emit: report”
    //      2) a channel of versions.yml under “emit: versions”
    //
    def (ch_fastqc_raw_reports, ch_fastqc_raw_versions) = FASTQC_RAW(raw_reads)



    //
    // d) RUN TRIMMOMATIC on raw reads → returns three channels:
    //      1) trimmed_reads      (meta, [r1_trimmed, r2_trimmed])
    //      2) trim_log           (meta, trim_log.txt)
    //      3) versions.yml       (meta, versions.yml)
    //
    def (ch_trimmed_reads, ch_trimmomatic_logs, ch_trimmomatic_versions) = TRIMMOMATIC(raw_reads)


    //
    // e) RUN FASTQC on the trimmed reads
    //    First, build the “trimmed_for_qc” channel from `ch_trimmed_reads`
    //
    def trimmed_for_qc = ch_trimmed_reads.map { meta, reads ->
        // `reads` is a list [r1_trimmed.fastq.gz, r2_trimmed.fastq.gz]
        def new_meta = meta.clone()
        new_meta.id = "${meta.id}_trimmed"
        return tuple(new_meta, reads)
    }

    // FASTQC_TRIM returns (report, versions)
    def (ch_fastqc_trim_reports, ch_fastqc_trim_versions) = FASTQC_TRIM(trimmed_for_qc)



    //
    // f) RUN Kraken2 on trimmed reads → returns two channels:
    //      1) kraken_report     (meta, kraken_report.tsv)
    //      2) versions.yml      (meta, versions.yml)
    //
    def kraken_input = ch_trimmed_reads.map { meta, reads ->
        tuple(meta, reads)
    }
    def (ch_kraken_reports, ch_kraken_versions) = KRAKEN(
        kraken_input,
        kraken2_db,
        false,
        true
    )


    //
    // g) MERGE all “report/log” files into a single channel for MultiQC
    //
    //    - ch_fastqc_raw_reports   → each emits (meta, fastqc_raw.zip or .html)
    //    - ch_fastqc_trim_reports  → each emits (meta, fastqc_trimmed.zip or .html)
    //    - ch_trimmomatic_logs     → each emits (meta, trimmomatic_log.txt)
    //    - ch_kraken_reports       → each emits (meta, kraken_report.tsv)
    //
    Channel
        .merge(
            ch_fastqc_raw_reports.map   { meta, path -> path },
            ch_fastqc_trim_reports.map  { meta, path -> path },
            ch_trimmomatic_logs.map     { meta, path -> path },
            ch_kraken_reports.map       { meta, path -> path }
        )
        .set { ch_multiqc_files }


    //
    // h) MERGE all “versions.yml” files into one channel (MultiQC can read multiple version files)
    //
    Channel
        .merge(
            ch_fastqc_raw_versions.map    { meta, path -> path },
            ch_trimmomatic_versions.map   { meta, path -> path },
            ch_fastqc_trim_versions.map   { meta, path -> path },
            ch_kraken_versions.map        { meta, path -> path }
        )
        .set { ch_software_versions }


    //
    // i) Create single‐value channels for the other 4 MultiQC inputs:
    //      • multiqc_config.yaml
    //      • extra_multiqc_config.yaml
    //      • logo.png
    //      • workflow_summary.txt
    //
    def ch_multiqc_config    = Channel.value( file("dummy_multiqc_config.yaml") )
    def ch_extra_config      = Channel.value( file("dummy_extra_config.yaml") )
    def ch_multiqc_logo      = Channel.value( file("dummy_logo.png") )
    def ch_workflow_summary  = Channel.value( file("dummy_summary.txt") )


    //
    // j) CALL MultiQC with exactly SIX inputs, in the same order as its `input:` block:
    //      1) ch_multiqc_files      (all report/log files)
    //      2) ch_multiqc_config     (base config)
    //      3) ch_extra_config       (extra overrides)
    //      4) ch_multiqc_logo       (logo)
    //      5) ch_software_versions  (one or more versions.yml)
    //      6) ch_workflow_summary   (summary.txt)
    //
    MULTIQC(
        ch_multiqc_files,
        ch_multiqc_config,
        ch_extra_config,
        ch_multiqc_logo,
        ch_software_versions,
        ch_workflow_summary
    )


    //
    // k) Print completion status
    //
    println("Workflow completed at: ${new Date()}")
    println("Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}")
    println("Output directory: ${params.outdir ?: 'N/A'}")
}
