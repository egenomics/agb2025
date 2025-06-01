nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW } from './modules/nf-core/fastqc/main.nf'
include { FASTQC as FASTQC_TRIM } from './modules/nf-core/fastqc/main.nf'
include { TRIMMOMATIC } from './modules/nf-core/trimmomatic/main.nf'


workflow {
    println("Timestamp: ${params.timestamp}")
    println("Output directory: ${params.outdir}")

    // Channel with paired-end fastq files
    Channel.fromFilePairs("runs/${params.run_id}/raw_data/*_{1,2}.fastq.gz", size: 2)
        .ifEmpty { error("No paired FASTQ files found in raw_data/") }
        .map { sample_id, reads ->
            def meta = [
                id: sample_id,
                single_end: false,
            ]
            return tuple(meta, reads)
        }
        .set { raw_reads }

    // Debug: view the channel contents (commented)
    // raw_reads.view { meta, reads -> "Sample: ${meta.id}, Files: ${reads}" }

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

    println("Workflow completed at: ${new Date()}")
    //println("Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}")
    println("Output directory: ${params.outdir}")
}
