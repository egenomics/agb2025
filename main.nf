#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Could create a subworkflow for fastqc so it's not needed to include the same process twice.

// Preprocessing and QC
include { FASTQC as FASTQC_RAW  }  from './modules/nf-core/fastqc/main.nf'
include { FASTQC as FASTQC_TRIM }  from './modules/nf-core/fastqc/main.nf'
include { TRIMMOMATIC }             from './modules/nf-core/trimmomatic/main.nf'

// Classification
include { KRAKEN2_KRAKEN2 as KRAKEN } from './modules/nf-core/kraken2/kraken2/main.nf'
include { MULTIQC               } from './modules/nf-core/multiqc/main.nf'

// Metadata operations
include { MERGE_METADATA_MULTIQC }       from './modules/local/merge_metadata.nf'
include { HOMOSAPINENS_CONTAMINATION }   from './modules/local/homo_contam.nf'
include { CLASSIFY_QUALITY_SAMPLES }     from './modules/local/classify_quality.nf'


// QIIME2 modules
include { IMPORT_READS          } from './modules/local/import_reads.nf'
include { DENOISE_DADA2         } from './modules/local/denoise_dada2.nf'
include { CLASSIFY_TAXONOMY     } from './modules/local/classify_taxonomy.nf'
include { BUILD_TREE            } from './modules/local/build_tree.nf'
include { SUMMARIZE_SEQS        } from './modules/local/summarize_seqs.nf'
include { SUMMARIZE_TABLE       } from './modules/local/summarize_table.nf'
include { EXPORT_FEATURE_TABLE  } from './modules/local/export_feature_table.nf'
include { EXPORT_TAXONOMY       } from './modules/local/export_taxonomy.nf'
include { EXPORT_TREE           } from './modules/local/export_tree.nf'
include { EXPORT_REP_SEQS       } from './modules/local/export_rep_seqs.nf'
include { RAREFACTION_THRESHOLD } from './modules/local/rarefaction_threshold.nf'
include { ALPHA_DIVERSITY } from './modules/local/alpha_diversity.nf'
include { EXPORT_ALPHAPLOT } from './modules/local/export_alphaplot.nf'
include { CREATE_RESULTS_SUMMARY } from './modules/local/create_results_summary.nf'
include { SUGGEST_TRUNCATION_LENGTHS } from './modules/local/suggest_truncation_lengths.nf'
include { SUMMARIZE_DEMUX } from './modules/local/summarize_demux.nf'


workflow {
    println("Output directory: ${params.outdir}")

    // Check some files exist
    if ( !file(params.metadata).exists() ) {
        error "Metadata file not found: ${params.metadata}"
    }
    if ( !file(params.classifier_db).exists() ) {
        error "Classifier database not found: ${params.classifier_db}"
    }
    if ( !file(params.kraken2_db).exists() ) {
        error "Kraken2 database not found: ${params.kraken2_db}"
    }

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
        // Channel with metadata
        ch_metadata = Channel.fromPath(params.metadata, checkIfExists: true)
        // Channel with classifier
        ch_classifier = Channel.fromPath(params.classifier_db, checkIfExists: true)


    // 1. Pre-processing
    // 1.1. Run FastQC on raw reads
    def fastqc_raw_outputs = FASTQC_RAW(raw_reads)
    def ch_fastqc_raw_reports = fastqc_raw_outputs[0]
    def ch_fastqc_raw_versions = fastqc_raw_outputs[1]

    // 1.2. Run Trimmomatic
    def (ch_trimmed_reads, ch_trimmomatic_logs, ch_trimmomatic_versions) = TRIMMOMATIC(raw_reads)

    // 1.3. Run FastQC on trimmed reads
    def trimmed_for_qc = ch_trimmed_reads.map { meta, reads ->
        def new_meta = meta.clone()
        new_meta.id = "${meta.id}_trimmed"
        return tuple(new_meta, reads)
    }
    def (ch_fastqc_trim_reports, ch_fastqc_trim_versions) = FASTQC_TRIM(trimmed_for_qc)


    // 1.4. Run Kraken2
    def kraken_input = ch_trimmed_reads.map { meta, reads ->
        def fixedReads = reads.collect { file ("runs/${params.run_id}/trimmed_reads/${it.getName()}") }
        tuple(meta, fixedReads)
    }

    def (ch_kraken_reports, ch_kraken_versions) = KRAKEN(
        kraken_input,
        file(params.kraken2_db),
        false,
        true
    )

    ch_multiqc_trigger = Channel
        .empty()
        .mix(ch_fastqc_raw_reports, ch_fastqc_trim_reports, ch_trimmomatic_logs, ch_kraken_reports)
        .collect()

    // 1.5. Run MultiQC
    def (multiqc_output) = MULTIQC(ch_multiqc_trigger.map { file("runs/${params.run_id}/") })

    TRIMMOMATIC.out.trimmed_reads
    .map { it[1] }
    .flatten()
    .collect()
    .set { all_trimmed_files }

    // Metadata handling
    ch_metadata   = Channel.fromPath(params.metadata, checkIfExists: true)
    ch_classifier = Channel.fromPath(params.classifier_db, checkIfExists: true)

    def metadata_tsv_path      = "runs/${params.run_id}/metadata/metadata.tsv"
    def multiqc_path           = "runs/${params.run_id}/multiqc/multiqc_data/multiqc_fastqc.txt"

    metadata_sample_ch = Channel.fromPath(metadata_tsv_path, checkIfExists: true)
    kraken_reports_ch = KRAKEN.out.report.map { meta, report -> report }

    // Merge Metadata and MultiQC (if available)
    merged_csv_ch = MERGE_METADATA_MULTIQC(
        metadata_sample_ch,
        multiqc_output.map { file(multiqc_path) }
    )

    // Augment with Human Contamination Data
    human_augmented_ch = HOMOSAPINENS_CONTAMINATION(
        merged_csv_ch,
        kraken_reports_ch
    )
    
    // Perform Quality Classification on Augmented Data
    CLASSIFY_QUALITY_SAMPLES(human_augmented_ch)
    
    // 2. QIIME
    // 2.1. Import reads
    IMPORT_READS( all_trimmed_files )

    // 2.2. Summarize Demux artifact
    SUMMARIZE_DEMUX( IMPORT_READS.out.demux_qza )

    // 2.3. Generate Truncation suggestion report
    SUGGEST_TRUNCATION_LENGTHS( SUMMARIZE_DEMUX.out.demux_qzv )

    // 2.4. Denoise
    // Let the user know if they are using default truncation values
    if (params.trunc_len_f == 0 || params.trunc_len_r == 0) {
        SUGGEST_TRUNCATION_LENGTHS.out.view()
        error "DADA2 requires truncation lengths. No values for --trunc_len_f or --trunc_len_r were provided. A suggestion report has been generated in the results directory. Please inspect it and re-run with appropriate values."
    }
    DENOISE_DADA2( IMPORT_READS.out.demux_qza )
    ch_denoised_table = DENOISE_DADA2.out.table
    ch_denoised_reps =  DENOISE_DADA2.out.rep_seqs

    // 2.5. Feature Table Summaries
    SUMMARIZE_TABLE( ch_denoised_table, ch_metadata)
    SUMMARIZE_SEQS( ch_denoised_reps )

    // 2.6. Export Feature Table to TSV
    EXPORT_FEATURE_TABLE( ch_denoised_table )

    // 2.7. Export Representative Sequences to FASTA
    EXPORT_REP_SEQS( ch_denoised_reps )

    // 2.8. Taxonomic Classification
    CLASSIFY_TAXONOMY( ch_denoised_reps, ch_classifier )

    // 2.9. Export Taxonomy to TSV
    EXPORT_TAXONOMY( CLASSIFY_TAXONOMY.out.taxonomy )

    // 2.10. Phylogenetic Tree Construction
    BUILD_TREE( ch_denoised_reps )

    EXPORT_TREE( BUILD_TREE.out.rooted_tree )

    // 2.11. Calculate rarefaction threshold first
    def threshold_channel
    if (params.auto_rarefaction) {
        RAREFACTION_THRESHOLD(ch_denoised_table, ch_metadata)
        threshold_channel = RAREFACTION_THRESHOLD.out.threshold_file
    } else {
        threshold_channel = Channel.value(params.sampling_depth)
    }

    // 2.12. Run alpha diversity only after threshold is calculated
    ALPHA_DIVERSITY(
        ch_denoised_table,
        BUILD_TREE.out.rooted_tree,
        ch_metadata,
        threshold_channel
    )

    if (params.auto_rarefaction) {
        RAREFACTION_THRESHOLD.out.summary.view {
            "\n=== RAREFACTION THRESHOLD CALCULATED ===\n${it.text}\n========================================\n"
        }
        ch_rarefaction_summary = RAREFACTION_THRESHOLD.out.summary
    } else {
        ch_rarefaction_summary = Channel.empty()
    }

    // 2.13. Export alpha visualization to a plot
    EXPORT_ALPHAPLOT(
        ALPHA_DIVERSITY.out[1]  // This should be the .qzv file from alpha diversity
    )

    // 2.14. Create Final Results Summary
    CREATE_RESULTS_SUMMARY(
        EXPORT_FEATURE_TABLE.out.feature_table_tsv,
        EXPORT_REP_SEQS.out.rep_seqs_fasta,
        EXPORT_TAXONOMY.out.taxonomy_tsv,
        EXPORT_TREE.out.tree_newick,
        ch_metadata,
        ch_rarefaction_summary.ifEmpty(file("NO_FILE"))
    )

}