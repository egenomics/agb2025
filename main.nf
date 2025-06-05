#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Could create a subworkflow for fastqc so it's not needed to include the same process twice.
include { FASTQC as FASTQC_RAW  } from './modules/nf-core/fastqc/main.nf'
include { FASTQC as FASTQC_TRIM } from './modules/nf-core/fastqc/main.nf'
include { TRIMMOMATIC           } from './modules/nf-core/trimmomatic/main.nf'
include { KRAKEN2_KRAKEN2 as KRAKEN } from './modules/nf-core/kraken2/kraken2/main.nf'

// QIIME2 modules
include { IMPORT_READS          } from './modules/local/import_reads.nf'
include { DENOISE_DADA2         } from './modules/local/denoise_dada2.nf'
include { DENOISE_DEBLUR        } from './modules/local/denoise_deblur.nf'
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
include { CREATE_RESULTS_SUMMARY } from './modules/local/create_results_summary.nf'


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
    FASTQC_RAW(raw_reads)

    // 1.2. Run Trimmomatic
    TRIMMOMATIC(raw_reads)

    // 1.3. Run FastQC on trimmed reads
    FASTQC_TRIM(
        TRIMMOMATIC.out.trimmed_reads.map { meta, reads ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}_trimmed"
            return tuple(new_meta, reads)
        }
    )

    KRAKEN(
        TRIMMOMATIC.out.trimmed_reads
            .map { meta, reads ->
                tuple(meta, reads)
            },
            file(params.kraken2_db),
            false,
            true
    )

    TRIMMOMATIC.out.trimmed_reads
    .map { it[1] }
    .flatten()
    .collect()
    .set { all_trimmed_files }


    // 2. QIIME
    // 2.1. Import reads
    IMPORT_READS( all_trimmed_files )

    // 2.2. Denoise
    if (params.denoiser == 'dada2') {
        DENOISE_DADA2( IMPORT_READS.out.demux_qza )
        ch_denoised_table = DENOISE_DADA2.out.table
        ch_denoised_reps =  DENOISE_DADA2.out.rep_seqs
    } else if (params.denoiser == 'deblur') {
        DENOISE_DEBLUR( IMPORT_READS.out.demux_qza )
        ch_denoised_table = DENOISE_DEBLUR.out.table
        ch_denoised_reps =DENOISE_DEBLUR.out.rep_seqs
    } else {
        error "Invalid denoiser option: ${params.denoiser}. Choose 'dada2' or 'deblur'."
    }

    // 2.3. Feature Table Summaries
    SUMMARIZE_TABLE( ch_denoised_table, ch_metadata)
    SUMMARIZE_SEQS( ch_denoised_reps )

    // 2.4. Export Feature Table to TSV
    EXPORT_FEATURE_TABLE( ch_denoised_table )

    // 2.5. Export Representative Sequences to FASTA
    EXPORT_REP_SEQS( ch_denoised_reps )

    // 2.6. Taxonomic Classification
    CLASSIFY_TAXONOMY( ch_denoised_reps, ch_classifier )

    // 2.7. Export Taxonomy to TSV
    EXPORT_TAXONOMY( CLASSIFY_TAXONOMY.out.taxonomy )

    // 2.8. Phylogenetic Tree Construction
    BUILD_TREE( ch_denoised_reps )

    // 2.9. Export Phylogenetic Tree to Newick format
    EXPORT_TREE( BUILD_TREE.out.rooted_tree )

    // 2.10. Calculate rarefaction threshold first
    def threshold_channel
    if (params.auto_rarefaction) {
        RAREFACTION_THRESHOLD(ch_denoised_table, ch_metadata)
        threshold_channel = RAREFACTION_THRESHOLD.out.threshold_file
    } else {
        threshold_channel = Channel.value(params.sampling_depth)
    }

    // 2.11. Run alpha diversity only after threshold is calculated
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

    // 2.11. Create Final Results Summary
    CREATE_RESULTS_SUMMARY(
        EXPORT_FEATURE_TABLE.out.feature_table_tsv,
        EXPORT_REP_SEQS.out.rep_seqs_fasta,
        EXPORT_TAXONOMY.out.taxonomy_tsv,
        EXPORT_TREE.out.tree_newick,
        ch_metadata,
        ch_rarefaction_summary.ifEmpty(file("NO_FILE"))
    )

    // println("Workflow completed at: ${new Date()}")
    //println("Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}")
    }

