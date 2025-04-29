#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """\
         P R O J E C T  R U L E S  Group 2B - QIIME 2 Pipeline
         ===================================
         Input Reads      : ${params.reads}
         Output Directory : ${params.outdir}
         Denoiser         : ${params.denoiser}
         Classifier DB    : ${params.classifier_db}
         Sampling Depth   : ${params.sampling_depth}
         """

// == Define Parameters ==
params.reads = "group2_B/data/simulated_reads/*_R{1,2}.fastq.gz" // Default, override with --reads
params.outdir = "group2_B/results" // Default, override with --outdir
params.denoiser = "dada2" // Options: 'dada2', 'deblur'
params.classifier_db = "path/to/your/qiime2_classifier.qza" // <<< IMPORTANT: Specify path to your classifier (e.g., Silva/Greengenes)
params.sampling_depth = 1000 // <<< IMPORTANT: Adjust based on your data (check feature table summary)
params.trunc_len_f = 100 // DADA2: Forward read truncation length ### dpends on the demux report!!!!!!
params.trunc_len_r = 100 // DADA2: Reverse read truncation length
params.trim_length = 250 // Deblur: Read trim length

// == Input Channel ==
ch_reads = Channel.fromFilePairs(params.reads, flat: true)


// == Workflow ==
workflow {
    // 1. Import Reads
    IMPORT_READS( ch_reads )

    // 2. Denoise
    if (params.denoiser == 'dada2') {
        DENOISE_DADA2( IMPORT_READS.out.demux_qza )
        ch_denoised_table = DENOISE_DADA2.out.table
        ch_denoised_reps = DENOISE_DADA2.out.rep_seqs
    } else if (params.denoiser == 'deblur') {
        DENOISE_DEBLUR( IMPORT_READS.out.demux_qza )
        ch_denoised_table = DENOISE_DEBLUR.out.table
        ch_denoised_reps = DENOISE_DEBLUR.out.rep_seqs
    } else {
        error "Invalid denoiser option: ${params.denoiser}. Choose 'dada2' or 'deblur'."
    }

    // 3. Feature Table Summaries (Useful for picking sampling depth)
    SUMMARIZE_TABLE( ch_denoised_table )
    SUMMARIZE_SEQS( ch_denoised_reps )

    // 4. Taxonomic Classification
    CLASSIFY_TAXONOMY( ch_denoised_reps )

    // 5. Phylogenetic Tree Construction
    BUILD_TREE( ch_denoised_reps )

    // 6. Core Diversity Metrics
    CORE_DIVERSITY( ch_denoised_table, ch_denoised_reps, BUILD_TREE.out.rooted_tree )
}

// == Processes ==

// 1. Import Reads
process IMPORT_READS {
    publishDir "${params.outdir}/qiime2_artifacts/01_imported_reads", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_R1), path(reads_R2)
    
    output:
    path("${sample_id}_demux.qza", emit: demux_qza)
    
    script:
    """
    # Create a manifest file for QIIME2 import
    echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > manifest.tsv
    echo -e "${sample_id}\t\$PWD/${reads_R1}\t\$PWD/${reads_R2}" >> manifest.tsv
    
    qiime tools import \\
        --type 'SampleData[PairedEndSequencesWithQuality]' \\
        --input-path manifest.tsv \\
        --output-path ${sample_id}_demux.qza \\
        --input-format PairedEndFastqManifestPhred33V2
    
    """
    
    stub: // Minimal command for testing without actual execution
    """
    touch ${sample_id}_demux.qza
    """
}



// 2a. Denoise with DADA2
process DENOISE_DADA2 {
    publishDir "${params.outdir}/qiime2_artifacts/02_denoised_dada2", mode: 'copy'

    input:
    path(demux_qza)

    output:
    path("table.qza", emit: table)
    path("rep-seqs.qza", emit: rep_seqs)
    path("denoising-stats.qza", emit: stats)

    script:
    """
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs ${demux_qza} \
        --p-trunc-len-f ${params.trunc_len_f} \
        --p-trunc-len-r ${params.trunc_len_r} \
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats denoising-stats.qza \
        --p-n-threads ${task.cpus}
    """
    stub:
    """
    touch table.qza rep-seqs.qza denoising-stats.qza
    """
}

// 2b. Denoise with Deblur
process DENOISE_DEBLUR {
    publishDir "${params.outdir}/qiime2_artifacts/02_denoised_deblur", mode: 'copy'

    input:
    path(demux_qza)

    output:
    path("table.qza", emit: table)
    path("rep-seqs.qza", emit: rep_seqs)
    path("stats.qza", emit: stats)

    script:
    """
    # Deblur requires quality filtering first
    qiime quality-filter q-score \
     --i-demux ${demux_qza} \
     --o-filtered-sequences demux-filtered.qza \
     --o-filter-stats filter-stats.qza

    # Run Deblur
    qiime deblur denoise-16S \
      --i-demultiplexed-seqs demux-filtered.qza \
      --p-trim-length ${params.trim_length} \
      --o-representative-sequences rep-seqs.qza \
      --o-table table.qza \
      --p-sample-stats \
      --o-stats stats.qza \
      --p-jobs-to-start ${task.cpus}
    """
    stub:
    """
    touch table.qza rep-seqs.qza stats.qza
    """
}

// 3a. Summarize Feature Table
process SUMMARIZE_TABLE {
    publishDir "${params.outdir}/qiime2_visualizations/03_summaries", mode: 'copy'

    input:
    path(table_qza)

    output:
    path("table.qzv", emit: table_summary)

    script:
    """
    qiime feature-table summarize \
      --i-table ${table_qza} \
      --o-visualization table.qzv
      --m-sample-metadata-file path/to/metadata.tsv #define this!!!!
    """
    stub:
    """
    touch table.qzv
    """
}

// 3b. Summarize Representative Sequences
process SUMMARIZE_SEQS {
    publishDir "${params.outdir}/qiime2_visualizations/03_summaries", mode: 'copy'

    input:
    path(rep_seqs_qza)

    output:
    path("rep-seqs.qzv", emit: seqs_summary)

    script:
    """
    qiime feature-table tabulate-seqs \
      --i-data ${rep_seqs_qza} \
      --o-visualization rep-seqs.qzv
    """
    stub:
    """
    touch rep-seqs.qzv
    """
}


// 4. Taxonomic Classification
process CLASSIFY_TAXONOMY {
    publishDir "${params.outdir}/qiime2_artifacts/04_taxonomy", mode: 'copy'

    input:
    path(rep_seqs_qza)

    output:
    path("taxonomy.qza", emit: taxonomy)
    path("taxonomy.qzv", emit: taxonomy_viz) // Add visualization of taxonomy

    script:
    if (params.classifier_db == "path/to/your/qiime2_classifier.qza" || !file(params.classifier_db).exists()) {
        error "Classifier database not found or default path used: ${params.classifier_db}. Please specify a valid path using --classifier_db"
    }
    """
    qiime feature-classifier classify-sklearn \
      --i-classifier ${params.classifier_db} \
      --i-reads ${rep_seqs_qza} \
      --o-classification taxonomy.qza \
      --p-n-jobs ${task.cpus}

    qiime metadata tabulate \
      --m-input-file taxonomy.qza \
      --o-visualization taxonomy.qzv
    """
    stub:
    """
    touch taxonomy.qza taxonomy.qzv
    """
}

// 5. Build Phylogenetic Tree
process BUILD_TREE {
    publishDir "${params.outdir}/qiime2_artifacts/05_phylogeny", mode: 'copy'

    input:
    path(rep_seqs_qza)

    output:
    path("aligned-rep-seqs.qza")
    path("masked-aligned-rep-seqs.qza")
    path("unrooted-tree.qza")
    path("rooted-tree.qza", emit: rooted_tree)

    script:
    """
    qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences ${rep_seqs_qza} \
      --o-alignment aligned-rep-seqs.qza \
      --o-masked-alignment masked-aligned-rep-seqs.qza \
      --o-tree unrooted-tree.qza \
      --o-rooted-tree rooted-tree.qza \
      --p-n-threads ${task.cpus}
    """
    stub:
    """
    touch aligned-rep-seqs.qza masked-aligned-rep-seqs.qza unrooted-tree.qza rooted-tree.qza
    """
}

// 6. Core Diversity Metrics
process CORE_DIVERSITY {
    publishDir "${params.outdir}/qiime2_artifacts/06_diversity_core", mode: 'copy'
    publishDir "${params.outdir}/qiime2_visualizations/06_diversity_core", mode: 'copy', pattern: '*.qzv'

    input:
    path(table_qza)
    path(rep_seqs_qza) // Although not directly used, good practice to show dependency
    path(rooted_tree_qza)
    // Add metadata file input here if available: path(metadata)

    output:
    // Alpha diversity
    path("faith_pd_vector.qza")
    path("observed_features_vector.qza")
    path("shannon_vector.qza")
    path("evenness_vector.qza")
    path("faith-pd-group-significance.qzv")
    path("evenness-group-significance.qzv")
    path("shannon-group-significance.qzv")
    // Beta diversity
    path("unweighted_unifrac_distance_matrix.qza")
    path("weighted_unifrac_distance_matrix.qza")
    path("jaccard_distance_matrix.qza")
    path("bray_curtis_distance_matrix.qza")
    path("unweighted_unifrac_pcoa_results.qza")
    path("weighted_unifrac_pcoa_results.qza")
    path("jaccard_pcoa_results.qza")
    path("bray_curtis_pcoa_results.qza")
    path("unweighted_unifrac_emperor.qzv")
    path("weighted_unifrac_emperor.qzv")
    path("jaccard_emperor.qzv")
    path("bray_curtis_emperor.qzv")

    script:
    if (params.sampling_depth <= 0) {
        error "Sampling depth (--sampling_depth) must be a positive integer. Check the feature table summary (table.qzv) to choose an appropriate value."
    }
    // Metadata is required for group significance testing and PCoA plots with coloring
    // Add --m-metadata-file path/to/metadata.tsv if available
    """
    qiime diversity core-metrics-phylogenetic \
      --i-phylogeny ${rooted_tree_qza} \
      --i-table ${table_qza} \
      --p-sampling-depth ${params.sampling_depth} \
      --output-dir core-metrics-results \
      --p-n-jobs-or-threads ${task.cpus}
      # Add metadata if available: --m-metadata-file ${metadata}

    # Move results out of the directory created by qiime
    mv core-metrics-results/* .
    rmdir core-metrics-results
    """
    stub:
    """
    mkdir core-metrics-results
    touch core-metrics-results/faith_pd_vector.qza core-metrics-results/observed_features_vector.qza core-metrics-results/shannon_vector.qza core-metrics-results/evenness_vector.qza
    touch core-metrics-results/faith-pd-group-significance.qzv core-metrics-results/evenness-group-significance.qzv core-metrics-results/shannon-group-significance.qzv
    touch core-metrics-results/unweighted_unifrac_distance_matrix.qza core-metrics-results/weighted_unifrac_distance_matrix.qza core-metrics-results/jaccard_distance_matrix.qza core-metrics-results/bray_curtis_distance_matrix.qza
    touch core-metrics-results/unweighted_unifrac_pcoa_results.qza core-metrics-results/weighted_unifrac_pcoa_results.qza core-metrics-results/jaccard_pcoa_results.qza core-metrics-results/bray_curtis_pcoa_results.qza
    touch core-metrics-results/unweighted_unifrac_emperor.qzv core-metrics-results/weighted_unifrac_emperor.qzv core-metrics-results/jaccard_emperor.qzv core-metrics-results/bray_curtis_emperor.qzv
    mv core-metrics-results/* .
    rmdir core-metrics-results
    """
}
