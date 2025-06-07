#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """\
         P R O J E C T  R U L E S  Group 2B - QIIME 2 Pipeline
         ===================================
         Input Reads      : ${params.reads}
         Output Directory : ${params.outdir}
         Denoiser         : ${params.denoiser}
         Classifier DB    : ${params.classifier_db}
         Metadata         : ${params.metadata}
         """

// == Define Parameters ==
params.reads = "group2_B/data/*_{1,2}.fastq.gz" // override with --reads
params.outdir = "group2_B/results" // override with --outdir
params.denoiser = "dada2" // Options: 'dada2', 'deblur'
params.classifier_db = "group2_B/classifier/silva138_noEuk_AB_classifier.qza"
params.metadata = "group2_B/metadata/metadata.tsv"

// params.trunc_len_f = 240 // DADA2: Forward read truncation length ### depends on the demux report!!!!!!
// params.trunc_len_r = 180 // DADA2: Reverse read truncation length



// == Input Channel - Alternative approach ==
// Get the input directory and create channel from there
ch_input_dir = Channel.fromPath("group2_B/data", type: 'dir', checkIfExists: true)

ch_classifier = Channel.fromPath(params.classifier_db, checkIfExists: true)
ch_metadata = Channel.fromPath(params.metadata, checkIfExists: true)

// == Workflow ==
workflow {
    // 1. Import Reads
    IMPORT_READS( ch_input_dir )

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

    // 3. Feature Table Summaries
    SUMMARIZE_TABLE( ch_denoised_table, ch_metadata)
    SUMMARIZE_SEQS( ch_denoised_reps )

    // 4. Taxonomic Classification
    CLASSIFY_TAXONOMY( ch_denoised_reps, ch_classifier )

    // 5. Phylogenetic Tree Construction
    BUILD_TREE( ch_denoised_reps )
}

// == Processes ==

// 1. Import Reads
process IMPORT_READS {
    publishDir "${params.outdir}/qiime2_artifacts/01_imported_reads", mode: 'copy'
    
    input:
    path(input_dir)
    
    output:
    path("demux.qza"), emit: demux_qza
    path("demux.qzv"), emit: demux_qzv
    
    script:
    """
    # Create manifest file by finding all fastq.gz files
    echo -e "sample-id\\tforward-absolute-filepath\\treverse-absolute-filepath" > manifest.tsv
    
    # Find all R1 files and create manifest entries
    for r1_file in ${input_dir}/*_1.fastq.gz; do
        if [ -f "\$r1_file" ]; then
            # Extract sample name (remove path and _1.fastq.gz)
            sample_name=\$(basename "\$r1_file" _1.fastq.gz)
            
            # Construct R2 file path
            r2_file="${input_dir}/\${sample_name}_2.fastq.gz"
            
            # Check if R2 file exists
            if [ -f "\$r2_file" ]; then
                # Add to manifest with absolute paths
                echo -e "\${sample_name}\\t\$(realpath \$r1_file)\\t\$(realpath \$r2_file)" >> manifest.tsv
            else
                echo "Warning: R2 file not found for \$sample_name"
            fi
        fi
    done
    
    # Check if manifest has any data rows
    data_rows=\$(tail -n +2 manifest.tsv | wc -l)
    if [ "\$data_rows" -eq 0 ]; then
        echo "ERROR: No valid sample pairs found in ${input_dir}"
        echo "Looking for files matching pattern: *_1.fastq.gz and *_2.fastq.gz"
        exit 1
    fi
    
    # Import the demultiplexed data
    qiime tools import \\
        --type 'SampleData[PairedEndSequencesWithQuality]' \\
        --input-path manifest.tsv \\
        --output-path demux.qza \\
        --input-format PairedEndFastqManifestPhred33V2
    
    # Create visualization
    qiime demux summarize \\
        --i-data demux.qza \\
        --o-visualization demux.qzv
    """
    
    stub:
    """
    touch demux.qza demux.qzv
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
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats denoising-stats.qza \
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
    path(metadata_file)

    output:
    path("table.qzv", emit: table_summary)

    script:
    """
    qiime feature-table summarize \
      --i-table ${table_qza} \
      --o-visualization table.qzv \
      --m-sample-metadata-file ${metadata_file}

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
    path(classifier_db)

    output:
    path("taxonomy.qza", emit: taxonomy)
    path("taxonomy.qzv", emit: taxonomy_viz) // Add visualization of taxonomy

    script:
    if (params.classifier_db == "group2_B/proba_data/silva138_noEuk_AB_classifier.qza" || !file(params.classifier_db).exists()) {
        error "Classifier database not found or default path used: ${params.classifier_db}. Please specify a valid path using --classifier_db"
    }
    """
    qiime feature-classifier classify-sklearn \
      --i-classifier ${classifier_db} \
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
