process CLASSIFY_TAXONOMY {
    label 'qiime2'
    publishDir "${params.outdir}/artifacts/04_taxonomy", mode: 'copy'

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
