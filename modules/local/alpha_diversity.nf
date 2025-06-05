process ALPHA_DIVERSITY {
    label 'qiime2'
    publishDir "${params.outdir}/artifacts/06_alpha_diversity", mode: 'copy'
    
    input:
    path(table_qza)
    path(rooted_tree_qza)
    path(metadata_file)
    path(threshold_file)
    
    output:
    path("*.qza")
    path("*.qzv")
    
    script:
    """
    # Debugging: Print threshold file content
    echo "Threshold file content:"
    cat ${threshold_file}
    
    # Read the threshold from the file
    THRESHOLD=\$(cat ${threshold_file})
    
    echo "Using rarefaction threshold: \$THRESHOLD"
    
    # Calculate alpha diversity metrics
    qiime diversity alpha-rarefaction \\
      --i-table ${table_qza} \\
      --i-phylogeny ${rooted_tree_qza} \\
      --p-max-depth \$THRESHOLD \\
      --m-metadata-file ${metadata_file} \\
      --o-visualization alpha-rarefaction.qzv
    
    # Core alpha diversity metrics at the selected threshold
    qiime diversity core-metrics-phylogenetic \\
      --i-table ${table_qza} \\
      --i-phylogeny ${rooted_tree_qza} \\
      --p-sampling-depth \$THRESHOLD \\
      --m-metadata-file ${metadata_file} \\
      --output-dir core-metrics-results
    
    # Move results to current directory
    mv core-metrics-results/* .
    """
    
    stub:
    """
    echo "5000" > threshold_file
    touch alpha-rarefaction.qzv
    touch observed_features_vector.qza shannon_vector.qza evenness_vector.qza
    touch faith_pd_vector.qza rarefied_table.qza
    """
}