process EXPORT_FEATURE_TABLE {
    label 'qiime2'

    input:
    path(table_qza)

    output:
    path("feature_table.tsv", emit: feature_table_tsv)

    script:
    """
    qiime tools export \
      --input-path ${table_qza} \
      --output-path exported_table
    
    # Convert BIOM to TSV
    biom convert \
      -i exported_table/feature-table.biom \
      -o feature_table.tsv \
      --to-tsv
    """
    stub:
    """
    touch feature_table.tsv
    """
}
