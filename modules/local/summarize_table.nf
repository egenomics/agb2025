process SUMMARIZE_TABLE {
    label 'qiime2'

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
