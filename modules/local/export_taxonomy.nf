process EXPORT_TAXONOMY {
    label 'qiime2'
    publishDir "${params.outdir}/qiime_output/relevant_results", mode: 'copy'

    input:
    path(taxonomy_qza)

    output:
    path("taxonomy.tsv", emit: taxonomy_tsv)

    script:
    """
    qiime tools export \
      --input-path ${taxonomy_qza} \
      --output-path exported_taxonomy
    
    mv exported_taxonomy/taxonomy.tsv taxonomy.tsv
    """
    stub:
    """
    touch taxonomy.tsv
    """
}

