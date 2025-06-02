process EXPORT_TREE {
    label 'qiime2'
    publishDir "${params.outdir}/qiime_output/relevant_results", mode: 'copy'

    input:
    path(rooted_tree_qza)

    output:
    path("phylogenetic_tree.nwk", emit: tree_newick)

    script:
    """
    qiime tools export \
      --input-path ${rooted_tree_qza} \
      --output-path exported_tree
    
    mv exported_tree/tree.nwk phylogenetic_tree.nwk
    """
    stub:
    """
    touch phylogenetic_tree.nwk
    """
}

