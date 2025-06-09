process EXPORT_TREE {
    label 'qiime2'

    input:
    path(tree_qza)

    output:
    path("phylogenetic_tree.nwk", emit: tree_newick)

    script:
    """
    qiime tools export \
      --input-path ${tree_qza} \
      --output-path exported_tree
    
    mv exported_tree/tree.nwk phylogenetic_tree.nwk
    """
    stub:
    """
    touch phylogenetic_tree.nwk
    """
}

