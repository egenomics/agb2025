process BUILD_TREE {
    label 'qiime2'
    publishDir "${params.outdir}/qiime_output/artifacts/05_phylogeny", mode: 'copy'

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
