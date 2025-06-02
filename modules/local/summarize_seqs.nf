process SUMMARIZE_SEQS {
    label 'qiime2'
    publishDir "${params.outdir}/qiime_output/artifacts/03_summaries", mode: 'copy'

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

