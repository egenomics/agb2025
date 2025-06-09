process DENOISE_DADA2 {
    label 'qiime2'

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
        --p-trim-left-f ${params.trim_left_f} \
        --p-trim-left-r ${params.trim_left_r} \
        --p-trunc-len-f ${params.trunc_len_f} \
        --p-trunc-len-r ${params.trunc_len_r} \
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats denoising-stats.qza \
    """
    stub:
    """
    touch table.qza rep-seqs.qza denoising-stats.qza
    """
}

