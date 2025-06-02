process DENOISE_DADA2 {
    label 'qiime2'
    publishDir "${params.outdir}/artifacts/02_denoised_dada2", mode: 'copy'

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
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
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

