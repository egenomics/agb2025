process SUMMARIZE_DEMUX {
    label 'qiime2'

    input:
    path(demux_qza)

    output:
    path("demux.qzv"), emit: demux_qzv

    script:
    """
    qiime demux summarize \\
        --i-data ${demux_qza} \\
        --o-visualization demux.qzv
    """
    stub:
    """
    touch demux.qzv
    """
} 