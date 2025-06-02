process DENOISE_DEBLUR {
    label 'qiime2'
    publishDir "${params.outdir}/qiime_output/artifacts/02_denoised_deblur", mode: 'copy'

    input:
    path(demux_qza)

    output:
    path("table.qza", emit: table)
    path("rep-seqs.qza", emit: rep_seqs)
    path("stats.qza", emit: stats)

    script:
    """
    # Deblur requires quality filtering first
    qiime quality-filter q-score \
     --i-demux ${demux_qza} \
     --o-filtered-sequences demux-filtered.qza \
     --o-filter-stats filter-stats.qza

    # Run Deblur
    qiime deblur denoise-16S \
      --i-demultiplexed-seqs demux-filtered.qza \
      --p-trim-length ${params.trim_length} \
      --o-representative-sequences rep-seqs.qza \
      --o-table table.qza \
      --p-sample-stats \
      --o-stats stats.qza \
      --p-jobs-to-start ${task.cpus}
    """
    stub:
    """
    touch table.qza rep-seqs.qza stats.qza
    """
}

