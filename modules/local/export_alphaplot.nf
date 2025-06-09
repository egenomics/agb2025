process EXPORT_ALPHAPLOT {
    label 'qiime2'

    input:
    path(alpha_rarefaction_qzv)

    output:
    path("alpha_rarefaction/*")
    
    script:
    """
    qiime tools export \
      --input-path ${alpha_rarefaction_qzv} \
      --output-path alpha_rarefaction
    """

    stub:
    """
    mkdir -p alpha_rarefaction
    touch alpha_rarefaction/index.html
    touch alpha_rarefaction/data.tsv
    """
}