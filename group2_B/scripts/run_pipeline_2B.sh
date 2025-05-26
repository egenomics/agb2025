nextflow run group2_B/scripts/qiime2_pipeline.nf \
        --reads "group2_B/data/*_{1,2}.fastq.gz" \
        --outdir "group2_B/results/qiime_output" \
        --classifier_db "group2_B/classifier/silva-138-99-nb-classifier.qza" \
        --sampling_depth 5000 \
        --denoiser dada2 \
        -profile standard
