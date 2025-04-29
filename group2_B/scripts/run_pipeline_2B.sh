    nextflow run group2_B/scripts/qiime2_pipeline.nf \
        --reads "group2_B/data/simulated_fastqs/*_{R1,R2}.fastq.gz"\
        --outdir "group2_B/results/qiime_output" \
        --classifier_db "group2_B/databases/Greengenes_16S_2011_1.arb.gz" \
        --sampling_depth 5000 \
        --denoiser dada2 \
        -profile standard # Or specify a profile for your execution environment (e.g., docker, singularity)
