    nextflow run group2_B/scripts/qiime2_pipeline.nf \
        --reads "group2_B/data/simulated_reads/*_R{1,2}.fastq.gz" \
        --outdir "group2_B/results/qiime_output" \
        --classifier_db "/path/to/your/downloaded/silva-138-99-nb-classifier.qza" \
        --sampling_depth 5000 \
        -profile standard # Or specify a profile for your execution environment (e.g., docker, singularity)