process IMPORT_READS {
    label 'qiime2'

    input:
    collect(path(trimmed_reads))
    
    output:
    path("demux.qza", emit: demux_qza)
    path("demux.qzv", emit: demux_qzv)
    
    script:
    """
    echo -e "sample-id\\tforward-absolute-filepath\\treverse-absolute-filepath" > manifest.tsv
    
    # Procesar TODOS los archivos trimmeados
    for r1_file in *.paired.trim_1.fastq.gz; do
        if [ -f "\$r1_file" ]; then
            sample_name=\$(basename "\$r1_file" .paired.trim_1.fastq.gz)
            r2_file="\${sample_name}.paired.trim_2.fastq.gz"
            
            if [ -f "\$r2_file" ]; then
                echo -e "\${sample_name}\\t\$(realpath \$r1_file)\\t\$(realpath \$r2_file)" >> manifest.tsv
            fi
        fi
    done
    
    qiime tools import \\
        --type 'SampleData[PairedEndSequencesWithQuality]' \\
        --input-path manifest.tsv \\
        --output-path demux.qza \\
        --input-format PairedEndFastqManifestPhred33V2
    
    qiime demux summarize \\
        --i-data demux.qza \\
        --o-visualization demux.qzv
    """
}
