process EXPORT_REP_SEQS {
    label 'qiime2'

    input:
    path(rep_seqs_qza)

    output:
    path("representative_sequences.fasta", emit: rep_seqs_fasta)

    script:
    """
    qiime tools export \
      --input-path ${rep_seqs_qza} \
      --output-path exported_seqs
    
    mv exported_seqs/dna-sequences.fasta representative_sequences.fasta
    """
    stub:
    """
    touch representative_sequences.fasta
    """
}
