process CREATE_RESULTS_SUMMARY {
    label 'qiime2'
    publishDir "${params.outdir}/qiime_output/relevant_results", mode: 'copy'

    input:
    path(feature_table_tsv)
    path(rep_seqs_fasta)
    path(taxonomy_tsv)
    path(tree_newick)
    path(metadata_file)
    path(rarefaction_summary)

    output:
    path("analysis_summary.txt")
    path("feature_table_with_taxonomy.tsv")

    script:
    """
    # Create analysis summary
    echo "QIIME 2 Analysis Summary" > analysis_summary.txt
    echo "========================" >> analysis_summary.txt
    echo "Date: \$(date)" >> analysis_summary.txt
    echo "Denoiser used: ${params.denoiser}" >> analysis_summary.txt
    echo "Auto rarefaction: ${params.auto_rarefaction}" >> analysis_summary.txt
    if [ "${params.auto_rarefaction}" = "true" ]; then
        echo "Rarefaction threshold: ${params.auto_rarefaction}" >> analysis_summary.txt
    else
        echo "Manual sampling depth: ${params.sampling_depth}" >> analysis_summary.txt
    fi    
    echo "" >> analysis_summary.txt
    
    # Count samples and features
    n_samples=\$(head -1 ${feature_table_tsv} | awk -F'\t' '{print NF-1}')
    n_features=\$(tail -n +2 ${feature_table_tsv} | wc -l)
    n_seqs=\$(grep -c "^>" ${rep_seqs_fasta})
    
    echo "Number of samples: \$n_samples" >> analysis_summary.txt
    echo "Number of features (ASVs/OTUs): \$n_features" >> analysis_summary.txt
    echo "Number of representative sequences: \$n_seqs" >> analysis_summary.txt
    echo "" >> analysis_summary.txt

    # Add rarefaction summary if available
    if [ "${rarefaction_summary}" != "NO_FILE" ] && [ -f "${rarefaction_summary}" ]; then
        echo "Rarefaction Analysis:" >> analysis_summary.txt
        echo "--------------------" >> analysis_summary.txt
        cat ${rarefaction_summary} >> analysis_summary.txt
        echo "" >> analysis_summary.txt
    fi
    
    echo "Output files:" >> analysis_summary.txt
    echo "- feature_table.tsv: Feature abundance table" >> analysis_summary.txt
    echo "- representative_sequences.fasta: Representative sequences" >> analysis_summary.txt
    echo "- taxonomy.tsv: Taxonomic classifications" >> analysis_summary.txt
    echo "- phylogenetic_tree.nwk: Phylogenetic tree in Newick format" >> analysis_summary.txt
    echo "- feature_table_with_taxonomy.tsv: Combined feature table with taxonomy" >> analysis_summary.txt
    echo "" >> analysis_summary.txt

    # Create combined feature table with taxonomy using Python
    python3 << 'EOF'
import sys
import pandas as pd

# Load feature table
feature_table = pd.read_csv('${feature_table_tsv}', sep='\\t', index_col=0, skiprows=1)

# Load taxonomy
taxonomy = pd.read_csv('${taxonomy_tsv}', sep='\\t', index_col=0)

# Merge tables
combined = feature_table.join(taxonomy, how='left')

# Save combined table
combined.to_csv('feature_table_with_taxonomy.tsv', sep='\\t')
print("Combined feature table with taxonomy created successfully!")
EOF
    """
    stub:
    """
    touch analysis_summary.txt feature_table_with_taxonomy.tsv
    """
}
