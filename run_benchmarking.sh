#!/bin/bash

echo " MICROBIOME PIPELINE BENCHMARKING v15"
echo " Genus-Level Validation System"
echo "============================================"
echo ""

# Ensure we're in the right directory
cd /home/annie/Documents/Bioinfo/Third/AGB/Project/agb2025

echo "Working directory: $(pwd)"
echo "Pipeline output: runs/S01070625"
echo "Ground truth: group4/microbiome-benchmarking/benchmarking_output/synthetic_samples"
echo "ASV mapping: group4/microbiome-benchmarking/input_data/asv_metadata.tsv"
echo "Results will be saved to: group4/microbiome-benchmarking/benchmarking_results/"
echo ""

# Install required packages quietly
echo "Installing required Python packages..."
pip install matplotlib seaborn scikit-learn pandas numpy scipy --quiet

echo "Running genus-level benchmarking..."
echo ""

# Run the benchmarking script
python group4/microbiome-benchmarking/genus_benchmark_v15.py

# Check if it completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "BENCHMARKING COMPLETED SUCCESSFULLY!"
    echo ""
    echo "Results are available in:"
    echo "   group4/microbiome-benchmarking/benchmarking_results/"
    echo ""
    echo "Quick check of results:"
    ls -la group4/microbiome-benchmarking/benchmarking_results/ 2>/dev/null || echo "Results directory not found"
else
    echo ""
    echo "BENCHMARKING FAILED!"
    echo "Check the error messages above."
fi
