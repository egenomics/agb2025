#!/bin/bash

# 30-Sample Microbiome Pipeline Benchmarking Script
# Uses ISS abundance generation instead of predefined read counts

# Set variables
GENOMES="input_data/library_gut.fasta"
CPUS=8
BASE_SEED=42

echo "Starting generation of 30 synthetic microbiome samples for benchmarking..."
echo "Estimated total time: 30-45 minutes"

# Check if input files exist
if [ ! -f "$GENOMES" ]; then
    echo "ERROR: Cannot find $GENOMES"
    echo "Make sure you're running from the main project directory"
    exit 1
fi

echo "Found input file:"
echo "  Genomes: $GENOMES"

# Get absolute paths
GENOMES_PATH=$(realpath "$GENOMES")

# Create output directory in organized structure
OUTPUT_DIR="benchmarking_output/synthetic_samples"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# GROUP 1: Standard Replicates (8 samples) - 100K reads each
echo "Generating Group 1: Standard Replicates (8 samples)..."
for i in {1..8}; do
    echo "  Generating standard replicate $i/8..."
    iss generate \
        --genomes "$GENOMES_PATH" \
        --n_reads 100000 \
        --abundance lognormal \
        --model miseq \
        --gc_bias \
        --output standard_rep_${i} \
        --seed $((BASE_SEED + i)) \
        --cpus $CPUS \
        --quiet
done

# GROUP 2: Different Error Models (6 samples) - 100K reads each
echo "Generating Group 2: Different Error Models (6 samples)..."
models=("miseq" "miseq-24" "miseq-28")
for model in "${models[@]}"; do
    for rep in 1 2; do
        echo "  Generating ${model} replicate $rep/2..."
        iss generate \
            --genomes "$GENOMES_PATH" \
            --n_reads 100000 \
            --abundance lognormal \
            --model $model \
            --gc_bias \
            --output ${model}_rep_${rep} \
            --seed $((BASE_SEED + 100 + rep)) \
            --cpus $CPUS \
            --quiet
    done
done

# GROUP 3: Depth Testing (8 samples) - Different read counts
echo "Generating Group 3: Depth Testing (8 samples)..."
declare -A depth_reads=(
    ["0.25x"]=25000
    ["0.5x"]=50000
    ["2.0x"]=200000
    ["5.0x"]=500000
)

for depth in "${!depth_reads[@]}"; do
    n_reads=${depth_reads[$depth]}
    for rep in 1 2; do
        echo "  Generating depth ${depth} replicate $rep/2 (${n_reads} reads)..."
        iss generate \
            --genomes "$GENOMES_PATH" \
            --n_reads $n_reads \
            --abundance lognormal \
            --model miseq \
            --gc_bias \
            --output depth_${depth}_rep_${rep} \
            --seed $((BASE_SEED + 200 + rep)) \
            --cpus $CPUS \
            --quiet
    done
done

# GROUP 4: No GC Bias Controls (3 samples) - 100K reads each
echo "Generating Group 4: No GC Bias Controls (3 samples)..."
for i in {1..3}; do
    echo "  Generating no-GC-bias sample $i/3..."
    iss generate \
        --genomes "$GENOMES_PATH" \
        --n_reads 100000 \
        --abundance lognormal \
        --model miseq \
        --output no_gc_bias_${i} \
        --seed $((BASE_SEED + 300 + i)) \
        --cpus $CPUS \
        --quiet
done

# GROUP 5: Basic Error Model (3 samples) - 100K reads each
echo "Generating Group 5: Basic Error Model (3 samples)..."
for i in {1..3}; do
    echo "  Generating basic error model sample $i/3..."
    iss generate \
        --genomes "$GENOMES_PATH" \
        --n_reads 100000 \
        --abundance lognormal \
        --model miseq \
        --mode basic \
        --output basic_error_${i} \
        --seed $((BASE_SEED + 400 + i)) \
        --cpus $CPUS \
        --quiet
done

# GROUP 6: High-Quality Replicates (2 samples) - 100K reads each
echo "Generating Group 6: High-Quality Replicates (2 samples)..."
for i in {1..2}; do
    echo "  Generating high-quality sample $i/2..."
    iss generate \
        --genomes "$GENOMES_PATH" \
        --n_reads 100000 \
        --abundance lognormal \
        --model miseq-36 \
        --gc_bias \
        --output high_quality_${i} \
        --seed $((BASE_SEED + 500 + i)) \
        --cpus $CPUS \
        --quiet
done

# Generate summary report
echo "Generating sample summary report..."
cat > sample_summary.txt << EOF
MICROBIOME BENCHMARKING DATASET SUMMARY
=====================================
Total Samples: 30
Generation Date: $(date)
Base Seed: $BASE_SEED
Output Location: benchmarking_output/synthetic_samples/
Abundance Generation: ISS log-normal distribution

GROUP BREAKDOWN:
1. Standard Replicates (n=8): 100K reads each
   - standard_rep_1 to standard_rep_8
   - Model: MiSeq with GC bias
   - Purpose: Test reproducibility and baseline performance

2. Different Error Models (n=6): 100K reads each
   - MiSeq replicates (n=2): miseq_rep_1, miseq_rep_2
   - MiSeq-24 replicates (n=2): miseq-24_rep_1, miseq-24_rep_2  
   - MiSeq-28 replicates (n=2): miseq-28_rep_1, miseq-28_rep_2
   - Purpose: Test platform-specific effects

3. Depth Testing (n=8): CRITICAL FOR DIVERSITY ANALYSIS
   - Low depth 25K reads (n=2): depth_0.25x_rep_1, depth_0.25x_rep_2
   - Medium depth 50K reads (n=2): depth_0.5x_rep_1, depth_0.5x_rep_2
   - High depth 200K reads (n=2): depth_2.0x_rep_1, depth_2.0x_rep_2
   - Very high depth 500K reads (n=2): depth_5.0x_rep_1, depth_5.0x_rep_2
   - Purpose: Test rarefaction effects and depth-dependent diversity

4. No GC Bias Controls (n=3): 100K reads each
   - no_gc_bias_1 to no_gc_bias_3
   - Model: MiSeq without GC bias
   - Purpose: Assess GC bias impact

5. Basic Error Model (n=3): 100K reads each
   - basic_error_1 to basic_error_3
   - Model: Basic error model (no indels)
   - Purpose: Test with simpler error structure

6. High-Quality Replicates (n=2): 100K reads each
   - high_quality_1 to high_quality_2
   - Model: MiSeq-36 (high quality) with GC bias
   - Purpose: Test with higher-quality reads

GROUND TRUTH:
- Each sample has its own abundance file generated by ISS
- Diversity metrics calculated from actual ISS-generated abundances
- More realistic than predefined read counts
EOF

# Go back to main directory
cd ../..

echo ""
echo "=========================================="
echo "BENCHMARKING DATASET GENERATION COMPLETE!"
echo "=========================================="
echo "Total samples generated: 30"
echo "Output directory: benchmarking_output/synthetic_samples/"
echo "Summary report: benchmarking_output/synthetic_samples/sample_summary.txt"
