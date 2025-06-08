#!/bin/bash

# Script to create a run folder for synthetic data
# This script should be saved in group4/microbiome-benchmarking/

# Get the script's directory and project root
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/../.." && pwd )"

# Change to project root to ensure correct paths
cd "$PROJECT_ROOT"

# Create run_id for synthetic data
DATE=$(date +%d%m%y)
RUN_ID="S01${DATE}"  # S for synthetic

# Create the output directory structure
RUN_DIR="runs/${RUN_ID}"
mkdir -p "${RUN_DIR}/raw_data"
mkdir -p "${RUN_DIR}/metadata"

echo "Created run directory: ${RUN_DIR}"

# Copy synthetic samples to raw_data
# Path relative to project root
SYNTHETIC_DIR="group4/microbiome-benchmarking/benchmarking_output/synthetic_samples"

if [ -d "$SYNTHETIC_DIR" ]; then
    echo "Copying synthetic samples from ${SYNTHETIC_DIR}..."
    
    # Check for R1/R2 fastq files and compress them
    if ls ${SYNTHETIC_DIR}/*_R1.fastq 1> /dev/null 2>&1; then
        echo "Found uncompressed .fastq files, compressing..."
        
        # Compress R1 files
        for file in ${SYNTHETIC_DIR}/*_R1.fastq; do
            filename=$(basename "$file" .fastq)
            # Convert R1 to _1 naming convention
            newname="${filename//_R1/_1}"
            echo "Compressing: $file -> ${newname}.fastq.gz"
            gzip -c "$file" > "${RUN_DIR}/raw_data/${newname}.fastq.gz"
        done
        
        # Compress R2 files
        for file in ${SYNTHETIC_DIR}/*_R2.fastq; do
            filename=$(basename "$file" .fastq)
            # Convert R2 to _2 naming convention
            newname="${filename//_R2/_2}"
            echo "Compressing: $file -> ${newname}.fastq.gz"
            gzip -c "$file" > "${RUN_DIR}/raw_data/${newname}.fastq.gz"
        done
        
        echo "Compressed and renamed all .fastq files"
    fi
else
    echo "ERROR: Synthetic samples directory not found at ${SYNTHETIC_DIR}"
    exit 1
fi

# Create a basic metadata file
echo "Creating basic metadata file..."
cat > "${RUN_DIR}/metadata/metadata.tsv" << EOF
sample-id	condition	replicate
EOF

# Add entries for each sample found
cd "${RUN_DIR}/raw_data"
for file in *_1.fastq.gz; do
    if [ -f "$file" ]; then
        sample_name="${file%_1.fastq.gz}"
        echo -e "${sample_name}\tsynthetic\t1" >> "../metadata/metadata.tsv"
    fi
done
cd - > /dev/null

echo "Run setup complete!"
echo "Run ID: ${RUN_ID}"
echo "Files copied to: ${RUN_DIR}/raw_data/"
echo ""
echo "To create metadata file, run:"
echo "  python3 group4/microbiome-benchmarking/create_synthetic_metadata.py runs/${RUN_ID}"
echo ""
echo "After creating metadata, run the pipeline with:"
echo "  nextflow run main.nf --run_id ${RUN_ID} -profile docker"
