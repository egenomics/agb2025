#!/bin/bash

# Path to the CSV file relative to where the script is run
CSV_PATH="group4/Dataset/dev.csv"

# Create run_id to make new directory inside runs/
DATE=$(date +%d%m%y)
LAST_RUN=$(ls -1d runs/R*${DATE} 2>/dev/null | sed "s/runs\/R\([0-9]\{2\}\)${DATE}/\1/" | sort -n | tail -1)

# If no runs that day, start with 01, otherwise increment by 1
if [[ -z "$LAST_RUN" ]]; then
    NEXT_RUN="01"
else
    NEXT_RUN=$(printf "%02d" $((10#$LAST_RUN + 1)))
fi

# Create the output directory
RUN_DIR="runs/R${NEXT_RUN}${DATE}/raw_data"
mkdir -p "$RUN_DIR"
echo "Created: $RUN_DIR"

# Extract run_accession column and download FASTQ files
tail -n +2 "$CSV_PATH" | cut -d, -f2 | while read -r sample_id; do
  # Clean the sample_id by stripping quotes and whitespace
  sample_id=$(echo "$sample_id" | tr -d '"' | xargs)

  if [[ -n "$sample_id" ]]; then
    prefix=${sample_id:0:6}
    length=${#sample_id}

    # Decide folder name depending on sample_id length
    if [[ $length -eq 11 ]]; then
      # Use the last two digits and prepend a '0'
      subdir="0${sample_id: -2}"
    elif [[ $length -eq 10 ]]; then
      # Use the last digit and prepend '00'
      subdir="00${sample_id: -1}"
    else
      echo "Skipping sample ID $sample_id due to unsupported length."
      continue
    fi

    echo "Processing FASTQ files for ${sample_id}..."

    for read in 1 2; do
      outfile="${sample_id}_${read}.fastq.gz"
      url="https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/${subdir}/${sample_id}/${sample_id}_${read}.fastq.gz"

      if [[ -f "$outfile" ]]; then
        echo "File $outfile already exists. Skipping download."
      else
        echo "Downloading: $url"
        wget -P $RUN_DIR "$url"
      fi
    done
  fi
done
