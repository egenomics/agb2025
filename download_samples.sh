#!/bin/bash

# Path to the CSV file relative to where the script is run
CSV_PATH="group4/Dataset/dev.csv"

# Create the output directory
mkdir -p raw_data

# Extract run_accession column and download FASTQ files
tail -n +2 "$CSV_PATH" | cut -d, -f2 | while read -r run_id; do
  # Clean the run_id by stripping quotes and whitespace
  run_id=$(echo "$run_id" | tr -d '"' | xargs)

  if [[ -n "$run_id" ]]; then
    prefix=${run_id:0:6}
    length=${#run_id}

    # Decide folder name depending on run_id length
    if [[ $length -eq 11 ]]; then
      # Use the last two digits and prepend a '0'
      subdir="0${run_id: -2}"
    elif [[ $length -eq 10 ]]; then
      # Use the last digit and prepend '00'
      subdir="00${run_id: -1}"
    else
      echo "Skipping run ID $run_id due to unsupported length."
      continue
    fi

    echo "Processing FASTQ files for ${run_id}..."

    for read in 1 2; do
      outfile="raw_data/${run_id}_${read}.fastq.gz"
      url="https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/${subdir}/${run_id}/${run_id}_${read}.fastq.gz"

      if [[ -f "$outfile" ]]; then
        echo "File $outfile already exists. Skipping download."
      else
        echo "Downloading: $url"
        wget -P raw_data "$url"
      fi
    done
  fi
done