#!/bin/bash
# filepath: /Users/radostina.kisleva/agb2025/download_samples.sh
#!/bin/bash

# Path to the CSV file relative to where the script is run
CSV_PATH="group4/Development_Dataset/development_dataset.csv"

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
      # Fallback (or skip if unexpected length)
      echo "Skipping run ID $run_id due to unsupported length."
      continue
    fi

    echo "Downloading FASTQ files for ${run_id}..."

    wget -P raw_data "https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/${subdir}/${run_id}/${run_id}_1.fastq.gz"
    wget -P raw_data "https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/${subdir}/${run_id}/${run_id}_2.fastq.gz"
  fi
done