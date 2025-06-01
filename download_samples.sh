#!/bin/bash

# Usage message
usage() {
  cat <<EOF
# Usage:
#   ./download_samples.sh [--dev|--test|--healthy-controls]
#
#   --dev              => CSV_PATH="group4/Dataset/dev.csv" and output folder "raw_data"
#   --test             => CSV_PATH="group4/Dataset/test.csv" and output folder "test_data"
#   --healthy-controls => CSV_PATH="group4/Healthy_Controls/final_healthy_controls.csv" and output folder "healthy_controls"
EOF
}

CSV_PATH=""
OUTPUT_DIR=""
COL_NUM=2  # Default column for dev/test

# Parse command-line flags
if [ $# -eq 0 ]; then
  usage
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dev)
      CSV_PATH="group4/Dataset/dev.csv"
      OUTPUT_DIR="raw_data"
      COL_NUM=2
      shift
      ;;
    --test)
      CSV_PATH="group4/Dataset/test.csv"
      OUTPUT_DIR="test_data"
      COL_NUM=2
      shift
      ;;
    --healthy-controls)
      CSV_PATH="group4/Healthy_Controls/final_healthy_controls.csv"
      OUTPUT_DIR="healthy_controls"
      COL_NUM=6
      shift
      ;;
    *)
      echo "Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

# Ensure a CSV path was provided
if [[ -z "$CSV_PATH" ]]; then
  usage
  exit 1
fi

# Create the output directory if needed
mkdir -p "$OUTPUT_DIR"

# For healthy_controls (column 6) and dev/test (column 2).
# Skip the header line with tail -n +2
tail -n +2 "$CSV_PATH" | cut -d, -f"$COL_NUM" | while read -r run_id; do
  # Clean the run_id by stripping quotes and whitespace
  run_id=$(echo "$run_id" | tr -d '"' | xargs)

  if [[ -n "$run_id" ]]; then
    prefix=${run_id:0:6}
    length=${#run_id}

    # Decide folder name depending on run_id length
    if [[ $length -eq 11 ]]; then
      subdir="0${run_id: -2}"
    elif [[ $length -eq 10 ]]; then
      subdir="00${run_id: -1}"
    else
      echo "Skipping run ID $run_id due to unsupported length."
      continue
    fi

    for read in 1 2; do
      outfile="${OUTPUT_DIR}/${run_id}_${read}.fastq.gz"
      url="https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/${subdir}/${run_id}/${run_id}_${read}.fastq.gz"

      if [[ -f "$outfile" ]]; then
        echo "File $outfile already exists. Skipping download."
      else
        echo "Downloading: $url"
        wget -P "$OUTPUT_DIR" "$url"
      fi
    done
  fi
done