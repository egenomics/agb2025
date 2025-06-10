#!/bin/bash

# Path to the CSV file relative to where the script is run
CSV_PATH="group4/datasets/dev.csv"

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
METADATA_DIR="runs/R${NEXT_RUN}${DATE}/metadata"
mkdir -p "$RUN_DIR"
mkdir -p "$METADATA_DIR"
echo "Created: $RUN_DIR"
echo "Created: $METADATA_DIR"

# Generate the metadata file.
# This temporary file will hold the sample IDs that were downloaded.
SAMPLE_IDS_FILE="${METADATA_DIR}/sample_ids.txt"
> "$SAMPLE_IDS_FILE"  # Empty the file if it exists

while read -r sample_id; do
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
      echo "Skipping Sample_ID $sample_id due to unsupported length."
      continue
    fi

    echo "Processing FASTQ files for ${sample_id}..."

    # Append the sample_id the list of samples before downloading
    echo "$sample_id" >> "$SAMPLE_IDS_FILE"

    for read in 1 2; do
      outfile="${RUN_DIR}/${sample_id}_${read}.fastq.gz"
      url="https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/${subdir}/${sample_id}/${sample_id}_${read}.fastq.gz"

      if [[ -f "$outfile" ]]; then
        echo "File $outfile already exists. Skipping download."
      else
        echo "Downloading: $url"
        wget -P "$RUN_DIR" "$url"
      fi
    done
  fi
done < <(tail -n +2 "$CSV_PATH" | cut -d, -f2)

echo -e "\033[31mFASTQ files download finished.\033[0m"

# Generate the metadata.tsv for this run.
CURATED_METADATA="metadata/run_development_dataset/curated/metadata_cleaned.csv"
RUN_METADATA="${METADATA_DIR}/metadata.tsv"

if [[ ! -f "$CURATED_METADATA" ]]; then
  echo "Curated metadata file not found: $CURATED_METADATA"
  exit 1
fi

# Extract header and then filter sample rows using the collected sample_ids,
# converting commas to tabs for TSV output.
header=$(head -n 1 "$CURATED_METADATA" | tr ',' '\t')
echo -e "$header" > "$RUN_METADATA"
tail -n +2 "$CURATED_METADATA" \
  | grep -F -f "$SAMPLE_IDS_FILE" \
  | tr ',' '\t' >> "$RUN_METADATA"

echo -e "\033[31mGenerated metadata file: $RUN_METADATA\033[0m"

# Remove the temporary sample_ids file
rm -f "$SAMPLE_IDS_FILE"

echo -e "\033[31mRun ID for this run: R${NEXT_RUN}${DATE}\033[0m"
