#!/usr/bin/env bash

# Download the Kraken2 db

set -euo pipefail
mkdir -p databases
DB_FILENAME="databases/k2_Human_20230629.tar.gz"
TMP_FILENAME="${DB_FILENAME}.part"
EXTRACTED_DIR="databases/k2_Human_20230629"
DB_URL="https://zenodo.org/records/8339700/files/k2_Human_20230629.tar.gz?download=1"

if [[ -f "$DB_FILENAME" ]]; then
  echo "[INFO] File '$DB_FILENAME' already exists. Skipping download."
else
  echo "[INFO] Downloading Kraken2 DB..."

  # Download to a temporary file
  wget --quiet --show-progress \
       --no-check-certificate \
       --output-document="$TMP_FILENAME" \
       "$DB_URL"

  # Rename to final filename
  mv "$TMP_FILENAME" "$DB_FILENAME"
  echo -e "\n[INFO] Download complete. Saved as '$DB_FILENAME'."
fi

if [[ -d "$EXTRACTED_DIR" ]]; then
  echo "[INFO] Directory '$EXTRACTED_DIR' already exists. Skipping extraction."
else
  echo "[INFO] Extracting '$DB_FILENAME'..."
  mkdir "$EXTRACTED_DIR"
  tar -xzf "$DB_FILENAME" -C "$EXTRACTED_DIR"
  echo "[INFO] Extraction complete. Database available at '$EXTRACTED_DIR'."
  # Optional: remove the compressed file to save space
  rm -f "$DB_FILENAME"
  echo "[INFO] Removed compressed file '$DB_FILENAME' after extraction."
fi


# Make sure the databases folder exists
mkdir -p databases

# Download the classifier into databases/
wget -O databases/silva-138-99-nb-classifier.qza \
  https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza


