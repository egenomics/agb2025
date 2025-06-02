#!/usr/bin/env bash

set -euo pipefail

DB_FILENAME="k2_Human_20230629.tar.gz"
TMP_FILENAME="${DB_FILENAME}.part"
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