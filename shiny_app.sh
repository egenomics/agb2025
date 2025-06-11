#!/bin/bash
# Change to the directory of this script
cd "$(dirname "$0")"

echo "Opening Shiny app... Please wait."

# Run the Shiny app and log all output
Rscript results/shiny_dashboard_results_app.R > results/shiny_log.txt 2>&1

# Check if Rscript was successful
if [ $? -ne 0 ]; then
  echo "❌ ERROR: Failed to launch the Shiny app. Check results/shiny_log.txt for details."
  exit 1
else
  echo "Shiny app launched successfully!"
fi


