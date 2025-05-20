# Define input and output file paths at the top
input_path <- "../metadata/merged_metadata/metadata.csv"
output_path <- "../metadata/merged_metadata/metadata_reordered.csv"

# Load necessary libraries
library(readr)
library(dplyr)
library(stringr)

# Function to reorder columns and add new columns
process_metadata <- function(input_file, output_file) {
  # Read the CSV file
  data <- read_csv(input_file, show_col_types = FALSE)
  
  # Print column names for verification
  cat("Original column names in the CSV file:\n")
  print(colnames(data))
  
  # Get all column names
  all_cols <- colnames(data)
  
  # Step 1: Rename Notes column to Notes_Samples
  if ("Notes" %in% all_cols) {
    data <- data %>% rename(Notes_Samples = Notes)
    # Update all_cols to reflect the name change
    all_cols[all_cols == "Notes"] <- "Notes_Samples"
  }
  
  # Find the index of columns
  run_id_col <- "Run_ID"
  technician_name_col <- "Technician_name"
  sample_id_col <- "Sample_ID"
  
  run_id_index <- which(all_cols == run_id_col)
  technician_name_index <- which(all_cols == technician_name_col)
  sample_id_index <- which(all_cols == sample_id_col)
  
  # If the columns are not found, stop with an error
  if (length(run_id_index) == 0 || length(technician_name_index) == 0 || length(sample_id_index) == 0) {
    stop(paste("Could not find", run_id_col, ",", technician_name_col, "or", sample_id_col, "columns in the CSV file."))
  }
  
  # Step 2: Create the RunHM column based on the formula
  # RunHM = R + sample_ID[2:3] + Run_ID[-5:] + sample_ID[4:5] + sequence_number
  
  # Create a function to generate sequence number for each Run_ID
  generate_sequence <- function(run_id) {
    # Group by Run_ID and create a sequence number for each row
    data %>%
      group_by(!!sym(run_id_col)) %>%
      mutate(seq_num = sprintf("%02d", row_number())) %>%
      pull(seq_num)
  }
  
  # Generate sequence numbers
  seq_numbers <- generate_sequence(run_id_col)
  
  # Create the RunHM column
  data <- data %>%
    mutate(RunHM = paste0(
      "R",
      substr(!!sym(sample_id_col), 2, 3),  # 2nd and 3rd letters of Sample_ID
      substr(!!sym(run_id_col), nchar(!!sym(run_id_col)) - 4, nchar(!!sym(run_id_col))),  # Last 5 chars of Run_ID
      "_",  # First underscore
      substr(!!sym(sample_id_col), 4, 5),  # 4th and 5th position of Sample_ID
      "_",  # Second underscore
      seq_numbers  # Sequence number
    ))
  
  # Identify columns to move (from Run_ID to Technician_name, inclusive)
  cols_to_move <- all_cols[run_id_index:technician_name_index]
  
  # Identify columns that remain in place (after Technician_name)
  cols_after <- all_cols[(technician_name_index+1):length(all_cols)]
  
  # Create the new column order: first cols_after, then RunHM, then cols_to_move
  new_col_order <- c(cols_after, "RunHM", cols_to_move)
  
  # Reorder the data frame columns
  data_reordered <- data %>% select(all_of(new_col_order))
  
  # Step 3: Add Notes_Runs column at the end
  data_reordered <- data_reordered %>% 
    mutate(Notes_Runs = NA_character_)
  
  # Print the final column names for verification
  cat("\nFinal column names in the reordered CSV file:\n")
  print(colnames(data_reordered))
  
  # Write the reordered data to a new CSV file
  write_csv(data_reordered, output_file)
  
  cat("\nProcessing complete! File has been saved to:", output_file, "\n")
  return(data_reordered)
}

# Execute the function
tryCatch({
  processed_data <- process_metadata(input_path, output_path)
}, error = function(e) {
  cat("Error:", conditionMessage(e), "\n")
  cat("Please check column names and data structure\n")
})