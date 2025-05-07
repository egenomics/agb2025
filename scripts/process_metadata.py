import os
import pandas as pd
import glob
import re
from pathlib import Path

def clean_header(df):
    if 'Marca de temps' in df.columns:
        df = df.drop('Marca de temps', axis=1)
    
    # Replace spaces with underscores in column names
    df.columns = [col.replace(' ', '_') for col in df.columns]
    
    return df

def rename_sample_columns(df):
    column_mapping = {
        'If_you_are_smoker,_specify_the_number_of_cigarettes_per_day': 'Daily_cigarettes',
        'If_you_have_stop_smoking': 'Smoking_cessation_date',
        'For_those_who_consume_daily,_please_specify_the_number_of_drinks_per_day.': 'Daily_alcohol_drinks'
    }
    
    # Apply the column renaming only for columns that exist in the dataframe
    for old_col, new_col in column_mapping.items():
        if old_col in df.columns:
            df = df.rename(columns={old_col: new_col})
    
    return df

def merge_run_and_sample_metadata(folder_path):
    # Find all run and sample metadata files
    run_files = glob.glob(os.path.join(folder_path, "run_metadata*.csv"))
    sample_files = glob.glob(os.path.join(folder_path, "sample_metadata*.csv"))
    
    if not run_files or not sample_files:
        print(f"❌ Missing required files in folder {folder_path}")
        return None
    
    # Load all run metadata files
    run_dfs = []
    for run_file in run_files:
        try:
            df = pd.read_csv(run_file)
            run_dfs.append(df)
        except Exception as e:
            print(f"❌ Error reading {run_file}: {e}")
    
    # Load all sample metadata files
    sample_dfs = []
    for sample_file in sample_files:
        try:
            df = pd.read_csv(sample_file)
            sample_dfs.append(df)
        except Exception as e:
            print(f"❌ Error reading {sample_file}: {e}")
    
    if not run_dfs or not sample_dfs:
        print(f"❌ No valid data files in folder {folder_path}")
        return None
    
    # Combine all run files and all sample files
    combined_run_df = pd.concat(run_dfs, ignore_index=True)
    combined_sample_df = pd.concat(sample_dfs, ignore_index=True)
    
    # Clean headers and remove 'Marca de temps' column
    combined_run_df = clean_header(combined_run_df)
    combined_sample_df = clean_header(combined_sample_df)
    
    # Rename specific columns in sample_metadata
    combined_sample_df = rename_sample_columns(combined_sample_df)
    
    # Merge dataframes based on Run_ID
    merged_df = pd.merge(combined_run_df, combined_sample_df, on='Run_ID', how='inner')
    
    # Save merged file
    merged_file = os.path.join(folder_path, "merged.csv")
    merged_df.to_csv(merged_file, index=False)
    print(f"✅ Created merged file: {merged_file}")
    
    return merged_file

def combine_all_merged_files(merged_files, output_path):
    # Read and combine all dataframes
    dfs = []
    for file in merged_files:
        try:
            df = pd.read_csv(file)
            dfs.append(df)
        except Exception as e:
            print(f"❌ Error reading {file}: {e}")
    
    if not dfs:
        print("❌ No valid merged files found")
        return
    
    # Concatenate all dataframes
    final_df = pd.concat(dfs, ignore_index=True)
    
    # Save the final combined file
    final_df.to_csv(output_path, index=False)
    print(f"✅ Created final metadata file: {output_path}")

def main():
    # Find all fetch_{timestamp} folders in the metadata/original_data directory
    base_directory = "../metadata/original_data"
    fetch_folders = glob.glob(os.path.join(base_directory, "fetch_*"))
    
    # Create merged_metadata directory if it doesn't exist
    merged_directory = "../metadata/merged_metadata"
    os.makedirs(merged_directory, exist_ok=True)
    
    # Process each fetch folder
    merged_files = []
    for folder in fetch_folders:
        result = merge_run_and_sample_metadata(folder)
        if result:
            merged_files.append(result)
    
    # Combine all merged files into a single metadata.csv
    output_file = os.path.join(merged_directory, "metadata.csv")
    combine_all_merged_files(merged_files, output_file)

if __name__ == "__main__":
    main()