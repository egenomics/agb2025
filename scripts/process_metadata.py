import os
import pandas as pd
import glob

def clean_header(df):
    if 'Marca de temps' in df.columns:
        df = df.drop('Marca de temps', axis=1)
    return df

def format_df(df):
    df = df.map(lambda x: str(x).replace(' ', '_').replace('-', '_').replace('/', '_'))
    df.columns = (df.columns.str.strip().str.lower().str.replace(r'[^\w]+', '_', regex=True))  # Replace non-word chars with underscore
    return df

def rename_sample_columns(df):
    column_mapping = {
        'if_you_are_smoker,_specify_the_number_of_cigarettes_per_day': 'daily_cigarettes',
        'if_you_have_stop_smoking': 'smoking_cessation_date',
        'for_those_who_consume_daily,_please_specify_the_number_of_drinks_per_day.': 'daily_alcohol_drinks'
    }
    
    # Apply the column renaming only for columns that exist in the dataframe
    for old_col, new_col in column_mapping.items():
        if old_col in df.columns:
            df = df.rename(columns={old_col: new_col})
    
    return df

def merge_run_and_sample_metadata(folder_path):
    # Find the run and sample metadata files from the fetch folder
    run_files = glob.glob(os.path.join(folder_path, "run_metadata*.csv"))
    sample_files = glob.glob(os.path.join(folder_path, "sample_metadata*.csv"))

    if not run_files: 
        print(f"❌ Missing required {run_files} in folder {folder_path}")
        return None
    
    elif not sample_files:
        print(f"❌ Missing required {sample_files} in folder {folder_path}")
        return None
    
    try:
        run_df_list = [pd.read_csv(f) for f in run_files]
        run_df = pd.concat(run_df_list, ignore_index=True).drop_duplicates()
    except Exception as e:
        print(f"❌ Error reading run metadata files: {e}")
        return None

    try:
        sample_df_list = [pd.read_csv(f) for f in sample_files]
        sample_df = pd.concat(sample_df_list, ignore_index=True).drop_duplicates()
    except Exception as e:
        print(f"❌ Error reading sample metadata files: {e}")
        return None

    # Clean headers and remove 'Marca de temps' column
    cleaned_run_df = clean_header(run_df)
    cleaned_sample_df = clean_header(sample_df)

    # Normalize column names
    cleaned_run_df.columns = cleaned_run_df.columns.str.lower()
    cleaned_sample_df.columns = cleaned_sample_df.columns.str.lower()

    # Format dataframes to replace spaces and special characters
    formatted_run_df = format_df(cleaned_run_df)
    formatted_sample_df = format_df(cleaned_sample_df)


    # Rename specific columns in sample_metadata
    formatted_sample_df = rename_sample_columns(formatted_sample_df)

    # Merge dataframes based on Run_ID
    merged_df = pd.merge(formatted_run_df, formatted_sample_df, on='run_id', how='inner')
    
    return merged_df

def split_by_run_id(merged_df, output_directory):
    """
    Split a merged metadata CSV file into separate files based on run_ID
    
    Args:
        merged_file (str): Path to the merged metadata CSV file
        output_directory (str): Directory where the split files will be saved
    """
    os.makedirs(output_directory, exist_ok=True)
    
    try:
        for run_id, group_df in merged_df.groupby('run_id'):
            output_file = os.path.join(output_directory, f"run_{run_id}_metadata.csv")
            group_df.to_csv(output_file, index=False)
            print(f"✅ Created metadata file for run {run_id}: {output_file}")
            
    except Exception as e:
        print(f"❌ Error processing {merged_df}: {e}")

def main():
    # Find all fetch_{timestamp} folders in the metadata/original_data directory
    base_directory = "agb2025/metadata/original_data"
    fetch_folders = glob.glob(os.path.join(base_directory, "fetched_*"))
    
    # Create output directory for split files
    split_directory = "agb2025/metadata/merged_data"
    os.makedirs(split_directory, exist_ok=True)
    
    # Process each fetch folder
    for folder in fetch_folders:
        merged_file = merge_run_and_sample_metadata(folder)
        if merged_file is not None:
            split_by_run_id(merged_file, split_directory)

if __name__ == "__main__":
    main()