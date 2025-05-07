import gspread
from google.oauth2.service_account import Credentials
import pandas as pd
from datetime import datetime
import os

# === CONFIGURATION ===

# Google API scopes
scopes = [
    'https://www.googleapis.com/auth/spreadsheets',
    'https://www.googleapis.com/auth/drive'
]

# Authenticate with service account
creds = Credentials.from_service_account_file('axial-feat-458915-g3-2b5be024ca29.json', scopes=scopes)
gc = gspread.authorize(creds)

# Spreadsheet ID (same for both forms)
spreadsheet_id = "1vjw4acVHFCqRz6lW_AzRxRUV3A7B9m61P6cB7ushjEo"

# Tab (worksheet) names and output file base names
tabs = [
    {"name": "Run_metadata", "output_base": "run_metadata"},
    {"name": "Sample_metadata", "output_base": "sample_metadata"}
]

# Timestamp for filenames and folder
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")

# Create unique folder for this run
base_folder = "../metadata/original_data"
output_folder = os.path.join(base_folder, f"fetch_{timestamp}")
os.makedirs(output_folder, exist_ok=True)

# === FUNCTION TO PROCESS EACH WORKSHEET ===

def process_tab(spreadsheet_id, worksheet_name, output_base):
    try:
        spreadsheet = gc.open_by_key(spreadsheet_id)
        worksheet = spreadsheet.worksheet(worksheet_name)

        records = worksheet.get_all_records()
        df = pd.DataFrame(records)

        # Output file path inside the timestamped folder
        output_csv = os.path.join(output_folder, f"{output_base}.csv")
        df.to_csv(output_csv, index=False)

        num_rows = len(worksheet.get_all_values())
        if num_rows > 1:
            worksheet.batch_clear([f"A2:Z{num_rows}"])

        print(f"✅ Saved: {output_csv}")
    except Exception as e:
        print(f"❌ Error processing '{worksheet_name}': {e}")

# === PROCESS BOTH TABS ===
for tab in tabs:
    process_tab(spreadsheet_id, tab["name"], tab["output_base"])
