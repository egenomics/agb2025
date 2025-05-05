# Final full version of the script with:
# - Google Sheets reading
# - One-to-many merge of run and sample metadata
# - Full cleaning (lowercase, underscores, alphanumeric only)
# - Logging to a .log file
# - Hard-stop on validation errors

import pandas as pd
import logging
import re
from datetime import datetime
from google.oauth2 import service_account
from googleapiclient.discovery import build

# --- CONFIG ---
SPREADSHEET_ID = 'YOUR_SPREADSHEET_ID_HERE'
RUN_SHEET = 'Form Responses 1'
SAMPLE_SHEET = 'Form Responses 2'
SERVICE_ACCOUNT_FILE = 'credentials.json'
LOG_FILENAME = f"run_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

# --- LOGGING SETUP ---
logging.basicConfig(
    filename=LOG_FILENAME,
    filemode='w',
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def log_and_print(message, level='info'):
    print(message)
    if level == 'info':
        logging.info(message)
    elif level == 'warning':
        logging.warning(message)
    elif level == 'error':
        logging.error(message)

# --- AUTHENTICATION ---
SCOPES = ['https://www.googleapis.com/auth/spreadsheets.readonly']
credentials = service_account.Credentials.from_service_account_file(
    SERVICE_ACCOUNT_FILE, scopes=SCOPES
)
sheets_api = build('sheets', 'v4', credentials=credentials).spreadsheets()

def read_sheet(sheet_name):
    result = sheets_api.values().get(spreadsheetId=SPREADSHEET_ID, range=sheet_name).execute()
    values = result.get('values', [])
    if not values:
        raise ValueError(f"No data found in sheet: {sheet_name}")
    headers = values[0]
    rows = values[1:]
    return pd.DataFrame(rows, columns=headers)

# --- CLEANING FUNCTION ---
def clean_value(val):
    try:
        if val is None or val == '':
            return val
        if isinstance(val, list):
            return [clean_value(v) for v in val]
        if isinstance(val, str):
            cleaned = val.strip().lower()
            cleaned = re.sub(r'[^a-z0-9]', '_', cleaned)
            cleaned = re.sub(r'_+', '_', cleaned)
            cleaned = cleaned.strip('_')
            return cleaned
        return val
    except Exception as e:
        logging.warning(f"Failed to clean value: {val} — {str(e)}")
        return val

try:
    # --- LOAD DATA ---
    run_df = read_sheet(RUN_SHEET)
    sample_df = read_sheet(SAMPLE_SHEET)

    # --- CLEAN COLUMN NAMES ---
    run_df.columns = [clean_value(col) for col in run_df.columns]
    sample_df.columns = [clean_value(col) for col in sample_df.columns]

    # --- CLEAN ALL VALUES ---
    run_df = run_df.applymap(clean_value)
    sample_df = sample_df.applymap(clean_value)

    # --- VALIDATION 1: Duplicate IDs in run metadata ---
    duplicate_ids = run_df['id'][run_df.duplicated('id')]
    if not duplicate_ids.empty:
        msg = f"⚠️ Duplicate run IDs found: {duplicate_ids.unique().tolist()}"
        log_and_print(msg, 'warning')
        run_df = run_df.drop_duplicates(subset='id')
        log_and_print("Deduplicated run metadata by keeping the first occurrence.")

    # --- MERGE: One run row to many sample rows by 'id' ---
    merged_df = sample_df.merge(run_df, on='id', how='left', indicator=True)

    # --- VALIDATION 2: Check for unmatched samples ---
    unmatched_samples = merged_df[merged_df['_merge'] == 'left_only']
    if not unmatched_samples.empty:
        msg = f"❌ ERROR: Unmatched samples found. IDs: {unmatched_samples['id'].tolist()}"
        log_and_print(msg, 'error')
        log_and_print("🔍 Full unmatched sample rows:\n" + unmatched_samples.to_string(), 'error')
        raise ValueError("Aborted: Some samples have no matching run metadata.")

    # --- FINALIZE AND EXPORT ---
    merged_df.drop(columns=['_merge'], inplace=True)
    output_file = "merged_metadata.csv"
    merged_df.to_csv(output_file, index=False)
    log_and_print(f"✅ Merged metadata saved to '{output_file}'")
    log_and_print(f"📄 Log file saved to '{LOG_FILENAME}'")

except Exception as e:
    log_and_print(f"❌ Script failed: {str(e)}", 'error')
    raise e
