process SUGGEST_TRUNCATION_LENGTHS {
    label 'qiime2' // We can reuse the qiime2 environment as it has python

    input:
    path(demux_qzv)

    output:
    path("truncation_suggestion_report.txt")

    script:
    """
    #!/usr/bin/env python
    import zipfile
    import csv
    from pathlib import Path

    demux_qzv = '${demux_qzv}'

    # --- Configuration ---
    # The quality score threshold. We use the 25th percentile (lower quartile).
    # If the quality of 25% of reads at a position drops below this, we suggest truncating before it.
    QUALITY_THRESHOLD = 25
    MIN_OVERLAP = 20 # Expected minimum overlap for DADA2 to merge reads

    # --- Script ---
    print(f"Analyzing {demux_qzv} for truncation suggestions...")
    report_lines = []
    trunc_f, trunc_r = 0, 0

    def find_truncation_position(qzv_path, read_direction):
        with zipfile.ZipFile(qzv_path, 'r') as z:
            # Find the seven-number-summaries.tsv file for the specified read direction
            summary_file_path = None
            for f in z.namelist():
                if f.endswith(f'{read_direction}-seven-number-summaries.tsv'):
                    summary_file_path = f
                    break
            
            if not summary_file_path:
                return None, [f"ERROR: Could not find '{read_direction}-seven-number-summaries.tsv' in the QZV file."]

            report = [f"Parsing quality data from: {summary_file_path}"]
            
            with z.open(summary_file_path) as csvfile:
                reader = csv.reader([line.decode('utf-8') for line in csvfile], delimiter='\\t')
                header = next(reader) # Skip header

                last_good_pos = 0
                for row in reader:
                    position = int(row[0])
                    # Column 2 is the '25%' (lower quartile) quality score
                    lower_quartile_quality = int(float(row[2]))

                    if lower_quartile_quality < QUALITY_THRESHOLD:
                        report.append(f"  - Quality drop detected at position {position}. Lower quartile quality ({lower_quartile_quality}) is below threshold ({QUALITY_THRESHOLD}).")
                        # We suggest truncating at the last position that was good.
                        return last_good_pos, report
                    
                    last_good_pos = position
            
            # If quality never drops, suggest using the full length
            report.append("  - No significant quality drop found. Suggesting full read length.")
            return last_good_pos, report

    # Analyze forward reads
    trunc_f, f_report = find_truncation_position(demux_qzv, 'forward')
    report_lines.extend(f_report)

    # Analyze reverse reads
    trunc_r, r_report = find_truncation_position(demux_qzv, 'reverse')
    report_lines.extend(r_report)

    # --- Generate Final Report ---
    with open("truncation_suggestion_report.txt", "w") as f:
        f.write("=====================================================\\n")
        f.write(" DADA2 Truncation Length Suggestion Report           \\n")
        f.write("=====================================================\\n\\n")

        f.write(f"Method: Find the last base position where the lower quartile (25th percentile) quality score is >= {QUALITY_THRESHOLD}.\\n\\n")
        f.write("--- Analysis Details ---\\n")
        for line in report_lines:
            f.write(f"{line}\\n")
        f.write("\\n")

        f.write("--- Suggestions ---\\n")
        f.write(f"Suggested --p-trunc-len-f: {trunc_f}\\n")
        f.write(f"Suggested --p-trunc-len-r: {trunc_r}\\n\\n")

        f.write("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\\n")
        f.write("!!!                  IMPORTANT                    !!!\\n")
        f.write("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\\n")
        f.write("1. These are only SUGGESTIONS. You are responsible for the final values.\\n")
        f.write(f"2. CRITICAL: Check if your reads will still overlap by at least {MIN_OVERLAP}bp after this truncation.\\n")
        f.write(f"   (e.g., if amplicon is 450bp, {trunc_f} + {trunc_r} > 450 + {MIN_OVERLAP}).\\n")
        f.write("3. Visually inspect the quality plots in the demux.qzv file yourself to confirm these suggestions make sense.\\n\\n")
        f.write("To use these, re-run your pipeline with the following flags (and adjust if needed):\\n")
        f.write(f"   --trunc_len_f {trunc_f} --trunc_len_r {trunc_r}\\n")

    print("Report 'truncation_suggestion_report.txt' created.")
    """
} 