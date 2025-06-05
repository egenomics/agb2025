process RAREFACTION_THRESHOLD {
    publishDir "${params.outdir}/rarefaction_analysis", mode: 'copy'
    
    input:
    path(feature_table_qza)
    path(metadata_file)
    
    output:
    path("rarefaction_threshold.txt"), emit: threshold_file
    path("rarefaction_summary.txt"), emit: summary
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import subprocess
    import tempfile
    import os
    
    def export_feature_table(qza_file):
        \"\"\"Export QIIME2 feature table to get sample depths\"\"\"
        with tempfile.TemporaryDirectory() as tmpdir:
            # Export feature table
            subprocess.run([
                'qiime', 'tools', 'export',
                '--input-path', qza_file,
                '--output-path', tmpdir
            ], check=True)
            
            # Convert BIOM to TSV
            subprocess.run([
                'biom', 'convert',
                '-i', f'{tmpdir}/feature-table.biom',
                '-o', f'{tmpdir}/feature-table.tsv',
                '--to-tsv'
            ], check=True)
            
            # Load the TSV file
            df = pd.read_csv(f'{tmpdir}/feature-table.tsv', sep='\\t', index_col=0, skiprows=1)
            return df
    
    def calculate_adaptive_threshold(sample_depths, min_retention=${params.min_samples_retained}):
        \"\"\"Calculate adaptive threshold based on sample distribution\"\"\"
        depths = np.array(sample_depths)
        
        # Remove extreme outliers (beyond 3 standard deviations)
        mean_depth = np.mean(depths)
        std_depth = np.std(depths)
        filtered_depths = depths[np.abs(depths - mean_depth) <= 3 * std_depth]
        
        # Try different percentiles and find the one that retains desired fraction
        best_threshold = None
        best_retention = 0
        
        for percentile in [5, 10, 15, 20, 25]:
            threshold = np.percentile(filtered_depths, percentile)
            retention_rate = np.mean(depths >= threshold)
            
            if retention_rate >= min_retention and retention_rate > best_retention:
                best_threshold = threshold
                best_retention = retention_rate
        
        # Fallback: if no threshold meets criteria, use one that retains closest to target
        if best_threshold is None:
            for percentile in [30, 35, 40]:
                threshold = np.percentile(filtered_depths, percentile)
                retention_rate = np.mean(depths >= threshold)
                if retention_rate >= min_retention * 0.9:  # Accept 90% of target
                    best_threshold = threshold
                    best_retention = retention_rate
                    break
        
        # Final fallback: 15th percentile
        if best_threshold is None:
            best_threshold = np.percentile(depths, 15)
            best_retention = np.mean(depths >= best_threshold)
        
        return int(best_threshold), best_retention
    
    def calculate_percentile_threshold(sample_depths, percentile=10):
        \"\"\"Simple percentile-based threshold\"\"\"
        threshold = np.percentile(sample_depths, percentile)
        retention_rate = np.mean(sample_depths >= threshold)
        return int(threshold), retention_rate
    
    def calculate_knee_threshold(sample_depths):
        \"\"\"Find knee point in sorted sample depths\"\"\"
        sorted_depths = np.sort(sample_depths)[::-1]  # Descending order
        
        # Calculate normalized differences
        diffs = np.diff(sorted_depths)
        norm_diffs = diffs / np.max(np.abs(diffs))
        
        # Find the point with maximum change (knee)
        knee_idx = np.argmax(np.abs(norm_diffs))
        threshold = sorted_depths[knee_idx + 1]  # Take the lower value
        retention_rate = np.mean(sample_depths >= threshold)
        
        return int(threshold), retention_rate
    
    # Main execution
    print("Loading feature table...")
    feature_table = export_feature_table('${feature_table_qza}')
    sample_depths = feature_table.sum(axis=0)
    
    print(f"Analyzing {len(sample_depths)} samples")
    print(f"Depth range: {sample_depths.min():.0f} - {sample_depths.max():.0f}")
    print(f"Median depth: {sample_depths.median():.0f}")
    
    # Calculate threshold based on selected method
    method = "${params.rarefaction_method}"
    
    if method == "adaptive":
        threshold, retention = calculate_adaptive_threshold(sample_depths)
        method_desc = "Adaptive threshold based on sample distribution"
    elif method == "percentile":
        threshold, retention = calculate_percentile_threshold(sample_depths, 15)
        method_desc = "15th percentile threshold"
    elif method == "knee":
        threshold, retention = calculate_knee_threshold(sample_depths)
        method_desc = "Knee detection in sample depth distribution"
    else:
        # Default to adaptive
        threshold, retention = calculate_adaptive_threshold(sample_depths)
        method_desc = "Adaptive threshold (default)"
    
    samples_retained = int(retention * len(sample_depths))
    samples_lost = len(sample_depths) - samples_retained
    
    # Write results
    with open('rarefaction_threshold.txt', 'w') as f:
        f.write(str(threshold))
    
    with open('rarefaction_summary.txt', 'w') as f:
        f.write(f"Rarefaction Threshold Analysis\\n")
        f.write(f"==============================\\n\\n")
        f.write(f"Method: {method_desc}\\n")
        f.write(f"Selected threshold: {threshold} reads\\n")
        f.write(f"Samples retained: {samples_retained}/{len(sample_depths)} ({retention:.1%})\\n")
        f.write(f"Samples lost: {samples_lost}\\n\\n")
        f.write(f"Sample depth statistics:\\n")
        f.write(f"  Min: {sample_depths.min():.0f}\\n")
        f.write(f"  Max: {sample_depths.max():.0f}\\n")
        f.write(f"  Mean: {sample_depths.mean():.0f}\\n")
        f.write(f"  Median: {sample_depths.median():.0f}\\n")
        f.write(f"  Std: {sample_depths.std():.0f}\\n")
    
    # Export threshold as environment variable for next processes
    print(f"RAREFACTION_THRESHOLD={threshold}")
    
    print(f"\\nSelected rarefaction threshold: {threshold}")
    print(f"This will retain {samples_retained}/{len(sample_depths)} samples ({retention:.1%})")
    """
    
    stub:
    """
    echo "5000" > rarefaction_threshold.txt
    echo "Stub mode - threshold set to 5000" > rarefaction_summary.txt
    """
}
