process RAREFACTION_THRESHOLD {
    label 'qiime2'

    input:
    path table_qza
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
    
    # Wrap the entire script in a try-except to catch and display errors
    try:
        def export_feature_table(qza_file):
            with tempfile.TemporaryDirectory() as tmpdir:
                subprocess.run([
                    'qiime', 'tools', 'export',
                    '--input-path', qza_file,
                    '--output-path', tmpdir
                ], check=True)
                
                subprocess.run([
                    'biom', 'convert',
                    '-i', f'{tmpdir}/feature-table.biom',
                    '-o', f'{tmpdir}/feature-table.tsv',
                    '--to-tsv'
                ], check=True)
                
                df = pd.read_csv(f'{tmpdir}/feature-table.tsv', sep='\\t', index_col=0, skiprows=1)
                return df

        print("Loading feature table...")
        feature_table = export_feature_table('${table_qza}')
        sample_depths = feature_table.sum(axis=0)
        
        print(f"Analyzing {len(sample_depths)} samples")
        
        # Simple percentile-based threshold
        threshold = int(np.percentile(sample_depths, 10))
        retention = np.mean(sample_depths >= threshold)
        
        # Write results
        with open('rarefaction_threshold.txt', 'w') as f:
            f.write(str(threshold))
        
        with open('rarefaction_summary.txt', 'w') as f:
            f.write(f"Rarefaction Threshold Analysis\\n")
            f.write(f"Selected threshold: {threshold}\\n")
            f.write(f"Samples retained: {int(retention * len(sample_depths))}/{len(sample_depths)}\\n")
        
        print(f"Successfully wrote threshold: {threshold}")
        
    except Exception as e:
        print(f"Error in Python script: {str(e)}")
        raise
    """

    stub:
    """
    echo "5000" > rarefaction_threshold.txt
    echo "Stub mode - threshold set to 5000" > rarefaction_summary.txt
    """
}