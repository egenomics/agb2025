#!/usr/bin/env python3

import os
import sys
import re

def parse_sample_name(filename):
    """Extract sample info from filename"""
    # Remove _1.fastq.gz or _2.fastq.gz
    base = filename.replace('_1.fastq.gz', '').replace('_2.fastq.gz', '')
    
    # Parse different sample types
    if 'standard_rep' in base:
        condition = 'standard'
        rep = re.search(r'rep_(\d+)', base).group(1) if re.search(r'rep_(\d+)', base) else '1'
    elif 'basic_error' in base:
        condition = 'basic_error'
        rep = base.split('_')[-1]
    elif 'high_quality' in base:
        condition = 'high_quality'
        rep = base.split('_')[-1]
    elif 'no_gc_bias' in base:
        condition = 'no_gc_bias'
        rep = base.split('_')[-1]
    elif 'depth_' in base:
        depth_match = re.search(r'depth_([\d.]+)x_rep_(\d+)', base)
        if depth_match:
            condition = f'depth_{depth_match.group(1)}x'
            rep = depth_match.group(2)
        else:
            condition = 'unknown_depth'
            rep = '1'
    elif 'miseq' in base:
        if 'miseq-24' in base:
            condition = 'miseq_24'
        elif 'miseq-28' in base:
            condition = 'miseq_28'
        else:
            condition = 'miseq_standard'
        rep_match = re.search(r'rep_(\d+)', base)
        rep = rep_match.group(1) if rep_match else '1'
    else:
        condition = 'unknown'
        rep = '1'
    
    return base, condition, rep

def main():
    if len(sys.argv) != 2:
        print("Usage: python create_synthetic_metadata.py <run_directory>")
        print("Example: python group4/microbiome-benchmarking/create_synthetic_metadata.py runs/S01071224")
        sys.exit(1)
    
    # Get the script directory and project root
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
    
    # Change to project root to handle relative paths correctly
    os.chdir(project_root)
    
    run_dir = sys.argv[1]
    raw_data_dir = os.path.join(run_dir, 'raw_data')
    metadata_file = os.path.join(run_dir, 'metadata', 'metadata.tsv')
    
    # Check if raw_data directory exists
    if not os.path.exists(raw_data_dir):
        print(f"Error: Directory {raw_data_dir} does not exist!")
        sys.exit(1)
    
    # Get all _1 files (forward reads)
    samples = []
    for file in sorted(os.listdir(raw_data_dir)):
        if file.endswith('_1.fastq.gz'):
            sample_id, condition, rep = parse_sample_name(file)
            samples.append((sample_id, condition, rep))
    
    if not samples:
        print("Error: No _1.fastq.gz files found in raw_data directory!")
        sys.exit(1)
    
    # Write metadata
    with open(metadata_file, 'w') as f:
        # QIIME2 expects tab-separated format with sample-id as first column
        f.write("sample-id\tcondition\treplicate\n")
        
        # Sample rows
        for sample_id, condition, rep in samples:
            f.write(f"{sample_id}\t{condition}\t{rep}\n")
    
    print(f"Created metadata file with {len(samples)} samples")
    print(f"Metadata saved to: {metadata_file}")
    
    # Print summary
    conditions = {}
    for _, condition, _ in samples:
        conditions[condition] = conditions.get(condition, 0) + 1
    
    print("\nSample summary by condition:")
    for condition, count in sorted(conditions.items()):
        print(f"  {condition}: {count} samples")

if __name__ == "__main__":
    main()
