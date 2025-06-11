#!/usr/bin/env python3
"""
Per-Sample Ground Truth Calculator
Calculates diversity metrics from each synthetic sample's abundance file
This creates the TRUE ground truth for benchmarking your pipeline
"""

import pandas as pd
import numpy as np
import json
import os
from glob import glob
from datetime import datetime

def calculate_diversity_metrics(counts):
    """Calculate diversity metrics from count data"""
    # Remove zero counts
    counts = counts[counts > 0]
    
    if len(counts) == 0:
        return {
            'shannon_diversity': 0,
            'simpson_diversity': 0,
            'inverse_simpson': 0,
            'observed_richness': 0,
            'pielou_evenness': 0,
            'simpson_dominance': 0,
            'total_reads': 0
        }
    
    # Calculate relative abundances
    total_reads = counts.sum()
    rel_abundance = counts / total_reads
    
    # Shannon diversity: H' = -Σ(pi × ln(pi))
    shannon = -np.sum(rel_abundance * np.log(rel_abundance))
    
    # Simpson diversity: D = 1 - Σ(pi²)
    simpson = 1 - np.sum(rel_abundance**2)
    
    # Inverse Simpson (Simpson's reciprocal index)
    inv_simpson = 1 / np.sum(rel_abundance**2)
    
    # Observed richness (number of ASVs with count > 0)
    richness = len(counts)
    
    # Evenness (Pielou's evenness): J = H' / ln(S)
    evenness = shannon / np.log(richness) if richness > 1 else 0
    
    # Dominance (Simpson's dominance index)
    dominance = np.sum(rel_abundance**2)
    
    return {
        'shannon_diversity': float(shannon),
        'simpson_diversity': float(simpson),
        'inverse_simpson': float(inv_simpson),
        'observed_richness': int(richness),
        'pielou_evenness': float(evenness),
        'simpson_dominance': float(dominance),
        'total_reads': int(total_reads)
    }

def parse_insilicoseq_abundance_file(abundance_file):
    """Parse InSilicoSeq abundance file to extract read counts per ASV"""
    print(f"  Parsing {abundance_file}...")
    
    try:
        # InSilicoSeq abundance files are tab-separated with columns:
        # genome_id, abundance, coverage, nb_reads
        df = pd.read_csv(abundance_file, sep='\t')
        print(f"    Columns found: {list(df.columns)}")
        print(f"    Shape: {df.shape}")
        
        # Look for the read count column
        if 'nb_reads' in df.columns:
            asv_ids = df.iloc[:, 0]  # First column should be ASV/genome IDs
            read_counts = df['nb_reads']
        elif df.shape[1] >= 4:
            # Assume last column is read counts
            asv_ids = df.iloc[:, 0]
            read_counts = df.iloc[:, -1]
        else:
            raise ValueError(f"Cannot find read count column in {abundance_file}")
        
        # Create clean dataframe
        result_df = pd.DataFrame({
            'ASV_ID': asv_ids,
            'read_count': read_counts
        })
        
        # Remove zero counts
        result_df = result_df[result_df['read_count'] > 0]
        
        print(f"    Found {len(result_df)} ASVs with reads > 0")
        print(f"    Total reads: {result_df['read_count'].sum():,}")
        
        return result_df
        
    except Exception as e:
        print(f"    ERROR: {e}")
        # Show first few lines for debugging
        try:
            with open(abundance_file, 'r') as f:
                print("    First 5 lines of file:")
                for i, line in enumerate(f):
                    if i < 5:
                        print(f"      {line.strip()}")
                    else:
                        break
        except:
            pass
        raise

def main():
    print("Calculating per-sample ground truth diversity metrics...")
    print("=" * 60)
    
    # Find all abundance files
    abundance_files = glob('*_abundance.txt')
    
    if not abundance_files:
        print("ERROR: No abundance files (*_abundance.txt) found in current directory")
        print("Make sure you're in the benchmarking_samples/ directory")
        return
    
    print(f"Found {len(abundance_files)} abundance files:")
    for f in sorted(abundance_files):
        print(f"  - {f}")
    print()
    
    all_sample_metrics = {}
    sample_dataframes = {}
    
    for abundance_file in sorted(abundance_files):
        sample_name = abundance_file.replace('_abundance.txt', '')
        print(f"Processing sample: {sample_name}")
        
        try:
            # Parse the abundance file
            sample_df = parse_insilicoseq_abundance_file(abundance_file)
            sample_dataframes[sample_name] = sample_df
            
            # Calculate diversity metrics
            metrics = calculate_diversity_metrics(sample_df['read_count'])
            
            # Add metadata
            metrics['sample_name'] = sample_name
            metrics['abundance_file'] = abundance_file
            metrics['calculation_date'] = datetime.now().isoformat()
            
            all_sample_metrics[sample_name] = metrics
            
            print(f"  Results:")
            print(f"    Shannon diversity: {metrics['shannon_diversity']:.4f}")
            print(f"    Simpson diversity: {metrics['simpson_diversity']:.4f}")
            print(f"    Observed richness: {metrics['observed_richness']}")
            print(f"    Total reads: {metrics['total_reads']:,}")
            
            # Save individual sample ground truth files
            with open(f'{sample_name}_ground_truth.json', 'w') as f:
                json.dump(metrics, f, indent=2)
            
            # Save individual sample ASV composition
            sample_df['rel_abundance'] = sample_df['read_count'] / sample_df['read_count'].sum()
            sample_df.to_csv(f'{sample_name}_ground_truth_asvs.csv', index=False)
            
            print(f"    Saved: {sample_name}_ground_truth.json")
            print(f"    Saved: {sample_name}_ground_truth_asvs.csv")
            print()
            
        except Exception as e:
            print(f"  ERROR processing {abundance_file}: {e}")
            print()
            continue
    
    if not all_sample_metrics:
        print("ERROR: No samples were successfully processed")
        return
    
    # Save combined results
    with open('all_samples_ground_truth.json', 'w') as f:
        json.dump(all_sample_metrics, f, indent=2)
    
    # Create summary CSV
    summary_data = []
    for sample_name, metrics in all_sample_metrics.items():
        summary_data.append({
            'sample_name': sample_name,
            'shannon_diversity': metrics['shannon_diversity'],
            'simpson_diversity': metrics['simpson_diversity'],
            'inverse_simpson': metrics['inverse_simpson'],
            'observed_richness': metrics['observed_richness'],
            'pielou_evenness': metrics['pielou_evenness'],
            'simpson_dominance': metrics['simpson_dominance'],
            'total_reads': metrics['total_reads']
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv('ground_truth_summary.csv', index=False)
    
    # Print summary statistics
    print("=" * 60)
    print("SUMMARY STATISTICS ACROSS ALL SAMPLES:")
    print("=" * 60)
    numeric_cols = ['shannon_diversity', 'simpson_diversity', 'observed_richness', 'total_reads']
    for col in numeric_cols:
        values = [metrics[col] for metrics in all_sample_metrics.values()]
        print(f"{col:20}: mean={np.mean(values):.4f}, std={np.std(values):.4f}, min={np.min(values):.4f}, max={np.max(values):.4f}")
    
    print(f"\nFILES CREATED:")
    print(f"  Per-sample ground truth:")
    for sample_name in sorted(all_sample_metrics.keys()):
        print(f"    - {sample_name}_ground_truth.json (diversity metrics)")
        print(f"    - {sample_name}_ground_truth_asvs.csv (ASV abundances)")
    print(f"  Summary files:")
    print(f"    - all_samples_ground_truth.json (all metrics combined)")
    print(f"    - ground_truth_summary.csv (summary table)")
    
    print(f"\nFOR BENCHMARKING:")
    print(f"  1. Run your pipeline on each sample's FASTQ files")
    print(f"  2. For each sample, compare pipeline results vs corresponding _ground_truth files")
    print(f"  3. Calculate accuracy metrics (sensitivity, precision, correlation, etc.)")

if __name__ == "__main__":
    main()
