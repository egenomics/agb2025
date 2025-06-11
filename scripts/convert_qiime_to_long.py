#!/usr/bin/env python3

import pandas as pd
import sys
import argparse
import re

def extract_species_name(taxonomy_string):
    """
    Extract species name in 'Genus species' format from taxonomy string.
    
    Args:
        taxonomy_string (str): Full taxonomy string (e.g., "k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_fragilis")
    
    Returns:
        str: Species name in 'Genus species' format or original if can't parse
    """
    if pd.isna(taxonomy_string) or taxonomy_string == 'Unassigned':
        return 'Unassigned'
    
    try:
        # Split by semicolon and clean up each level
        levels = [level.strip() for level in taxonomy_string.split(';')]
        
        genus = None
        species = None
        
        # Extract genus and species from the taxonomy levels
        for level in levels:
            if level.startswith('g__'):
                genus = level.replace('g__', '').strip()
                # Remove any trailing underscores or special characters
                genus = re.sub(r'[_\s]+$', '', genus)
                
            elif level.startswith('s__'):
                species_full = level.replace('s__', '').strip()
                # Handle different species naming conventions
                if '_' in species_full:
                    # Format: "Bacteroides_fragilis" -> "Bacteroides fragilis"
                    parts = species_full.split('_', 1)
                    if len(parts) == 2:
                        genus_from_species = parts[0]
                        species_epithet = parts[1]
                        # Use genus from species if genus wasn't found earlier
                        if not genus or genus == '':
                            genus = genus_from_species
                        species = species_epithet
                elif ' ' in species_full:
                    # Format: "Bacteroides fragilis" -> already correct
                    parts = species_full.split(' ', 1)
                    if len(parts) == 2:
                        genus_from_species = parts[0]
                        species_epithet = parts[1]
                        if not genus or genus == '':
                            genus = genus_from_species
                        species = species_epithet
                else:
                    # Single word species, might be just epithet
                    species = species_full
        
        # Construct final species name
        if genus and species:
            # Clean up genus and species names
            genus = re.sub(r'[^\w\s-]', '', genus).strip()
            species = re.sub(r'[^\w\s-]', '', species).strip()
            
            if genus and species:
                return f"{genus} {species}"
        
        # Fallback: if we have genus but no species
        if genus and not species:
            genus = re.sub(r'[^\w\s-]', '', genus).strip()
            if genus:
                return f"{genus} sp."
        
        # If we can't parse properly, return the original
        return taxonomy_string
        
    except Exception as e:
        print(f"Warning: Could not parse taxonomy '{taxonomy_string}': {e}")
        return taxonomy_string

def combine_taxonomy_and_features(taxonomy_file, features_file, output_file):
    """
    Combine taxonomy and feature table files into a single output file
    with Sample_ID, abundance, and taxonomy columns.
    
    Args:
        taxonomy_file (str): Path to taxonomy.tsv file
        features_file (str): Path to features_table.tsv file  
        output_file (str): Path for output file
    """
    
    # Read taxonomy file
    print("Reading taxonomy file...")
    taxonomy_df = pd.read_csv(taxonomy_file, sep='\t')
    print(f"Loaded {len(taxonomy_df)} taxonomy entries")
    
    # Read features table
    print("Reading features table...")
    
    # Manual parsing approach to handle BIOM-derived format correctly
    with open(features_file, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    
    # Find the header line (the one that starts with #OTU ID or similar)
    header_line = None
    data_start_idx = 0
    
    for i, line in enumerate(lines):
        if line.startswith('#OTU ID') or (line.startswith('#') and '\t' in line and not line.startswith('# ')):
            header_line = line
            data_start_idx = i + 1
            print(f"Found header at line {i+1}: {line}")
            break
        elif line.startswith('#'):
            print(f"Skipping comment line {i+1}: {line}")
        else:
            # This is a data line, so previous line should have been header
            if i > 0:
                header_line = lines[i-1] if lines[i-1].startswith('#') else None
                data_start_idx = i
            break
    
    if header_line is None:
        raise ValueError("Could not find header line in features table")
    
    # Parse header - remove # if present
    header_parts = header_line.lstrip('#').split('\t')
    feature_id_col = header_parts[0]
    sample_names = header_parts[1:]
    
    print(f"Feature ID column: '{feature_id_col}'")
    print(f"Found {len(sample_names)} samples: {sample_names}")
    
    # Parse data lines
    data = []
    feature_ids = []
    
    for i, line in enumerate(lines[data_start_idx:], data_start_idx + 1):
        if not line or line.startswith('#'):
            continue
            
        parts = line.split('\t')
        if len(parts) != len(header_parts):
            print(f"Warning: Line {i} has {len(parts)} parts, expected {len(header_parts)}")
            print(f"Line content: {line}")
            continue
        
        feature_ids.append(parts[0])
        # Convert abundances to float
        try:
            abundances = [float(x) for x in parts[1:]]
        except ValueError as e:
            print(f"Warning: Could not convert abundances on line {i}: {e}")
            abundances = [0.0] * len(sample_names)
        
        data.append(abundances)
    
    # Create DataFrame
    features_df = pd.DataFrame(data, index=feature_ids, columns=sample_names)
    
    print(f"Successfully parsed features table:")
    print(f"  - {len(features_df)} features")
    print(f"  - {len(features_df.columns)} samples")
    print(f"  - Sample names: {list(features_df.columns)}")
    print(f"  - First few feature IDs: {list(features_df.index[:3])}")
    
    # Create taxonomy lookup dictionary with species name extraction
    print("Processing taxonomy data...")
    taxonomy_lookup = {}
    for _, row in taxonomy_df.iterrows():
        feature_id = row['Feature ID']
        full_taxonomy = row['Taxon']
        species_name = extract_species_name(full_taxonomy)
        taxonomy_lookup[feature_id] = species_name
    
    # Show some examples of the conversion
    print("Example taxonomy conversions:")
    for i, (feature_id, species_name) in enumerate(list(taxonomy_lookup.items())[:5]):
        original_taxonomy = taxonomy_df[taxonomy_df['Feature ID'] == feature_id]['Taxon'].iloc[0]
        print(f"  {feature_id}: '{original_taxonomy}' -> '{species_name}'")
    
    # Prepare output data
    output_data = []
    
    print("Processing data...")
    # Iterate through each sample (column) in features table
    for sample_id in features_df.columns:
        print(f"Processing sample: {sample_id}")
        
        # Get abundances for this sample
        sample_abundances = features_df[sample_id]
        
        # Only include features with non-zero abundance
        non_zero_features = sample_abundances[sample_abundances > 0]
        
        # Add each feature to output
        for feature_id, abundance in non_zero_features.items():
            species_name = taxonomy_lookup.get(feature_id, 'Unassigned')
            
            output_data.append({
                'Sample_ID': sample_id,
                'Species': species_name,
                'Abundance': abundance
            })
    
    # Create output dataframe
    output_df = pd.DataFrame(output_data)
    
    # Sort by Sample_ID and Abundance (descending)
    output_df = output_df.sort_values(['Sample_ID', 'Abundance'], ascending=[True, False])
    
    # Save to file
    print(f"Saving {len(output_df)} records to {output_file}")
    output_df.to_csv(output_file, sep='\t', index=False)
    
    # Print summary statistics
    print("\nSummary:")
    print(f"Total records: {len(output_df)}")
    print(f"Unique samples: {output_df['Sample_ID'].nunique()}")
    print(f"Unique species: {output_df['Species'].nunique()}")
    print(f"Average abundance per record: {output_df['Abundance'].mean():.2f}")
    
    # Show top species by abundance  
    print("\nTop 10 species by total abundance:")
    species_totals = output_df.groupby('Species')['Abundance'].sum().sort_values(ascending=False)
    for species, total_abundance in species_totals.head(10).items():
        print(f"  {species}: {total_abundance:.2f}")
    
    return output_df

def combine_taxonomy_and_features_with_zeros(taxonomy_file, features_file, output_file):
    """
    Version that includes zero abundances - use with caution as this creates much larger files
    """
    
    # Read taxonomy file
    print("Reading taxonomy file...")
    taxonomy_df = pd.read_csv(taxonomy_file, sep='\t')
    
    # Read features file with comment handling
    print("Reading features table...")
    try:
        features_df = pd.read_csv(features_file, sep='\t', comment='#', index_col=0)
    except:
        # Manual parsing if automatic fails
        with open(features_file, 'r') as f:
            lines = [line.strip() for line in f.readlines() if not line.strip().startswith('#')]
        
        header = lines[0].split('\t')
        sample_names = header[1:]
        
        data = []
        feature_ids = []
        
        for line in lines[1:]:
            parts = line.split('\t')
            feature_ids.append(parts[0])
            abundances = [float(x) if x.replace('.', '').replace('-', '').isdigit() else 0.0 for x in parts[1:]]
            data.append(abundances)
        
        features_df = pd.DataFrame(data, index=feature_ids, columns=sample_names)
    
    # Create taxonomy lookup with species name extraction
    print("Processing taxonomy data...")
    taxonomy_lookup = {}
    for _, row in taxonomy_df.iterrows():
        feature_id = row['Feature ID']
        full_taxonomy = row['Taxon']
        species_name = extract_species_name(full_taxonomy)
        taxonomy_lookup[feature_id] = species_name
    
    # Melt the features dataframe to long format
    features_melted = features_df.reset_index().melt(
        id_vars=[features_df.index.name or 'feature_id'], 
        var_name='Sample_ID',
        value_name='abundance'
    )
    
    # Rename the feature ID column for consistency
    feature_id_col_name = features_df.index.name or 'feature_id'
    if feature_id_col_name != 'feature_id':
        features_melted = features_melted.rename(columns={feature_id_col_name: 'feature_id'})
    
    # Add species information
    features_melted['Species'] = features_melted['feature_id'].map(taxonomy_lookup).fillna('Unassigned')
    
    # Reorder columns - only the three you want in the correct order
    output_df = features_melted[['Sample_ID', 'Species', 'abundance']].copy()
    
    # Rename abundance column to match the desired format
    output_df = output_df.rename(columns={'abundance': 'Abundance'})
    
    # Sort by Sample_ID and Abundance
    output_df = output_df.sort_values(['Sample_ID', 'Abundance'], ascending=[True, False])
    
    # Save to file
    print(f"Saving {len(output_df)} records to {output_file}")
    output_df.to_csv(output_file, sep='\t', index=False)
    
    return output_df

def main():
    parser = argparse.ArgumentParser(description='Combine taxonomy and feature table files')
    parser.add_argument('taxonomy_file', help='Path to taxonomy.tsv file')
    parser.add_argument('features_file', help='Path to features_table.tsv file')
    parser.add_argument('output_file', help='Path for output file')
    parser.add_argument('--include-zeros', action='store_true', 
                       help='Include features with zero abundance (default: exclude)')
    
    args = parser.parse_args()
    
    try:
        # Modify function to handle zero abundances if requested
        if args.include_zeros:
            combine_taxonomy_and_features_with_zeros(args.taxonomy_file, 
                                                   args.features_file, 
                                                   args.output_file)
        else:
            combine_taxonomy_and_features(args.taxonomy_file, 
                                        args.features_file, 
                                        args.output_file)
            
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

# Example usage:
# python combine_taxonomy_features.py taxonomy.tsv features_table.tsv combined_output.tsv
# 
# Or to include zero abundances:
# python combine_taxonomy_features.py taxonomy.tsv features_table.tsv combined_output.tsv --include-zeros