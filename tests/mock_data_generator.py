#!/usr/bin/env python3
"""
Enhanced Mock Data Generator for 16S Pipeline Stress Testing
Extends the existing mockdata_creation.py
"""

import gzip
import os
import random
import sys
from pathlib import Path
import argparse

class MockDataGenerator:
    def __init__(self, base_dir="stress_tests"):
        self.base_dir = Path(base_dir)
        self.runs_dir = self.base_dir / "runs"
        
        # Create base runs directory
        self.runs_dir.mkdir(parents=True, exist_ok=True)
    
    def generate_16s_sequence(self, length=250, region="V4"):
        """Generate realistic 16S rRNA sequences with conserved regions"""
        conserved_regions = {
            "V4": [
                "GTGCCAGCMGCCGCGGTAA",   # 515F primer region
                "GGACTACHVGGGTWTCTAAT",  # 806R primer region
                "CCTACGGGNGGCWGCAG",     # Common V4 region
            ],
            "V3V4": [
                "CCTACGGGAGGCAGCAG",     # 341F region
                "GACTACHVGGGTATCTAATCC", # 785R region
            ]
        }
        
        bases = ['A', 'T', 'G', 'C']
        sequence = ""
        
        if region in conserved_regions:
            sequence += random.choice(conserved_regions[region])
        
        while len(sequence) < length:
            if random.random() < 0.1:  
                motifs = ["GAGTTT", "CTGGCT", "GGATCC", "AAGCTT", "GGTACC"]
                sequence += random.choice(motifs)
            else:
                sequence += random.choice(bases)
        
        return sequence[:length]
    
    def generate_quality_scores(self, length, quality_profile="medium", degradation_start=180):
        """Generate quality scores based on profile with degradation starting at specified position"""
        profiles = {
            "high": (35, 40),
            "medium": (25, 35), 
            "low": (10, 25),
            "very_low": (5, 15)
        }
        
        min_qual, max_qual = profiles.get(quality_profile, profiles["medium"])
        quals = []
        
        for i in range(length):
            if i < degradation_start:
                # High quality for positions before degradation_start
                qual = int(random.uniform(min_qual, max_qual))
            else:
                # Quality degrades after degradation_start position
                position_past_threshold = i - degradation_start
                remaining_length = length - degradation_start
                
                if remaining_length > 0:
                    degradation_factor = max(0.3, 1 - (position_past_threshold / remaining_length) * 0.7)
                else:
                    degradation_factor = 0.3
                
                qual = int(random.uniform(min_qual * degradation_factor, max_qual * degradation_factor))
            
            quals.append(chr(min(max(qual, 0), 40) + 33))
        
        return ''.join(quals)
    
    def create_paired_fastq(self, run_id, sample_id, num_reads, quality_profile="medium"):
        """Create paired-end FASTQ files for a sample"""
        run_dir = self.runs_dir / run_id / "raw_data"
        run_dir.mkdir(parents=True, exist_ok=True)
        
        r1_file = run_dir / f"{sample_id}_1.fastq.gz"
        r2_file = run_dir / f"{sample_id}_2.fastq.gz"
        
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        
        with gzip.open(r1_file, 'wt') as f1, gzip.open(r2_file, 'wt') as f2:
            for i in range(num_reads):
                # R1
                seq_id = f"@{sample_id}_{i+1}"
                r1_seq = self.generate_16s_sequence(250)
                r1_qual = self.generate_quality_scores(250, quality_profile, degradation_start=180)
                
                f1.write(f"{seq_id}/1\n{r1_seq}\n+\n{r1_qual}\n")
                
                # R2 
                r2_seq = ''.join([complement.get(base, base) for base in r1_seq[::-1]])
                
                # Realistic mutations (2-5%)
                if quality_profile in ["low", "very_low"]:
                    mutation_rate = random.uniform(0.03, 0.08)
                else:
                    mutation_rate = random.uniform(0.01, 0.03)
                
                r2_list = list(r2_seq)
                num_mutations = int(len(r2_seq) * mutation_rate)
                for _ in range(num_mutations):
                    pos = random.randint(0, len(r2_list) - 1)
                    r2_list[pos] = random.choice(['A', 'T', 'G', 'C'])
                
                r2_seq = ''.join(r2_list)
                r2_qual = self.generate_quality_scores(250, quality_profile, degradation_start=180)
                
                f2.write(f"{seq_id}/2\n{r2_seq}\n+\n{r2_qual}\n")
    
    def create_metadata(self, run_id, filename, samples):
        """Create metadata file for samples"""
        metadata_dir = self.runs_dir / run_id / "metadata"
        metadata_dir.mkdir(parents=True, exist_ok=True)
        metadata_path = metadata_dir / filename
        
        with open(metadata_path, 'w') as f:
            f.write("sample-id\ttreatment\tsubject\ttimepoint\tph\ttemperature\n")
            
            for i, sample in enumerate(samples, 1):
                treatment = random.choice(["Control", "Treatment_A", "Treatment_B", "Treatment_C"])
                subject = f"Subject{i:03d}"
                timepoint = random.choice(["T0", "T1", "T2", "T3"])
                ph = round(random.uniform(6.5, 8.5), 2)
                temp = round(random.uniform(20, 37), 1)
                
                f.write(f"{sample}\t{treatment}\t{subject}\t{timepoint}\t{ph}\t{temp}\n")
    
    
    def generate_stress_test_datasets(self):
        """Generate all stress test datasets"""
        print("Generating stress test datasets...")
        
        # Test 1: Low quality reads
        print("Creating low quality dataset...")
        run_id = "stress_test_low_quality"
        samples = [f"Sample{i:03d}" for i in range(1, 11)]  # 10 samples
        for sample in samples:
            self.create_paired_fastq(run_id, sample, 1000, "very_low")
        self.create_metadata(run_id, "metadata.tsv", samples)
        
        # Test 2: Very low read count
        print("Creating low read count dataset...")
        run_id = "stress_test_low_reads"
        samples = [f"Sample{i:03d}" for i in range(1, 6)]  # 5 samples
        for sample in samples:
            self.create_paired_fastq(run_id, sample, 50, "medium")
        self.create_metadata(run_id, "metadata.tsv", samples)
        
        # Test 3: High sample count
        print("Creating high sample count dataset...")
        run_id = "stress_test_high_samples"
        samples = [f"Sample{i:03d}" for i in range(1, 151)]  # 150 samples
        for sample in samples:
            self.create_paired_fastq(run_id, sample, 5000, "medium")
        self.create_metadata(run_id, "metadata.tsv", samples)
        
        # Test 4: High read depth
        print("Creating high read depth dataset...")
        run_id = "stress_test_high_depth"
        samples = [f"Sample{i:03d}" for i in range(1, 11)]  # 10 samples
        for sample in samples:
            self.create_paired_fastq(run_id, sample, 100000, "high")
        self.create_metadata(run_id, "metadata.tsv", samples)
        
        # Test 5: Mixed quality samples
        print("Creating mixed quality dataset...")
        run_id = "stress_test_mixed_quality"
        samples = [f"Sample{i:03d}" for i in range(1, 21)]  # 20 samples
        qualities = ["high", "medium", "low", "very_low"]
        for i, sample in enumerate(samples):
            quality = qualities[i % len(qualities)]
            reads = random.randint(1000, 10000)
            self.create_paired_fastq(run_id, sample, reads, quality)
        self.create_metadata(run_id, "metadata.tsv", samples)
        
        # Test 6: Standard dataset 
        print("Creating standard dataset...")
        run_id = "stress_test_standard"
        samples = [f"Sample{i:03d}" for i in range(1, 21)]  # 20 samples
        for sample in samples:
            self.create_paired_fastq(run_id, sample, 5000, "medium")
        self.create_metadata(run_id, "metadata.tsv", samples)
        
        # Test 7: Single sample
        print("Creating single sample dataset...")
        run_id = "stress_test_single_sample"
        self.create_paired_fastq(run_id, "Sample001", 10000, "high")
        self.create_metadata(run_id, "metadata.tsv", ["Sample001"])
        
        # Test 8: Medium samples for memory test
        print("Creating medium sample dataset...")
        run_id = "stress_test_limited_memory"
        samples = [f"Sample{i:03d}" for i in range(1, 51)]  # 50 samples
        for sample in samples:
            self.create_paired_fastq(run_id, sample, 5000, "medium")
        self.create_metadata(run_id, "metadata.tsv", samples)
        
        # Test 9: Maximum scale dataset
        print("Creating maximum scale dataset...")
        run_id = "stress_test_max_scale"
        samples = [f"Sample{i:03d}" for i in range(1, 501)]  # 500 samples
        for sample in samples[:10]:  # Only create first 10 for demo
            self.create_paired_fastq(run_id, sample, 20000, "high")
        self.create_metadata(run_id, "metadata.tsv", samples)
        
        self.create_corrupted_files()
        
        print("All stress test datasets created successfully!")
    
    def create_corrupted_files(self):
        """Create corrupted and special test files based on original mockdata_creation.py"""
        # Corrupted files
        corrupted_run = "stress_test_corrupted"
        corrupted_raw_dir = self.runs_dir / corrupted_run / "raw_data"
        corrupted_raw_dir.mkdir(parents=True, exist_ok=True)
        
        # RNA files
        rna_run = "stress_test_rna"
        rna_raw_dir = self.runs_dir / rna_run / "raw_data"
        rna_raw_dir.mkdir(parents=True, exist_ok=True)
        
        # Format validation files
        format_run = "stress_test_format_validation"
        format_raw_dir = self.runs_dir / format_run / "raw_data"
        format_raw_dir.mkdir(parents=True, exist_ok=True)
        
        corrupted_files = [
            ("missing_quality", self.create_missing_quality_file),
            ("missing_sequence", self.create_missing_sequence_file),
            ("malformed_header", self.create_malformed_header_file),
            ("no_separator", self.create_no_separator_file),
        ]
        
        for file_prefix, create_func in corrupted_files:
            create_func(corrupted_raw_dir / f"{file_prefix}_1.fastq.gz")
            create_func(corrupted_raw_dir / f"{file_prefix}_2.fastq.gz")
        
        # RNA sequences
        self.create_rna_files(rna_raw_dir)
        
        # Format validation files
        self.create_format_validation_files(format_raw_dir)
        
        # Create metadata files for each run
        corrupted_samples = ["missing_quality", "missing_sequence", "malformed_header", "no_separator"]
        self.create_metadata(corrupted_run, "metadata.tsv", corrupted_samples)
        self.create_metadata(rna_run, "metadata.tsv", ["RNA_Sample001"])
        self.create_metadata(format_run, "metadata.tsv", ["not_compressed", "incorrect_char"])
    
    def create_missing_quality_file(self, filepath):
        with gzip.open(filepath, 'wt') as f:
            for i in range(50):
                f.write(f"@read{i}\n")
                f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
                f.write("+\n")
                # Missing quality line
    
    def create_missing_sequence_file(self, filepath):
        with gzip.open(filepath, 'wt') as f:
            for i in range(50):
                f.write(f"@read{i}\n")
                # Missing sequence line
                f.write("+\n")
                f.write("A" * 50 + "\n")
    
    def create_malformed_header_file(self, filepath):
        with gzip.open(filepath, 'wt') as f:
            for i in range(50):
                f.write(f"read{i}\n")  # Missing @ symbol
                f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("A" * 50 + "\n")
    
    def create_no_separator_file(self, filepath):
        with gzip.open(filepath, 'wt') as f:
            for i in range(50):
                f.write(f"@read{i}\n")
                f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
                # Missing + separator
                f.write("A" * 50 + "\n")
    
    def create_rna_files(self, rna_dir):
        """Create RNA sequence files"""
        for suffix in ["1", "2"]:
            with gzip.open(rna_dir / f"RNA_Sample001_{suffix}.fastq.gz", 'wt') as f:
                for i in range(100):
                    seq = self.generate_16s_sequence(250).replace('T', 'U')  # Convert to RNA
                    qual = self.generate_quality_scores(250, "medium", degradation_start=180)
                    f.write(f"@RNA_read{i}\n{seq}\n+\n{qual}\n")
    
    def create_format_validation_files(self, format_dir):
    """Create files for format validation testing"""
    # Uncompressed files (should be .fastq.gz but aren't compressed)
    for suffix in ["1", "2"]:
        with open(format_dir / f"not_compressed_{suffix}.fastq.gz", 'w') as f:
            for i in range(50):
                f.write(f"@read{i}/{suffix}\n")
                f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("I" * 50 + "\n")
    
    # Files with incorrect characters in sequence
    for suffix in ["1", "2"]:
        with gzip.open(format_dir / f"incorrect_char_{suffix}.fastq.gz", 'wt') as f:
            for i in range(50):
                f.write(f"@read{i}/{suffix}\n")
                f.write("TVPLMTVPLMTVPLMTVPLMTVPLMTVPLMTVPLMTVPLMTVPLMTVPLM\n")  # Non-DNA chars
                f.write("+\n")
                f.write("I" * 50 + "\n")

def main():
    parser = argparse.ArgumentParser(description="Generate mock data for 16S pipeline stress testing")
    parser.add_argument("--base-dir", default="stress_tests", help="Base directory for test data")
    
    args = parser.parse_args()
    
    generator = MockDataGenerator(args.base_dir)
    generator.generate_stress_test_datasets()

if __name__ == "__main__":
    main()
