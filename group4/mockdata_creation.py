import gzip
import os
import random

#bash iss generate     --genomes silva_16s_subset_dna.fasta     --output mock_16s     --n_reads 50000     --model miseq     --abundance lognormal     --gc_bias   

if os.path.exists("mock_data") is False:
    os.makedirs("mock_data")

options = ['A', 'C', 'G', 'T']
random_sequence = ''.join(random.choices(options, k=50))
#Low number of reads
with gzip.open("mock_data/low_read_count.fastq.gz", "wt") as f:
    for i in range(5):  # Create 5 reads
        f.write(f"@read{i}\n")
        f.write(random_sequence + "\n")
        f.write("+\n")
        f.write("I" * 50 + "\n")  # High quality

options2 = ['!', '#', '$', '%']
low_quality = ''.join(random.choices(options, k=50))
#Low quality score overall
with gzip.open("mock_data/low_quality_reads.fastq.gz", "wt") as f:
    for i in range(50):
        f.write(f"@read{i}\n")
        f.write(random_sequence + "\n") 
        f.write("+\n")
        f.write(low_quality+ "\n")  # Lowest qualities
        
random_sequence = ''.join(random.choices(options, k=50))
#Mismatch in number of bases to number of quality scores by 1
with gzip.open("mock_data/slight_mismatch.fastq.gz", "wt") as f:
    for i in range(50):
        f.write(f"@read{i}\n")
        f.write(random_sequence + "\n")  
        f.write("+\n")
        f.write("A" * 49 + "\n")  # One score less

#Sequence is formed by non-dna characters (not ATGC)
with gzip.open("mock_data/incorrect_char_seq.fastq.gz", "wt") as f:
    for i in range(50):
        f.write(f"@read{i}\n")
        f.write("TVPLM" * 10 + "\n")  #Incorrect characters for sequence
        f.write("+\n")
        f.write("A" * 50 + "\n") 

random_sequence = ''.join(random.choices(options, k=50))
#Sequence is RNA (contains U)
with gzip.open("mock_data/rna_seq.fastq.gz", "wt") as f:
    for i in range(50):
        f.write(f"@read{i}\n")
        f.write(random_sequence.replace("T","U") + "\n")  #Sequence contains U
        f.write("+\n")
        f.write("A" * 50 + "\n") 

random_sequence = ''.join(random.choices(options, k=50))
#Corrupted file: Missing quality line
with gzip.open("mock_data/missing_quality.fastq.gz", "wt") as f:
    for i in range(50):
        f.write(f"@read{i}\n")
        f.write(random_sequence + "\n") 
        f.write("+\n")

#Corrupted file: Missing sequence line
with gzip.open("mock_data/missing_sequence.fastq.gz", "wt") as f:
    for i in range(50):
        f.write(f"@read{i}\n")
        f.write("+\n")
        f.write("A" * 50 + "\n") 

random_sequence = ''.join(random.choices(options, k=50))
#Corrupted file: Malformed header (missing @)
with gzip.open("mock_data/malformed_header.fastq.gz", "wt") as f:
    for i in range(50):
        f.write(f"read{i}\n")
        f.write(random_sequence + "\n") 
        f.write("+\n")
        f.write("A" * 50 + "\n") 

random_sequence = ''.join(random.choices(options, k=50))
#Corrupted file: Missing + separator line
with gzip.open("mock_data/no_separator.fastq.gz", "wt") as f:
    for i in range(50):
        f.write(f"@read{i}\n")
        f.write(random_sequence + "\n") 
        f.write("A" * 50 + "\n") 

random_sequence = ''.join(random.choices(options, k=50))
#File is not compressed
with open("mock_data/not_compressed.fastq.gz", "wt") as f:
    for i in range(50):
        f.write(f"@read{i}\n")
        f.write(random_sequence + "\n") 
        f.write("+\n")
        f.write("A" * 50 + "\n") 



