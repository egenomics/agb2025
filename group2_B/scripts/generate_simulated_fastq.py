import random
import argparse
import gzip

def generate_dna_sequence(length):
    """Generates a random DNA sequence of a given length."""
    return ''.join(random.choice('ACGT') for _ in range(length))

def generate_quality_scores(length):
    """Generates a high-quality Phred score string."""
    # Using characters corresponding to high Phred scores (e.g., ASCII 70-74 -> 'F' to 'J')
    # Phred+33 encoding
    return ''.join(random.choice('FGHIJ') for _ in range(length))

def write_fastq_read(file_handle, read_id, sequence, quality):
    """Writes a single FASTQ read to the file handle."""
    file_handle.write(f"@{read_id}\n".encode('utf-8'))
    file_handle.write(f"{sequence}\n".encode('utf-8'))
    file_handle.write("+\n".encode('utf-8'))
    file_handle.write(f"{quality}\n".encode('utf-8'))

def main():
    parser = argparse.ArgumentParser(
        description="Generate simulated, cleaned FASTQ files for microbiome analysis input."
    )
    parser.add_argument(
        "-n", "--num_reads", type=int, required=True,
        help="Number of reads to generate per file."
    )
    parser.add_argument(
        "-l", "--read_length", type=int, default=150,
        help="Length of each simulated read (default: 150)."
    )
    parser.add_argument(
        "-o1", "--output1", required=True,
        help="Output file path for FASTQ R1 (e.g., sample_R1.fastq.gz)."
    )
    parser.add_argument(
        "-o2", "--output2",
        help="Output file path for FASTQ R2 (paired-end). If omitted, only R1 is generated."
    )
    parser.add_argument(
        "--sample_id", default="sim_sample",
        help="Sample ID prefix for read names (default: sim_sample)."
    )

    args = parser.parse_args()

    print(f"Generating {args.num_reads} reads of length {args.read_length}...")

    # Open output files (use gzip for compression)
    file1 = gzip.open(args.output1, 'wb')
    file2 = gzip.open(args.output2, 'wb') if args.output2 else None

    try:
        for i in range(1, args.num_reads + 1):
            read_id_base = f"{args.sample_id}_read_{i}"
            
            # Generate R1
            seq1 = generate_dna_sequence(args.read_length)
            qual1 = generate_quality_scores(args.read_length)
            write_fastq_read(file1, f"{read_id_base}/1", seq1, qual1)

            # Generate R2 if requested
            if file2:
                seq2 = generate_dna_sequence(args.read_length)
                qual2 = generate_quality_scores(args.read_length)
                write_fastq_read(file2, f"{read_id_base}/2", seq2, qual2)

            if i % 10000 == 0:
                print(f"  Generated {i} reads...")

    finally:
        file1.close()
        if file2:
            file2.close()

    print(f"Successfully generated {args.num_reads} reads.")
    print(f"Output R1: {args.output1}")
    if args.output2:
        print(f"Output R2: {args.output2}")

if __name__ == "__main__":
    main() 