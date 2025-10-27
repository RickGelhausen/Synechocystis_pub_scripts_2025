#!/usr/bin/env python3
"""
Extract sequences from genome based on locus tags with customizable windows.

This script:
1. Reads locus tags and read counts from Excel files
2. Extracts CDS information from GFF annotation
3. Retrieves sequences from genome FASTA file
4. Outputs sequences with customizable windows around start codons
"""

import sys

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_gff(gff_file):
    """
    Parse GFF file and extract CDS features.

    Returns:
        dict: Dictionary mapping locus_tag to feature information
    """
    cds_features = {}

    with open(gff_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            if feature_type != 'CDS':
                continue

            seqid = parts[0]
            start = int(parts[3])  # GFF uses 1-based coordinates
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            # Extract locus_tag from attributes
            locus_tag = None
            for attr in attributes.split(';'):
                if 'locus_tag' in attr:
                    locus_tag = attr.split('=')[-1].strip()
                    break

            if locus_tag:
                cds_features[locus_tag] = {
                    'seqid': seqid,
                    'start': start,
                    'end': end,
                    'strand': strand
                }

    return cds_features


def load_genome(fasta_file):
    """
    Load genome sequences from FASTA file.

    Returns:
        dict: Dictionary mapping sequence IDs to SeqRecord objects
    """
    genome = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        genome[record.id] = record
    return genome


def extract_sequence_with_window(genome, feature, upstream=30, downstream=10):
    """
    Extract sequence with window around start codon.

    Args:
        genome: Dictionary of genome sequences
        feature: Feature dictionary with seqid, start, end, strand
        upstream: Number of nucleotides upstream of start codon
        downstream: Number of nucleotides downstream of start codon

    Returns:
        str: Extracted sequence or None if extraction fails
    """
    seqid = feature['seqid']
    strand = feature['strand']

    if seqid not in genome:
        return None

    seq_record = genome[seqid]
    seq = str(seq_record.seq)

    if strand == '+':
        # For forward strand, start is at feature['start']
        start_codon_pos = feature['start'] - 1  # Convert to 0-based (this is the 'A' of ATG)
        extract_start = max(0, start_codon_pos - upstream)
        extract_end = min(len(seq), start_codon_pos + downstream + 1)  # +1 to include the start codon position itself
        extracted_seq = seq[extract_start:extract_end]
    else:
        # For reverse strand, start codon is at the end (the 'A' of ATG on reverse complement)
        start_codon_pos = feature['end'] - 1  # Convert to 0-based
        extract_start = max(0, start_codon_pos - downstream)
        extract_end = min(len(seq), start_codon_pos + upstream + 1)  # +1 to include the start codon position itself
        extracted_seq = seq[extract_start:extract_end]
        # Reverse complement for minus strand
        extracted_seq = str(Seq(extracted_seq).reverse_complement())

    return extracted_seq


def main():
    """
    Main function to parse arguments and execute sequence extraction.
    """
    parser = argparse.ArgumentParser(
        description='Extract sequences from genome based on locus tags with customizable windows',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -x input.xlsx -g annotation.gff -f genome.fa -o output.fa
  %(prog)s -x input.xlsx -g annotation.gff -f genome.fa -o output.fa -u 50 -d 20 -c 100
        """
    )

    parser.add_argument('-x', '--xlsx', required=True,
                        help='Input Excel file with Locus_tag and Read_count columns')
    parser.add_argument('-g', '--gff', required=True,
                        help='GFF annotation file')
    parser.add_argument('-f', '--fasta', required=True,
                        help='Genome FASTA file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output FASTA file')
    parser.add_argument('-u', '--upstream', type=int, default=30,
                        help='Nucleotides upstream of start codon (default: 30)')
    parser.add_argument('-d', '--downstream', type=int, default=10,
                        help='Nucleotides downstream of start codon (default: 10)')
    parser.add_argument('-c', '--cutoff', type=float, default=0,
                        help='Minimum read count cutoff (default: 0)')

    args = parser.parse_args()

    # Load Excel file
    print(f"Reading Excel file: {args.xlsx}")
    try:
        df = pd.read_excel(args.xlsx)
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        sys.exit(1)

    # Check for required columns
    if 'Locus_tag' not in df.columns or 'Read_count' not in df.columns:
        print("Error: Excel file must contain 'Locus_tag' and 'Read_count' columns")
        print(f"Found columns: {list(df.columns)}")
        sys.exit(1)

    # Filter by read count
    df_filtered = df[df['Read_count'] >= args.cutoff]
    print(f"Found {len(df_filtered)} locus tags with read count >= {args.cutoff}")

    # Parse GFF file
    print(f"Parsing GFF file: {args.gff}")
    cds_features = parse_gff(args.gff)
    print(f"Found {len(cds_features)} CDS features in GFF")

    # Load genome
    print(f"Loading genome: {args.fasta}")
    genome = load_genome(args.fasta)
    print(f"Loaded {len(genome)} sequences from genome")

    # Extract sequences
    print(f"\nExtracting sequences with window: -{args.upstream}nt to +{args.downstream}nt")
    output_records = []
    found = 0
    not_found = 0

    for _, row in df_filtered.iterrows():
        locus_tag = row['Locus_tag']
        read_count = row['Read_count']

        if locus_tag not in cds_features:
            print(f"Warning: Locus tag '{locus_tag}' not found in GFF")
            not_found += 1
            continue

        feature = cds_features[locus_tag]
        sequence = extract_sequence_with_window(
            genome, feature,
            upstream=args.upstream,
            downstream=args.downstream
        )

        if sequence is None:
            print(f"Warning: Could not extract sequence for '{locus_tag}'")
            not_found += 1
            continue

        # Create SeqRecord
        record = SeqRecord(
            Seq(sequence),
            id=locus_tag,
            description=f"read_count={read_count} strand={feature['strand']} window=-{args.upstream}:+{args.downstream}"
        )
        output_records.append(record)
        found += 1

    # Write output
    print(f"\nWriting {found} sequences to: {args.output}")
    SeqIO.write(output_records, args.output, 'fasta')

    print("\nSummary:")
    print(f"  Total locus tags in Excel: {len(df)}")
    print(f"  After read count filter: {len(df_filtered)}")
    print(f"  Successfully extracted: {found}")
    print(f"  Not found/failed: {not_found}")
    print("\nDone!")


if __name__ == '__main__':
    main()