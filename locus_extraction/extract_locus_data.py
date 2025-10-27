#!/usr/bin/env python3
"""
Script to extract locus tag data from nested JSON structure and export to Excel.

Usage:
    python extract_locus_data.py -i <json_file> -s <sample> -r <read_length> -p <position> -o <output_excel> [--sort-by <sort_option>]

Example:
    python extract_locus_data.py -i data.json -s TIS-WT-2 -r 27 -p -97 -o output.xlsx
    python extract_locus_data.py -i data.json -s TIS-WT-2 -r 27 -p -97 -o output.xlsx --sort-by count
"""

from pathlib import Path
import traceback

import json
import sys
import argparse
import pandas as pd



def extract_locus_data(json_file, sample, read_length, position):
    """
    Extract locus tag data from JSON file based on specified parameters.

    Args:
        json_file: Path to the JSON file
        sample: Sample name (e.g., 'TIS-WT-2')
        read_length: Read length as string (e.g., '27')
        position: Position as string (e.g., '-97')

    Returns:
        Tuple of (merged_data_dict, chromosome_stats_dict)
        - merged_data: Dictionary with locus_tag as key and read_count as value
        - chromosome_stats: Dictionary with chromosome stats (loci count, read count)
    """
    # Load JSON data
    with open(json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    # Convert parameters to strings to match JSON keys
    read_length = str(read_length)
    position = str(position)

    # Collect data from all chromosomes
    merged_data = {}
    chromosome_stats = {}

    # Iterate through all chromosomes
    for chromosome, chrom_data in data.items():
        # Check if sample exists in this chromosome
        if sample in chrom_data:
            # Check if read_length exists
            if read_length in chrom_data[sample]:
                # Check if position exists
                if position in chrom_data[sample][read_length]:
                    # Get the locus_tag dictionary
                    locus_dict = chrom_data[sample][read_length][position]

                    # Track chromosome statistics
                    num_loci = len(locus_dict)
                    total_reads = sum(locus_dict.values())
                    chromosome_stats[chromosome] = {
                        'num_loci': num_loci,
                        'total_reads': total_reads
                    }

                    # Merge with existing data (sum counts if locus_tag appears in multiple chromosomes)
                    for locus_tag, count in locus_dict.items():
                        if locus_tag in merged_data:
                            merged_data[locus_tag] += count
                        else:
                            merged_data[locus_tag] = count

    return merged_data, chromosome_stats


def parse_gff(gff_file):
    """
    Parse GFF file to extract locus information.

    Args:
        gff_file: Path to GFF file

    Returns:
        Dictionary with structure: {chromosome: {locus_tag: {info_dict}}}
    """
    locus_info = {}

    with open(gff_file, 'r', encoding='utf-8') as f:
        for line in f:
            # Skip comments and empty lines
            if line.startswith('#') or not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chromosome = fields[0]
            feature_type = fields[2]
            start = fields[3]
            stop = fields[4]
            strand = fields[6]
            attributes = fields[8]

            # Only process CDS, gene, or similar features
            if feature_type not in ['CDS', 'gene', 'mRNA', 'tRNA', 'rRNA']:
                continue

            # Parse attributes to extract locus_tag and gene name
            attr_dict = {}
            for attr in attributes.split(';'):
                attr = attr.strip()
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            locus_tag = attr_dict.get('locus_tag', attr_dict.get('ID', ''))
            gene_name = attr_dict.get('gene', attr_dict.get('Name', ''))

            if not locus_tag:
                continue

            # Initialize chromosome if not exists
            if chromosome not in locus_info:
                locus_info[chromosome] = {}

            # Store locus information (update if exists with more complete info)
            if locus_tag not in locus_info[chromosome]:
                locus_info[chromosome][locus_tag] = {
                    'Chromosome': chromosome,
                    'Start': start,
                    'Stop': stop,
                    'Strand': strand,
                    'Locus_tag': locus_tag,
                    'Gene_name': gene_name
                }
            elif gene_name and not locus_info[chromosome][locus_tag]['Gene_name']:
                # Update gene name if it was empty before
                locus_info[chromosome][locus_tag]['Gene_name'] = gene_name

    return locus_info


def create_gff_excel(locus_data, chromosome_stats, gff_info, output_file):
    """
    Create Excel file with locus information from GFF, filtered by loci found in data.
    Creates one sheet per chromosome with relevant loci.

    Args:
        locus_data: Dictionary with locus_tag as key and read_count as value
        chromosome_stats: Dictionary with chromosome statistics
        gff_info: Dictionary with GFF information per chromosome
        output_file: Path to output Excel file
    """
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for chromosome in sorted(chromosome_stats.keys()):
            if chromosome not in gff_info:
                print(f"Warning: Chromosome {chromosome} not found in GFF file")
                continue

            # Get loci for this chromosome that are in our data
            chromosome_loci = []
            for locus_tag in locus_data.keys():
                if locus_tag in gff_info[chromosome]:
                    info = gff_info[chromosome][locus_tag].copy()
                    info['Read_count'] = locus_data[locus_tag]
                    chromosome_loci.append(info)

            if not chromosome_loci:
                print(f"Warning: No matching loci found for chromosome {chromosome}")
                continue

            # Create DataFrame
            df = pd.DataFrame(chromosome_loci)

            # Reorder columns
            column_order = ['Chromosome', 'Start', 'Stop', 'Strand', 'Locus_tag', 'Gene_name', 'Read_count']
            df = df[column_order]

            # Sort by start position
            df['Start'] = pd.to_numeric(df['Start'])
            df = df.sort_values('Start').reset_index(drop=True)

            # Use a shortened sheet name if chromosome name is too long (Excel limit is 31 chars)
            sheet_name = chromosome if len(chromosome) <= 31 else chromosome[:31]

            # Write to Excel
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    return True


def save_to_excel(data_dict, output_file, sort_by='Locus_tag', ascending=True):
    """
    Save the locus tag data to an Excel file.

    Args:
        data_dict: Dictionary with locus_tag as key and read_count as value
        output_file: Path to output Excel file
        sort_by: Column to sort by ('Locus_tag' or 'Read_count')
        ascending: Sort order (True for ascending, False for descending)
    """
    # Convert dictionary to DataFrame
    df = pd.DataFrame(list(data_dict.items()), columns=['Locus_tag', 'Read_count'])

    # Sort by specified column
    df = df.sort_values(sort_by, ascending=ascending).reset_index(drop=True)

    # Save to Excel
    df.to_excel(output_file, index=False, engine='openpyxl')

    return df


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Extract locus tag data from nested JSON structure and export to Excel.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Basic usage (creates simple Locus_tag + Read_count output)
  python extract_locus_data.py -i data.json -s TIS-WT-2 -r 27 -p -97 -o output.xlsx

  # With GFF annotation (creates BOTH the simple output AND a detailed multi-sheet output)
  python extract_locus_data.py -i data.json -s TIS-WT-2 -r 27 -p -97 -o output.xlsx -g annotation.gff --gff-output annotated.xlsx

  # With sorting and verbose output
  python extract_locus_data.py -i data.json -s TIS-WT-2 -r 27 -p -97 -o output.xlsx --sort-by count -v

JSON Structure:
  The script expects a JSON file with the following hierarchy:
  Chromosome -> Sample -> Read-length -> Position -> {Locus_tag: Count}

Output Files:
  -o/--output: Always created. Contains Locus_tag and Read_count columns
  --gff-output: Optional. When specified with --gff, creates a separate multi-sheet Excel
                where each chromosome has its own sheet with:
                Chromosome, Start, Stop, Strand, Locus_tag, Gene_name, Read_count
        """
    )

    parser.add_argument(
        '-i', '--input',
        required=True,
        type=str,
        metavar='FILE',
        help='Path to input JSON file'
    )

    parser.add_argument(
        '-s', '--sample',
        required=True,
        type=str,
        metavar='SAMPLE',
        help='Sample name (e.g., TIS-WT-2)'
    )

    parser.add_argument(
        '-r', '--read-length',
        required=True,
        type=str,
        metavar='LENGTH',
        help='Read length (e.g., 27)'
    )

    parser.add_argument(
        '-p', '--position',
        required=True,
        type=str,
        metavar='POS',
        help='Position (e.g., -97)'
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        type=str,
        metavar='FILE',
        help='Path to output Excel file'
    )

    parser.add_argument(
        '-g', '--gff',
        type=str,
        metavar='FILE',
        help='Optional: Path to GFF file for detailed locus annotation'
    )

    parser.add_argument(
        '--gff-output',
        type=str,
        metavar='FILE',
        help='Optional: Path to output Excel file with GFF annotations (requires --gff). Creates multi-sheet Excel with chromosome info'
    )

    parser.add_argument(
        '--sort-by',
        type=str,
        choices=['locus', 'count'],
        default='locus',
        metavar='SORT',
        help='Sort output by "locus" (alphabetical, default) or "count" (descending)'
    )

    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )

    return parser.parse_args()


def main():
    """Main function to handle command line arguments and execute the extraction."""
    args = parse_arguments()

    # Determine sorting parameters
    if args.sort_by == 'count':
        sort_by = 'Read_count'
        ascending = False  # Descending for counts (highest first)
    else:
        sort_by = 'Locus_tag'
        ascending = True  # Ascending for locus tags (alphabetical)

    # Check if input file exists
    if not Path(args.input).exists():
        print(f"Error: JSON file '{args.input}' not found!")
        sys.exit(1)

    # Check if GFF file exists (if provided)
    if args.gff and not Path(args.gff).exists():
        print(f"Error: GFF file '{args.gff}' not found!")
        sys.exit(1)

    # Check if gff-output is provided without gff
    if args.gff_output and not args.gff:
        print("Error: --gff-output requires --gff to be specified!")
        sys.exit(1)

    try:
        print(f"Extracting data from: {args.input}")
        print("Parameters:")
        print(f"  Sample: {args.sample}")
        print(f"  Read length: {args.read_length}")
        print(f"  Position: {args.position}")
        print(f"  Output: {args.output}")
        if args.gff:
            print(f"  GFF file: {args.gff}")
        if args.gff_output:
            print(f"  GFF output: {args.gff_output}")
        print(f"  Sorting by: {sort_by} ({'descending' if not ascending else 'ascending'})")
        print("=" * 70)

        # Extract data
        locus_data, chromosome_stats = extract_locus_data(
            args.input, args.sample, args.read_length, args.position
        )

        if not locus_data:
            print("\nWarning: No data found for the specified parameters!")
            print("Please check that the sample, read_length, and position exist in the JSON file.")
            sys.exit(1)

        # Display chromosome statistics
        print("\nChromosome Statistics:")
        print("-" * 70)
        total_loci = 0
        total_reads = 0
        for chromosome in sorted(chromosome_stats.keys()):
            stats = chromosome_stats[chromosome]
            num_loci = stats['num_loci']
            num_reads = stats['total_reads']
            total_loci += num_loci
            total_reads += num_reads
            print(f"{chromosome:20s} | Loci: {num_loci:6d} | Reads: {num_reads:10d}")

        print("-" * 70)
        print(f"{'TOTAL':20s} | Loci: {total_loci:6d} | Reads: {total_reads:10d}")
        print("=" * 70)

        # Always create the original output file
        df = save_to_excel(locus_data, args.output, sort_by=sort_by, ascending=ascending)

        print(f"\nSuccess! Found {len(locus_data)} unique locus tags across {len(chromosome_stats)} chromosome(s).")
        print(f"Total read count: {total_reads}")
        print(f"Results saved to: {args.output}")

        if args.verbose:
            print("\nPreview of results (top 10):")
            print(df.head(10).to_string(index=False))
            if len(df) > 10:
                print(f"\n... and {len(df) - 10} more rows (see Excel file for complete data)")

        # Optionally create GFF-annotated output
        if args.gff and args.gff_output:
            print(f"\nParsing GFF file: {args.gff}")
            gff_info = parse_gff(args.gff)
            print(f"Found {sum(len(loci) for loci in gff_info.values())} loci in GFF file")

            create_gff_excel(locus_data, chromosome_stats, gff_info, args.gff_output)
            print("\nGFF-annotated Excel file created!")
            print(f"Results saved to: {args.gff_output}")
            print(f"  - {len(chromosome_stats)} chromosome sheet(s) created")
            print("  - Each sheet contains: Chromosome, Start, Stop, Strand, Locus_tag, Gene_name, Read_count")

    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON file - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()