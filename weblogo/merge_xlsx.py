#!/usr/bin/env python3
"""
Merge two Excel files containing Locus_tag and Read_count columns.

For each locus_tag:
- Takes the union of all locus_tags from both files
- Keeps the higher read_count when a locus_tag appears in both files
"""

import sys
import argparse
import pandas as pd


def merge_replicates(file1, file2, output_file):
    """
    Merge two Excel files, keeping union of locus_tags and max read_counts.

    Args:
        file1: Path to first Excel file
        file2: Path to second Excel file
        output_file: Path to output Excel file
    """

    # Read both Excel files
    print(f"Reading file 1: {file1}")
    try:
        df1 = pd.read_excel(file1)
    except Exception as e:
        print(f"Error reading {file1}: {e}")
        sys.exit(1)

    print(f"Reading file 2: {file2}")
    try:
        df2 = pd.read_excel(file2)
    except Exception as e:
        print(f"Error reading {file2}: {e}")
        sys.exit(1)

    # Check for required columns in both files
    for df, filename in [(df1, file1), (df2, file2)]:
        if 'Locus_tag' not in df.columns or 'Read_count' not in df.columns:
            print(f"Error: {filename} must contain 'Locus_tag' and 'Read_count' columns")
            print(f"Found columns: {list(df.columns)}")
            sys.exit(1)

    print(f"\nFile 1: {len(df1)} entries")
    print(f"File 2: {len(df2)} entries")

    # Add source column to track where data comes from (for reporting)
    df1['_source'] = 'file1'
    df2['_source'] = 'file2'

    # Concatenate both dataframes
    df_combined = pd.concat([df1, df2], ignore_index=True)

    # Group by Locus_tag and aggregate
    # For Read_count: take the maximum
    # For _source: join to show which files contained each locus_tag
    df_merged = df_combined.groupby('Locus_tag').agg({
        'Read_count': 'max',
        '_source': lambda x: '+'.join(sorted(set(x)))
    }).reset_index()

    # Count overlaps and unique entries
    only_file1 = len(df_merged[df_merged['_source'] == 'file1'])
    only_file2 = len(df_merged[df_merged['_source'] == 'file2'])
    in_both = len(df_merged[df_merged['_source'] == 'file1+file2'])

    print("\nMerge statistics:")
    print(f"  Unique to file 1:    {only_file1}")
    print(f"  Unique to file 2:    {only_file2}")
    print(f"  Present in both:     {in_both}")
    print(f"  Total merged:        {len(df_merged)}")

    # Show some examples of merged entries (where both had the locus_tag)
    if in_both > 0:
        print("\nExample of merged entries (showing first 5 that were in both files):")
        both_mask = df_merged['_source'] == 'file1+file2'
        examples = df_merged[both_mask].head(5)

        for _, row in examples.iterrows():
            locus = row['Locus_tag']
            max_count = row['Read_count']
            count1 = df1[df1['Locus_tag'] == locus]['Read_count'].values[0]
            count2 = df2[df2['Locus_tag'] == locus]['Read_count'].values[0]
            print(f"  {locus}: file1={count1}, file2={count2} -> merged={max_count}")

    # Remove the source column before saving
    df_output = df_merged[['Locus_tag', 'Read_count']].copy()

    # Sort by read count (descending) for easier viewing
    df_output = df_output.sort_values('Read_count', ascending=False).reset_index(drop=True)

    # Save to Excel
    print(f"\nWriting merged file: {output_file}")
    df_output.to_excel(output_file, index=False)

    print(f"\nDone! Merged {len(df_output)} unique locus tags.")


def main():
    """Main function to parse arguments and execute merging of Excel files."""
    parser = argparse.ArgumentParser(
        description='Merge two Excel files with Locus_tag and Read_count columns',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script merges two replicate Excel files by:
1. Taking the union of all locus_tags
2. For locus_tags present in both files, keeping the higher read_count

Examples:
  %(prog)s -1 replicate1.xlsx -2 replicate2.xlsx -o merged.xlsx
  %(prog)s -1 sample_A.xlsx -2 sample_B.xlsx -o combined.xlsx
        """
    )

    parser.add_argument('-1', '--file1', required=True,
                        help='First Excel file (replicate 1)')
    parser.add_argument('-2', '--file2', required=True,
                        help='Second Excel file (replicate 2)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output merged Excel file')

    args = parser.parse_args()

    merge_replicates(args.file1, args.file2, args.output)


if __name__ == '__main__':
    main()