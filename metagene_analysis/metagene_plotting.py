#!/usr/bin/env python
"""
Standalone script to create metagene plots from existing JSON files.
No extraction needed - just reads the JSON and generates plots/tables.

This code was created using AI assistance, based on the original plotting functionality to provide a faster plotting .
"""

import argparse
import json
from pathlib import Path
import yaml

from lib import io
import lib.misc as misc
import lib.plotting as plotting


def parse_command_line() -> argparse.Namespace:
    """
    Parse the command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Create metagene plots from existing JSON files.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-c",
        "--config",
        dest="config_file",
        required=True,
        type=Path,
        help="Configuration file (YAML format)",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_dir",
        required=True,
        type=Path,
        help="Input directory containing JSON files",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        required=True,
        type=Path,
        help="Output directory for plots and tables",
    )
    return parser.parse_args()


def load_json_files(input_dir, mapping_methods):
    """
    Load JSON files for all mapping methods.
    Returns: {mapping_method: {'start': dict, 'stop': dict}}
    """
    results = {}

    for mapping_method in mapping_methods:
        json_file_start = input_dir / f"readcounts_start_{mapping_method}.json"
        json_file_stop = input_dir / f"readcounts_stop_{mapping_method}.json"

        if not json_file_start.exists() or not json_file_stop.exists():
            print(f"Warning: JSON files for {mapping_method} not found, skipping...")
            continue

        print(f"Loading JSON files for {mapping_method}...")

        with open(json_file_start, "r", encoding="utf-8") as f:
            start_dict = json.load(f)

        with open(json_file_stop, "r", encoding="utf-8") as f:
            stop_dict = json.load(f)

        # Extract sample names and calculate total counts
        total_counts = {}
        for chrom, sample_dict in start_dict.items():
            for sample, read_length_dict in sample_dict.items():
                if sample not in total_counts:
                    total_counts[sample] = {}
                if chrom not in total_counts[sample]:
                    total_counts[sample][chrom] = 0

                for read_length, position_dict in read_length_dict.items():
                    for position, locus_dict in position_dict.items():
                        total_counts[sample][chrom] += sum(locus_dict.values())

        results[mapping_method] = {
            'start': start_dict,
            'stop': stop_dict,
            'total_counts': total_counts
        }

    return results


def create_metagene_plots(results, config, output_dir):
    """
    Create metagene profile plots from coverage dictionaries.
    """
    print("Creating metagene plots...")

    read_length_list = io.parse_read_lengths(config.get("readLengths", "20-35"))
    color_list = config.get("colorList", None)
    normalization_methods = config.get("normalizationMethods", ["raw"])

    # Process each normalization method first (outer loop)
    for normalization_method in normalization_methods:
        normalization_method = normalization_method.lower()
        print(f"\nProcessing normalization method: {normalization_method}")

        # Create output subdirectory for this normalization
        norm_output_dir = output_dir / normalization_method
        norm_output_dir.mkdir(parents=True, exist_ok=True)

        # Collect figures across all mapping methods
        fig_list = []

        # Process each mapping mode
        for mapping_mode, data in results.items():
            print(f"  Processing mapping mode: {mapping_mode}")

            # Convert JSON dictionaries to dataframes
            df_start = misc.json_dict_to_dataframe(data['start'])
            df_stop = misc.json_dict_to_dataframe(data['stop'])

            # Aggregate by position (sum across genes)
            df_start_agg_all = misc.aggregate_coverage_dataframe(df_start, read_length_list=None)  # All read lengths
            df_stop_agg_all = misc.aggregate_coverage_dataframe(df_stop, read_length_list=None)   # All read lengths

            # Also create filtered version for plotting
            df_start_agg_filtered = misc.aggregate_coverage_dataframe(df_start, read_length_list)  # Only requested
            df_stop_agg_filtered = misc.aggregate_coverage_dataframe(df_stop, read_length_list)   # Only requested

            # Apply normalization
            total_counts_dict = data.get('total_counts', None)

            # Normalize ALL read lengths for Excel
            df_start_norm_all = misc.apply_normalization(df_start_agg_all, normalization_method, total_counts_dict)
            df_stop_norm_all = misc.apply_normalization(df_stop_agg_all, normalization_method, total_counts_dict)

            # Normalize FILTERED read lengths for plotting
            df_start_norm_filtered = misc.apply_normalization(df_start_agg_filtered, normalization_method, total_counts_dict)
            df_stop_norm_filtered = misc.apply_normalization(df_stop_agg_filtered, normalization_method, total_counts_dict)

            # Save Excel files with ALL read lengths (one sheet per chromosome)
            io.save_dataframes_to_excel(df_start_norm_all, norm_output_dir / f"{mapping_mode}_readcounts_start.xlsx")
            io.save_dataframes_to_excel(df_stop_norm_all, norm_output_dir / f"{mapping_mode}_readcounts_stop.xlsx")


            # Create plots for each chromosome using FILTERED read lengths
            for chrom in df_start_norm_filtered['chrom'].unique():
                df_start_chrom = df_start_norm_filtered[df_start_norm_filtered['chrom'] == chrom]
                df_stop_chrom = df_stop_norm_filtered[df_stop_norm_filtered['chrom'] == chrom]

                # Create combined plot
                fig, max_y = plotting.plot_metagene_from_dataframes(
                    df_start_chrom,
                    df_stop_chrom,
                    chrom,
                    read_length_list,
                    color_list
                )

                # Collect figures (don't write yet)
                fig_list.append((chrom, mapping_mode, fig))

        # Write all plots for this normalization method after collecting from all mapping methods
        output_formats = config.get("outputFormats", ["html"])
        include_plotlyjs = config.get("includePlotlyJS", True)

        io.write_plots_to_file(fig_list, output_formats, include_plotlyjs, "metagene", norm_output_dir)

        print(f"  Saved plots for {normalization_method} to {norm_output_dir}")

    print(f"\nAll plots saved to {output_dir}")


def main() -> None:
    """
    Main function to create plots from existing JSON files.
    """
    args = parse_command_line()

    # Load configuration
    print("Loading configuration...")
    with open(args.config_file, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    # Get mapping methods from config
    mapping_methods = config.get("mappingMethods", ["fiveprime"])

    # Load JSON files (contains ALL samples)
    results = load_json_files(args.input_dir, mapping_methods)

    if not results:
        print("Error: No valid JSON files found!")
        return

    # Extract all unique samples
    all_samples = set()
    for mapping_mode, data in results.items():
        for chrom, sample_dict in data['start'].items():
            all_samples.update(sample_dict.keys())

    print(f"Found samples: {sorted(all_samples)}")

    # Process EACH SAMPLE separately
    for sample_name in sorted(all_samples):
        print(f"\n{'='*60}")
        print(f"Processing sample: {sample_name}")
        print(f"{'='*60}")

        # Filter results to only include THIS sample
        sample_results = {}
        for mapping_mode, data in results.items():
            sample_start = {}
            sample_stop = {}

            # Extract only this sample's data
            for chrom in data['start'].keys():
                if sample_name in data['start'][chrom]:
                    if chrom not in sample_start:
                        sample_start[chrom] = {}
                    sample_start[chrom][sample_name] = data['start'][chrom][sample_name]

            for chrom in data['stop'].keys():
                if sample_name in data['stop'][chrom]:
                    if chrom not in sample_stop:
                        sample_stop[chrom] = {}
                    sample_stop[chrom][sample_name] = data['stop'][chrom][sample_name]

            # Get total counts for this sample
            sample_total_counts = {}
            if 'total_counts' in data and sample_name in data['total_counts']:
                sample_total_counts[sample_name] = data['total_counts'][sample_name]

            sample_results[mapping_mode] = {
                'start': sample_start,
                'stop': sample_stop,
                'total_counts': sample_total_counts
            }

        # Create output directory for this sample
        sample_output_dir = args.output_dir / sample_name
        sample_output_dir.mkdir(parents=True, exist_ok=True)

        # Call the plotting function with ONLY this sample's data
        create_metagene_plots(sample_results, config, sample_output_dir)

    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()