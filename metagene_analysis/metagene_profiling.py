#!/usr/bin/env python
"""
Unified metagene analysis pipeline:
1. Extract metagene windows from annotation with filtering
2. Generate read count JSON files
3. Create metagene profile plots
"""

import argparse
import json

from pathlib import Path

import yaml


from lib.bam_reader import LocusExtractor
from lib.bam_reader import IntervalReader
from lib import io
from lib import misc
from lib import annotation as ann
from lib import plotting
from lib import metagene as mg

def parse_command_line() -> argparse.Namespace:
    """
    Parse the command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Unified metagene analysis: filtering, extraction, and plotting.",
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
    return parser.parse_args()


def extract_read_counts(bam_file, sample_name, config, genome_length_dict):
    """
    Extract read counts for a bam file with its own filtering.
    """
    print(f"\n{'='*60}")
    print(f"Processing BAM file: {sample_name}")
    print(f"{'='*60}")

    # Get mapping modes
    mapping_methods = config.get("mappingMethods", ["fiveprime"])
    bam_file_path = Path(bam_file) if isinstance(bam_file, str) else bam_file

    print("Filtering annotation for this BAM file...")

    ir = IntervalReader(bam_file_path)
    read_intervals_dict, total_counts_dict = ir.output()

    results = {}
    for mapping_method in mapping_methods:
        print(f"  Mapping mode: {mapping_method}")

        # Step 1: Apply filtering specific to annotation
        cds_identifiers, locus_map = ann.filter_annotation_identifiers(
            read_intervals_dict,
            total_counts_dict,
            genome_length_dict,
            mapping_method,
            config
        )

        # Step 2: Extract metagene windows
        metagene_windows_start = mg.extract_start_metagene_windows(cds_identifiers, config)
        metagene_windows_stop = mg.extract_stop_metagene_windows(cds_identifiers, config)

        print(f"Found {len(metagene_windows_start)} start windows and {len(metagene_windows_stop)} stop windows")

        # Step 3: Create interlap dictionaries
        interlap_dict_start = mg.create_interlap_dict(metagene_windows_start, locus_map)
        interlap_dict_stop = mg.create_interlap_dict(metagene_windows_stop, locus_map)

        # Step 4: Extract read counts for each mapping mode
        tmp_dict_start = {}
        tmp_dict_stop = {}

        le = LocusExtractor(
            bam_file_path,
            mapping_method,
            interlap_dict_start,
            interlap_dict_stop,
            config,
            tmp_dict_start,
            tmp_dict_stop
        )
        tmp_dict_start, tmp_dict_stop = le.output()

        results[mapping_method] = {
            'start': tmp_dict_start,
            'stop': tmp_dict_stop,
            'total_counts': {sample_name: total_counts_dict}
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
                fig, _ = plotting.plot_metagene_from_dataframes(
                    df_start_chrom,
                    df_stop_chrom,
                    chrom,
                    read_length_list,
                    color_list
                )

                # Collect figures
                fig_list.append((chrom, mapping_mode, fig))

        # Write all plots for this normalization method after collecting from all mapping methods
        output_formats = config.get("outputFormats", ["html"])
        include_plotlyjs = config.get("includePlotlyJS", True)

        io.write_plots_to_file(fig_list, output_formats, include_plotlyjs, "Metagene Plots", norm_output_dir)

        print(f"  Saved plots for {normalization_method} to {norm_output_dir}")

    print(f"\nAll plots saved to {output_dir}")


def main() -> None:
    """
    Main function for unified metagene analysis.
    """
    args = parse_command_line()

    # Step 1: Load configuration
    print("Loading configuration...")
    with open(args.config_file, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    output_dir = Path(config["outputDirPath"])
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get BAM files and genome
    bam_files = io.parse_bam_file_names(Path(config["alignmentDirPath"]))
    genome_length_dict = io.parse_genome_lengths(config["genomeFilePath"])

    # Get mapping modes
    mapping_methods = config.get("mappingMethods", ["fiveprime"])

    # Initialize combined results for all samples (one per mapping method)
    combined_results = {}
    for mapping_method in mapping_methods:
        combined_results[mapping_method] = {
            'start': {},
            'stop': {},
            'total_counts': {}
        }

    # Step 2 & 3: Process each BAM file independently
    for sample_name, bam_file in bam_files.items():

        # Step 2: Extract read counts for this sample
        sample_results = extract_read_counts(bam_file, sample_name, config, genome_length_dict)

        # Merge this sample's results into combined results
        for mapping_method in mapping_methods:
            # Merge start counts: {chrom: {sample: {readlength: {position: {locus: count}}}}}
            for chrom, sample_data in sample_results[mapping_method]['start'].items():
                if chrom not in combined_results[mapping_method]['start']:
                    combined_results[mapping_method]['start'][chrom] = {}
                # sample_data already has the structure {sample_name: {readlength: ...}}
                combined_results[mapping_method]['start'][chrom].update(sample_data)

            # Merge stop counts
            for chrom, sample_data in sample_results[mapping_method]['stop'].items():
                if chrom not in combined_results[mapping_method]['stop']:
                    combined_results[mapping_method]['stop'][chrom] = {}
                combined_results[mapping_method]['stop'][chrom].update(sample_data)

            # Merge total counts
            combined_results[mapping_method]['total_counts'].update(
                sample_results[mapping_method]['total_counts']
            )

        # Step 3: Create plots and tables for this individual sample
        sample_output_dir = output_dir / sample_name
        sample_output_dir.mkdir(parents=True, exist_ok=True)
        create_metagene_plots(sample_results, config, sample_output_dir)

    # Step 4: Save combined JSON files (all samples together, one per mapping method)
    print("\nSaving combined JSON files...")
    for mapping_method, data in combined_results.items():
        json_file_start = output_dir / f"readcounts_start_{mapping_method}.json"
        json_file_stop = output_dir / f"readcounts_stop_{mapping_method}.json"

        with open(json_file_start, "w", encoding="utf-8") as f:
            json.dump(data['start'], f)

        with open(json_file_stop, "w", encoding="utf-8") as f:
            json.dump(data['stop'], f)

        print(f"  Saved combined JSON for mapping method: {mapping_method}")

    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)


if __name__ == "__main__":
    main()