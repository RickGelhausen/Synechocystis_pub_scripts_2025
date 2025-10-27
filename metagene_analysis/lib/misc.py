"""
Contains miscellaneous scripts for different calculations.
Author: Rick Gelhausen
"""
import pandas as pd
import numpy as np


def calculate_rpkm(gene_length, read_counts, total_counts):
    """
    Calculate rpkm for a gene
    """

    return (read_counts * 1000000000) / (total_counts * gene_length)

def count_reads(read_intervals_dict, chromosome, strand, beginning, end, mapping_method):
    """
    Select the right reads based on the mapping method.
    """
    read_intervals = read_intervals_dict[(chromosome,strand)].find((beginning, end))
    allowed_reads = []
    for read_interval in read_intervals:
        if mapping_method == "global":
            allowed_reads.append(read_interval)
        elif mapping_method == "threeprime":
            if strand == "+":
                if read_interval[1] <= end:
                    allowed_reads.append(read_interval)
            else:
                if read_interval[0] >= beginning:
                    allowed_reads.append(read_interval)
        elif mapping_method == "fiveprime":
            if strand == "+":
                if read_interval[0] >= beginning:
                    allowed_reads.append(read_interval)
            else:
                if read_interval[1] <= end:
                    allowed_reads.append(read_interval)
        elif mapping_method == "centered":
            center = round((read_interval[0] + read_interval[1]) / 2)
            if center >= beginning and center <= end:
                allowed_reads.append(read_interval)

    return len(allowed_reads)


def window_normalize_df(df, window_size):
    """
    Normalize every read length column by the total read length number and the window size.
    """

    columns = df.columns[1:].tolist()
    for column in columns:
        df[column] = df[column].div( (df[column].sum() / window_size))

    return df


def create_data_frame(metagene_dict, positions_out_ORF, positions_in_ORF, state):
    """
    Create a data frame containing the metagene profiling read counts.
    """

    dataframe_dict = {}
    if state == "start":
        coordinates = list(range(-positions_out_ORF, positions_in_ORF, 1))
    else:
        coordinates = list(range(-positions_in_ORF, positions_out_ORF, 1))

    for chrom in metagene_dict:
        if chrom not in dataframe_dict:
            dataframe_dict[chrom] = pd.DataFrame()
            dataframe_dict[chrom]["coordinates"] = coordinates
        for read_length in sorted(metagene_dict[chrom].keys()):
            dataframe_dict[chrom][f"{read_length}"] = metagene_dict[chrom][read_length]

    return dataframe_dict

def flatten_list(l):
    """
    Flatten a list of lists to a single list.
    """
    return [item for sublist in l for item in sublist]

def equalize_dictionary_keys(start_dict, stop_dict, positions_out_ORF, positions_in_ORF):
    """
    Ensure that both dictionaries have the same set of keys.
    Create new keys for missing values and initialize them with list of 0s.
    """
    metagene_area = np.full(positions_out_ORF + positions_in_ORF, 0)

    unique_keys = set(start_dict.keys()).union(set(stop_dict.keys()))

    start_list = flatten_list([[int(x) for x in start_dict[key].keys()] for key in unique_keys if key in start_dict])
    stop_list = flatten_list([[int(x) for x in stop_dict[key].keys()] for key in unique_keys if key in stop_dict])

    overall_min = min(start_list, stop_list)
    overall_max = max(start_list, stop_list)

    for key in unique_keys:
        if key not in start_dict:
            start_dict[key] = {}
            for read_length in range(overall_min, overall_max + 1):
                start_dict[key][read_length] = metagene_area

        else:
            for read_length in range(overall_min, overall_max + 1):
                if read_length not in start_dict[key]:
                    start_dict[key][read_length] = metagene_area

        if key not in stop_dict:
            stop_dict[key] = {}
            for read_length in range(overall_min, overall_max + 1):
                start_dict[key][read_length] = metagene_area

        else:
            for read_length in range(overall_min, overall_max + 1):
                if read_length not in stop_dict[key]:
                    stop_dict[key][read_length] = metagene_area

    return start_dict, stop_dict


def json_dict_to_dataframe(coverage_dict):
    """
    Convert nested JSON-style dictionary to a flat dataframe.
    Input: {chrom: {sample: {read_length: {position: {gene_id: count}}}}}
    Output: DataFrame with columns [chrom, sample, read_length, position, gene_id, count]
    """
    rows = []
    for chrom, sample_dict in coverage_dict.items():
        for sample, read_length_dict in sample_dict.items():
            for read_length, position_dict in read_length_dict.items():
                for position, gene_dict in position_dict.items():
                    for gene_id, count in gene_dict.items():
                        rows.append({
                            'chrom': chrom,
                            'sample': sample,
                            'read_length': read_length,
                            'position': position,
                            'gene_id': gene_id,
                            'count': count
                        })

    return pd.DataFrame(rows)

def aggregate_coverage_dataframe(df, read_length_list=None):
    """
    Aggregate counts by chromosome, sample, read_length, and position.
    Optionally filter by read lengths.
    """
    if read_length_list is not None:
        df = df[df['read_length'].isin(read_length_list)]

    # Sum counts across all genes for each position
    df_agg = df.groupby(['chrom', 'sample', 'read_length', 'position'])['count'].sum().reset_index()

    return df_agg


def apply_normalization(df, normalization_method, total_counts_dict=None):
    """
    Apply normalization to the dataframe.

    Parameters:
    - df: DataFrame with columns [chrom, sample, read_length, position, count]
    - normalization_method: 'raw', 'cpm', or 'window'
    - total_counts_dict: {sample: {chrom: total_count}} or {sample: total_count} for CPM
    """
    df = df.copy()

    if normalization_method == "raw":
        return df

    elif normalization_method == "cpm":
        if total_counts_dict is None:
            raise ValueError("CPM normalization requires total_counts_dict")

        # Apply CPM normalization per chromosome and sample
        def normalize_cpm(group):
            sample = group['sample'].iloc[0]
            chrom = group['chrom'].iloc[0]

            # Check if total_counts_dict has per-chromosome structure
            if isinstance(total_counts_dict.get(sample), dict):
                total_counts = total_counts_dict[sample].get(chrom, 1)
            else:
                # Fallback to total across all chromosomes
                total_counts = total_counts_dict.get(sample, 1)

            group['count'] = group['count'] / total_counts * 1e6
            return group

        df = df.groupby(['chrom', 'sample'], group_keys=False).apply(normalize_cpm)

    elif normalization_method == "window":
        # Normalize by total counts in the window per chromosome/sample/read_length
        def normalize_window(group):
            total = group['count'].sum()
            window_size = len(group['position'].unique())
            if total > 0:
                group['count'] = group['count'] / (total / window_size)
            return group

        df = df.groupby(['chrom', 'sample', 'read_length'], group_keys=False).apply(normalize_window)

    return df