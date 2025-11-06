"""
Contains scripts related to metagene window extraction and interlap dictionary creation.
Author: Rick Gelhausen
"""

import sys
import interlap

def create_interlap_dict(metagene_windows, locus_map) -> dict:
    """
    Create an interlap dictionary for the metagene windows.
    """
    tmp_dict = {}
    for chrom, start, stop, strand, identifier in metagene_windows:
        if identifier in locus_map:
            identifier = locus_map[identifier]

        if (chrom, strand) not in tmp_dict:
            tmp_dict[(chrom, strand)] = [(start, stop, identifier)]
        else:
            tmp_dict[(chrom, strand)].append((start, stop, identifier))

    interlap_dict = {}
    for (chrom, strand), values in tmp_dict.items():
        interlap_dict[(chrom, strand)] = interlap.InterLap(values)

    return interlap_dict

def extract_start_metagene_windows(annotation_identifiers, config) -> list:
    """
    Extract for each annotated CDS the window around each start codon.
    """
    windows = []
    window = [-config["positionsOutsideORF"], config["positionsInsideORF"]]

    if config["positionsOutsideORF"] < 0 or config["positionsInsideORF"] < 0:
        sys.exit("Error: positionsOutsideORF and positionsInsideORF must be non-negative integers.")

    for identifier in annotation_identifiers:
        chrom, start, stop, strand = identifier
        out_id = f"{chrom}:{start+1}-{stop+1}:{strand}"
        cds_length = stop - start

        # Check if positionsInsideORF is longer than the entire CDS
        if window[0] < 0 and abs(window[1]) > cds_length:
            continue

        if strand == "+":
            upstream_region = start + window[0]
            downstream_region = start + window[1]
        else:
            upstream_region = stop - window[1]
            downstream_region = stop - window[0]

        windows.append((chrom, upstream_region, downstream_region, strand, out_id))

    return windows

def extract_stop_metagene_windows(annotation_identifiers, config) -> list:
    """
    Extract for each annotated CDS the window around each stop codon.
    """
    windows = []
    window = [-config["positionsInsideORF"], config["positionsOutsideORF"]]

    if config["positionsOutsideORF"] < 0 or config["positionsInsideORF"] < 0:
        sys.exit("Error: positionsOutsideORF and positionsInsideORF must be non-negative integers.")

    for identifier in annotation_identifiers:
        chrom, start, stop, strand = identifier
        out_id = f"{chrom}:{start+1}-{stop+1}:{strand}"
        cds_length = stop - start

        # Check if downstream is longer than the entire CDS
        if window[0] < 0 and abs(window[0]) > cds_length:
            continue

        if strand == "+":
            upstream_region = stop + window[0]
            downstream_region = stop + window[1]
        else:
            upstream_region = start - window[1]
            downstream_region = start - window[0]

        windows.append((chrom, upstream_region, downstream_region, strand, out_id))

    return windows
