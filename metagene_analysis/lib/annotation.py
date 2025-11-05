"""
Contains scripts related to annotation processing.
Author: Rick Gelhausen
"""

import interlap
import pandas as pd

from lib import misc


def create_annotation_intervals_dict(annotation_df):
    """
    Create interlap instances for each chrom / strand in the annotation
    """

    annotation_intervals_dict = {}

    tmp_dict = {}
    for row in annotation_df.itertuples(index=False):
        chromosome = row[0]
        beginning = int(row[3]) - 1
        end = int(row[4]) - 1
        strand = row[6]
        feature = row[2]

        if feature.lower() != "cds":
            continue

        if (chromosome, strand) not in tmp_dict:
            tmp_dict[(chromosome, strand)] = [(beginning, end)]
        else:
            tmp_dict[(chromosome, strand)].append((beginning, end))

    for key, val in tmp_dict.items():
        inter = interlap.InterLap()
        inter.update(val)
        annotation_intervals_dict[key] = inter

    return annotation_intervals_dict


def filter_annotation_identifiers(
    read_intervals_dict, total_counts_dict, genome_length_dict, mapping_method, config
):
    """
    Filter annotation and return valid CDS identifiers.
    Returns the same format as parse_annotation_cds but filtered.

    In the publication this was done this way, but applying the overlap filter last might be more reasonable to keep more genes.
    """


    # Get filtering parameters
    filtering_methods = config.get("filteringMethods", [])

    # Parse annotation
    annotation_df = pd.read_csv(
        config["annotationFilePath"], sep="\t", comment="#", header=None
    )
    annotation_intervals_dict = create_annotation_intervals_dict(annotation_df)

    # Get filtering parameters
    rpkm_threshold = config.get("rpkmThreshold", 0)
    overlap_distance = config.get("neighboringGenesDistance", 0)
    length_cutoff = config.get("lengthCutoff", 50)  # min length
    positions_out_orf = config.get("positionsOutsideORF", 50)
    positions_in_orf = config.get("positionsInsideORF", 50)

    gene_types = {}
    filtered_identifiers = []
    locus_map = {}

    # First pass: get gene types (same as parse_annotation_cds)
    for row in annotation_df.itertuples(index=False):
        if row[2].lower() in ["gene", "pseudogene"]:
            attributes = row[8]
            gene_id = [
                attr.split("=")[1]
                for attr in attributes.split(";")
                if attr.startswith("ID=")
            ]
            if gene_id:
                gene_types[gene_id[0]] = row[2]

    excluded_genes = {
        "overlap": 0,
        "length": 0,
        "rpkm": 0,
        "type": 0,
        "boundary": 0,
        "pseudogene": 0,
    }
    included_genes = 0

    # Second pass: filter CDS features
    for row in annotation_df.itertuples(index=False):
        chromosome = row[0]
        feature_type = row[2]
        beginning = int(row[3]) - 1
        end = int(row[4]) - 1
        strand = row[6]
        attributes = row[8]

        # Only process CDS
        if feature_type.lower() != "cds":
            continue

        # Check for pseudogenes
        if "Parent=" in attributes:
            parent_id = [
                attr.split("=")[1]
                for attr in attributes.split(";")
                if attr.startswith("Parent=")
            ]
            if parent_id and gene_types.get(parent_id[0]) == "pseudogene":
                excluded_genes["pseudogene"] += 1
                continue

        # Check length divisible by 3
        gene_length = end - beginning + 1
        if gene_length % 3 != 0:
            excluded_genes["type"] += 1
            continue

        # Apply filtering
        if "overlap" in filtering_methods:
            if (chromosome, strand) in annotation_intervals_dict:
                overlapping = list(
                    annotation_intervals_dict[(chromosome, strand)].find(
                        (beginning - overlap_distance, end + overlap_distance)
                    )
                )
                if len(overlapping) > 1:
                    excluded_genes["overlap"] += 1
                    continue

        if "length" in filtering_methods:
            if gene_length < positions_in_orf or gene_length < length_cutoff:
                excluded_genes["length"] += 1
                continue

        if "rpkm" in filtering_methods:
            if read_intervals_dict is None or total_counts_dict is None:
                raise ValueError(
                    "RPKM filtering requires read_intervals_dict and total_counts_dict"
                )

            if (chromosome, strand) in read_intervals_dict:
                gene_read_counts = misc.count_reads(
                    read_intervals_dict,
                    chromosome,
                    strand,
                    beginning,
                    end,
                    mapping_method,
                )
                rpkm = misc.calculate_rpkm(
                    gene_length, gene_read_counts, total_counts_dict[chromosome]
                )
                if rpkm < rpkm_threshold:
                    excluded_genes["rpkm"] += 1
                    continue
            else:
                excluded_genes["rpkm"] += 1
                continue

        # Check boundary conditions
        if genome_length_dict:
            if (
                beginning - positions_out_orf < 0
                or end + positions_out_orf
                > genome_length_dict.get(chromosome, float("inf"))
            ):
                excluded_genes["boundary"] += 1
                continue


        # Gene passed all filters - add to identifiers
        filtered_identifiers.append((chromosome, beginning, end, strand))

        # Build locus_map
        locus_tag = [
            attr.split("=")[1]
            for attr in attributes.split(";")
            if attr.startswith("locus_tag=")
        ]
        if locus_tag:
            locus_map[f"{chromosome}:{beginning+1}-{end+1}:{strand}"] = locus_tag[0]

        included_genes += 1

    # Print filtering statistics
    print("\nAnnotation filtering results:")
    print(f"  Included genes: {included_genes}")
    print("  Excluded genes:")
    for key, count in excluded_genes.items():
        print(f"    - {key}: {count}")

    return filtered_identifiers, locus_map
