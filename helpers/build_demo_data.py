#!/usr/bin/env python3
"""
Generate minimal toy data for metagene profiling demonstration.
Creates: genome.fa, annotation.gff, sample1.bam, sample2.bam
"""

from pathlib import Path

import pysam


def create_genome_file(output_path):
    """Create a minimal genome FASTA file with proper start/stop codons"""
    print("Creating genome.fa...")

    # Start with mostly A's
    sequence = list("A" * 5000)

    # Gene 1: positions 500-799 (300 bp, divisible by 3)
    # Start codon ATG at position 500
    sequence[500:503] = list("ATG")
    # Stop codon TAA at position 797-799
    sequence[797:800] = list("TAA")

    # Gene 2: positions 1500-1898 (399 bp, divisible by 3)
    # Start codon ATG at position 1500
    sequence[1500:1503] = list("ATG")
    # Stop codon TAG at position 1896-1898
    sequence[1896:1899] = list("TAG")

    # Gene 3: positions 3000-3499 (500 bp, not divisible by 3, let's fix to 501)
    # Actually, let's make it 498 bp (divisible by 3)
    # Start codon ATG at position 3000
    sequence[3000:3003] = list("ATG")
    # Stop codon TGA at position 3495-3497
    sequence[3495:3498] = list("TGA")

    sequence_str = "".join(sequence)

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(">chr1 chromosome 1\n")
        # Write in 80-character lines
        for i in range(0, len(sequence_str), 80):
            f.write(sequence_str[i:i+80] + "\n")

    print(f"  Created: {output_path}")
    print("    Gene 1: start codon at 500 (ATG), stop codon at 797 (TAA), length=300bp")
    print("    Gene 2: start codon at 1500 (ATG), stop codon at 1896 (TAG), length=399bp")
    print("    Gene 3: start codon at 3000 (ATG), stop codon at 3495 (TGA), length=498bp")


def create_annotation_file(output_path):
    """Create a minimal GFF3 annotation file with 3 genes"""
    print("Creating annotation.gff...")

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("##gff-version 3\n")

        # Gene 1: positions 500-799 (300 bp = 100 codons)
        f.write("chr1\t.\tgene\t500\t799\t.\t+\t.\tID=gene1;Name=gene1\n")
        f.write("chr1\t.\tCDS\t500\t799\t.\t+\t0\tID=cds1;Parent=gene1\n")

        # Gene 2: positions 1500-1898 (399 bp = 133 codons)
        f.write("chr1\t.\tgene\t1500\t1898\t.\t+\t.\tID=gene2;Name=gene2\n")
        f.write("chr1\t.\tCDS\t1500\t1898\t.\t+\t0\tID=cds2;Parent=gene2\n")

        # Gene 3: positions 3000-3497 (498 bp = 166 codons)
        f.write("chr1\t.\tgene\t3000\t3497\t.\t+\t.\tID=gene3;Name=gene3\n")
        f.write("chr1\t.\tCDS\t3000\t3497\t.\t+\t0\tID=cds3;Parent=gene3\n")

    print(f"  Created: {output_path}")


def create_bam_file(output_path, sample_name, read_counts):
    """
    Create a BAM file with reads concentrated around start codons.

    Args:
        output_path: Path to output BAM file
        sample_name: Name of the sample
        read_counts: Dict mapping read_length to count for peak positions
                     e.g., {28: 20, 30: 25, 32: 30}
    """
    print(f"Creating {output_path.name}...")

    # BAM header
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 5000, 'SN': 'chr1'}],
        'RG': [{'ID': sample_name, 'SM': sample_name}]
    }

    # Gene start positions (0-based, start of ATG codon)
    gene_starts = [499, 1499, 2999]  # 0-based coordinates (GFF is 1-based)

    read_id = 0

    with pysam.AlignmentFile(output_path, "wb", header=header) as outf:
        for gene_start in gene_starts:
            # Create reads at various positions around each gene start
            # Positions relative to start codon: -50 to +50
            positions = list(range(-50, 51, 5))

            for rel_pos in positions:
                abs_pos = gene_start + rel_pos

                if abs_pos < 0:
                    continue

                # Create reads for different read lengths
                for read_length, peak_count in read_counts.items():
                    # Calculate read count based on distance from start
                    if -15 <= rel_pos <= 15:
                        # Strong peak at start codon
                        count = peak_count
                    elif -30 <= rel_pos <= 30:
                        # Medium coverage nearby
                        count = peak_count // 2
                    else:
                        # Low background
                        count = max(1, peak_count // 7)

                    for _ in range(count):
                        a = pysam.AlignedSegment()
                        a.query_name = f"{sample_name}_read_{read_id}"
                        a.query_sequence = "A" * read_length
                        a.flag = 0
                        a.reference_id = 0
                        a.reference_start = abs_pos
                        a.mapping_quality = 60
                        a.cigar = [(0, read_length)]  # M (match)
                        a.query_qualities = pysam.qualitystring_to_array("I" * read_length)

                        # Add required tags
                        a.set_tag('NH', 1)   # Unique mapping
                        a.set_tag('HI', 1)   # Hit index
                        a.set_tag('AS', read_length)  # Alignment score
                        a.set_tag('NM', 0)   # Edit distance
                        a.set_tag('RG', sample_name)  # Read group

                        outf.write(a)
                        read_id += 1

    # Index the BAM file
    pysam.index(str(output_path))

    print(f"  Created: {output_path} ({read_id} reads)")
    print(f"  Created: {output_path}.bai")

def create_config_file(output_path):
    """Create a sample configuration file (YAML)"""
    print("Creating config.yaml...")

    config_content = f"""\
    ## Configuration file for metagene profiling

# Input files
alignmentDirPath: "{output_path}/bam_files"
annotationFilePath: "{output_path}/annotation.gff"
genomeFilePath: "{output_path}/genome.fa"

# Output folder
outputDirPath: "{output_path}/metagene_output/"

### Metagene settings
# Window on which the metagene profiling is performed
# Number of nucleotides upstream of each annotated ORF that will be considered
positionsOutsideORF: 50
# Number of nucleotides in each annotated ORF that will be considered
positionsInsideORF: 100

# methods to filter the input annotation before metagene profiling is performed (overlap, length and/or rpkm)
filteringMethods: []

# Distance around each annotated genes that will be considered for removing overlapping genes
# neighboringGenesDistance 0 only removes directly overlapping genes.
neighboringGenesDistance: 50

# Filter ORFs based on RPKM threshold
rpkmThreshold: 10.0

# Minimum length of considered ORFs. ORFs smaller than positionsInORF are filtered automatically.
lengthCutoff: 50

# Mapping methods to be considered (fiveprime, threeprime, centered, global)
mappingMethods: ["threeprime"]

# Read lengths to be considered. Comma seperated string allows intervals denoted by "-" symbol (e.g. "22,23,27,34-35")
readLengths: "28-32"

# Normalization method (raw, cpm, window)
# raw: no normalization, raw read counts
# cpm: counts per million normalization (CPM)
# window: normalized data by total counts and considered window length.
normalizationMethods: ["raw"]

### Plotting
# Plot output types. Supported (svg, pdf, png, jpg, interactive)
outputFormats: ["interactive"]

# Interactive html files require more space as they require plotly.js integrated, online access or a local copy of plotly.js
# integrated: integrates plotly.js into each plot (+-3.7mb per file)
# online: references an online version of plotly.js (files can only be used with active internet connection)
# local: references a local plotly.js script (requires a local copy of plotly.js to open the plots)
includePlotlyJS: "integrated"

# Custom list of colors to be used for the plots. Must be atleast as many as the number of read lengths considered.
# If empty we try to match the number of read lengths trying to prioritize color-blind friendly colors.
colorList: []
    """

    with open(output_path / "config.yaml", 'w', encoding='utf-8') as f:
        f.write(config_content)

    print(f"  Created: {output_path / "config.yaml"}")


def main():
    """Generate all demo data files"""
    print("\n" + "="*60)
    print("Generating Metagene Profiling Demo Data")
    print("="*60 + "\n")

    # Create output directory
    output_dir = Path("demo_data")
    output_dir.mkdir(exist_ok=True)

    (output_dir / "bam_files").mkdir(exist_ok=True)

    # Generate files
    create_config_file(output_dir)
    create_genome_file(output_dir / "genome.fa")
    create_annotation_file(output_dir / "annotation.gff")
    create_bam_file(output_dir / "bam_files" / "TIS-test-1.bam", "TIS-test-1", {28: 20, 30: 25, 32: 30})
    create_bam_file(output_dir / "bam_files" / "TIS-test-2.bam", "TIS-test-2", {28: 18, 30: 23, 32: 28})




if __name__ == "__main__":
    main()