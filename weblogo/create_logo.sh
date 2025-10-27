#!/bin/bash

# Sequence extraction and WebLogo generation pipeline
# This script extracts sequences from genome and generates a sequence logo

# Default parameters
XLSX_FILE=""
GFF_FILE=""
FASTA_FILE=""
OUTPUT_FA="sequences.fa"
OUTPUT_SVG="weblogo.svg"
UPSTREAM=30
DOWNSTREAM=10
CUTOFF=0
LOGO_START=-30
LOGO_LENGTH=41

# Function to display usage
usage() {
    cat << EOF
Usage: $0 -x XLSX_FILE -g GFF_FILE -f FASTA_FILE [OPTIONS]

Required arguments:
  -x    Input Excel file with Locus_tag and Read_count columns
  -g    GFF annotation file
  -f    Genome FASTA file

Optional arguments:
  -o    Output FASTA filename (default: sequences.fa)
  -s    Output SVG filename (default: weblogo.svg)
  -u    Nucleotides upstream of start codon (default: 30)
  -d    Nucleotides downstream of start codon (default: 10)
  -c    Minimum read count cutoff (default: 0)
  -h    Display this help message

Examples:
  $0 -x data.xlsx -g annotation.gff -f genome.fa
  $0 -x data.xlsx -g annotation.gff -f genome.fa -u 50 -d 20 -c 100
  $0 -x data.xlsx -g annotation.gff -f genome.fa -o myseqs.fa -s mylogo.svg

EOF
    exit 1
}

# Parse command line arguments
while getopts "x:g:f:o:s:u:d:c:h" opt; do
    case $opt in
        x) XLSX_FILE="$OPTARG" ;;
        g) GFF_FILE="$OPTARG" ;;
        f) FASTA_FILE="$OPTARG" ;;
        o) OUTPUT_FA="$OPTARG" ;;
        s) OUTPUT_SVG="$OPTARG" ;;
        u) UPSTREAM="$OPTARG" ;;
        d) DOWNSTREAM="$OPTARG" ;;
        c) CUTOFF="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check required arguments
if [ -z "$XLSX_FILE" ] || [ -z "$GFF_FILE" ] || [ -z "$FASTA_FILE" ]; then
    echo "Error: Missing required arguments"
    echo ""
    usage
fi

# Check if files exist
if [ ! -f "$XLSX_FILE" ]; then
    echo "Error: Excel file '$XLSX_FILE' not found"
    exit 1
fi

if [ ! -f "$GFF_FILE" ]; then
    echo "Error: GFF file '$GFF_FILE' not found"
    exit 1
fi

if [ ! -f "$FASTA_FILE" ]; then
    echo "Error: FASTA file '$FASTA_FILE' not found"
    exit 1
fi

# Check if extract_sequences.py exists
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXTRACT_SCRIPT="$SCRIPT_DIR/extract_sequences.py"

if [ ! -f "$EXTRACT_SCRIPT" ]; then
    echo "Error: extract_sequences.py not found in $SCRIPT_DIR"
    exit 1
fi

# Calculate logo parameters based on window size
LOGO_LENGTH=$((UPSTREAM + DOWNSTREAM + 1))
LOGO_START=$((-UPSTREAM))

echo "=========================================="
echo "Sequence Extraction and WebLogo Pipeline"
echo "=========================================="
echo ""
echo "Parameters:"
echo "  Excel file:      $XLSX_FILE"
echo "  GFF file:        $GFF_FILE"
echo "  FASTA file:      $FASTA_FILE"
echo "  Output FASTA:    $OUTPUT_FA"
echo "  Output SVG:      $OUTPUT_SVG"
echo "  Upstream:        $UPSTREAM nt"
echo "  Downstream:      $DOWNSTREAM nt"
echo "  Read cutoff:     $CUTOFF"
echo "  Logo start pos:  $LOGO_START"
echo "  Logo length:     $LOGO_LENGTH"
echo ""

# Step 1: Extract sequences
echo "Step 1: Extracting sequences..."
echo "----------------------------------------"
python "$EXTRACT_SCRIPT" \
    -x "$XLSX_FILE" \
    -g "$GFF_FILE" \
    -f "$FASTA_FILE" \
    -o "$OUTPUT_FA" \
    -u "$UPSTREAM" \
    -d "$DOWNSTREAM" \
    -c "$CUTOFF"

if [ $? -ne 0 ]; then
    echo ""
    echo "Error: Sequence extraction failed"
    exit 1
fi

echo ""

# Step 2: Generate WebLogo
echo "Step 2: Generating WebLogo..."
echo "----------------------------------------"

# Check if weblogo is installed
if ! command -v weblogo &> /dev/null; then
    echo "Error: weblogo command not found"
    echo "Please install weblogo: pip install weblogo"
    exit 1
fi

weblogo -f "$OUTPUT_FA" \
    -F svg \
    -o "$OUTPUT_SVG" \
    -i "$LOGO_START" \
    -n "$LOGO_LENGTH" \
    -c classic \
    -s large

if [ $? -ne 0 ]; then
    echo ""
    echo "Error: WebLogo generation failed"
    exit 1
fi

echo ""
echo "=========================================="
echo "Pipeline completed successfully!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  Sequences: $OUTPUT_FA"
echo "  WebLogo:   $OUTPUT_SVG"
echo ""