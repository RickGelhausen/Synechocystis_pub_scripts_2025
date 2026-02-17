#!/bin/bash
# Common parameters
GFF="annotation_confident.gff"
GENOME="genome.fa"
INPUT_DIR="input"
OUTPUT_DIR="output"

# Define samples and cutoffs
SAMPLES=(
    "TIS-WT-1_4_32_f"
    "TIS-WT-1_2_32_f"
    "TIS-WT-1_19_32_f"
    "TIS-WT-2_4_32_f"
    "TIS-WT-2_2_32_f"
    "TIS-WT-2_19_32_f"
    "TIS-WT_4_32_f_merged"
    "TIS-WT_2_32_f_merged"
    "TIS-WT_19_32_f_merged"
)
CUTOFFS=(0 10)

# Function to determine sample type
get_sample_type() {
    local sample=$1
    if [[ $sample == TIS-WT-1_* ]]; then
        echo "TIS-WT-1"
    elif [[ $sample == TIS-WT-2_* ]]; then
        echo "TIS-WT-2"
    elif [[ $sample == TIS-WT_* ]]; then
        echo "TIS-WT"
    fi
}

# Process all samples with cutoffs 0 and 10
for sample in "${SAMPLES[@]}"; do
    sample_type=$(get_sample_type "$sample")

    for cutoff in "${CUTOFFS[@]}"; do
        # Create organized directory structure
        SAMPLE_DIR="$OUTPUT_DIR/${sample_type}/c${cutoff}"
        FASTA_DIR="$SAMPLE_DIR/fasta"
        PLOTS_DIR="$SAMPLE_DIR/plots"
        STARTCODONS_DIR="$SAMPLE_DIR/startcodons"

        mkdir -p "$FASTA_DIR" "$PLOTS_DIR" "$STARTCODONS_DIR"

        # Run create_logo.sh with organized output paths
        bash create_logo.sh -x "$INPUT_DIR/${sample}.xlsx" -g "$GFF" -f "$GENOME" \
            -o "$FASTA_DIR/${sample}.fa" \
            -s "$PLOTS_DIR/${sample}.svg" \
            -c "$cutoff"

    done
done