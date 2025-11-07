# Ribosome Profiling Analysis Tools

Code repository for the publication *Combined translation initiation and termination site ribosome profiling precisely defines novel small and alternative reading frames in a model cyanobacterium* ([DOI](TODO)).

## Overview

This repository provides three analysis tools:

- **`metagene_analysis`** — Scripts for performing ribosome density analysis
- **`locus_extraction`** — Scripts to extract locus information for specific positions in density plots
- **`weblogo`** — Scripts for generating weblogos from genome sequences

## Installation

### Using uv (Recommended)

We recommend [`uv`](https://github.com/astral-sh/uv) for dependency management:

```bash
uv sync
source .venv/bin/activate
```

Or run commands directly:
```bash
uv run <command>
```

### Using mamba

Alternatively, use [`mamba`](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html):

```bash
mamba create -f metagene_profiling.yml
mamba activate metagene_profiling
```

### Requirements

- Python ≥ 3.14.0
- pysam ≥ 0.23.3
- pandas ≥ 2.3.3
- numpy ≥ 2.3.4
- biopython ≥ 1.8.6
- interlap ≥ 0.2.7
- openpyxl ≥ 3.1.5
- weblogo ≥ 3.7.9
- xlsxwriter ≥ 3.2.9
- xlrd ≥ 2.0.2
- python-kaleido ≥ 1.2.0
- pytest ≥ 8.4.2
- pytest-mock ≥ 3.15.1
- pyyaml ≥ 6.0.3
- plotly ≥ 6.3.1

> Other versions may work but were not tested.

## Usage

### 1. Metagene Analysis

Generates per-position read-density plots.

**Process:**
1. Read input annotation files and extract CDS features
2. Apply user-defined filters (overlap, length, minimum coverage/RPKM)
3. Define windows based on start and stop codon positions
4. Count aligned reads per read_length, position, and locus
5. Collect read counts in dataframes and output as Excel tables
6. Generate plots from dataframes

**Run the analysis:**

```bash
# Create config file from template
cp config_sample.yml my_config.yml

# Run analysis
python3 metagene_profiling.py -c my_config.yml
```

**Fast replotting** (using pre-processed loci.json files):

```bash
python3 metagene_plotting.py -c config.yaml -i path_to_profiling_json_dir -o path_to_plot_folder
```

> ⏱️ Runtime can be several minutes depending on the number and size of BAM files.

#### Configuration Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `alignmentDirPath` | Path to folder containing BAM alignment files | — |
| `annotationFilePath` | Path to GFF annotation file | — |
| `genomeFilePath` | Path to genome reference file | — |
| `outputDirPath` | Path to output folder | — |
| `positionsOutsideORF` | Nucleotides upstream (downstream for stop) of each ORF | 50 |
| `positionsInsideORF` | Nucleotides within each ORF | 50 |
| `filteringMethods` | Filter methods: `overlap`, `length`, `rpkm` | — |
| `neighboringGenesDistance` | Distance to consider for overlapping genes (0 = direct overlaps only) | 50 |
| `rpkmThreshold` | Minimum RPKM threshold | 10.0 |
| `lengthCutoff` | Minimum ORF length (ORFs smaller than `positionsInORF` auto-filtered) | 50 |
| `mappingMethods` | Read mapping methods: `fiveprime`, `threeprime`, `centered`, `global` | — |
| `readLengths` | Read lengths to consider (e.g., "16-18,23-25") | — |
| `normalizationMethods` | Normalization: `raw`, `cpm`, `window` | — |
| `outputFormats` | Output formats: `svg`, `pdf`, `png`, `jpg`, `interactive` | — |
| `includePlotlyJS` | Plotly.js inclusion: `integrated` (~3.7MB/file), `online` (requires internet), `local` | — |
| `colorList` | Custom colors (auto-generated if empty) | — |

### 2. Locus Extraction

Extract loci contributing to specific peaks in density plots.

> **Why this matters:** Peaks can result from artifacts (size selection, cDNA synthesis, sequencing) rather than genuine biological signals. This tool helps investigate peak composition.

**Usage:**

```bash
python3 locus_extraction/extract_locus_data.py \
    -i=/path/to/metagene/files/readcounts_start_threeprime.json \
    -s <sample_name> \
    -r <read_length> \
    -p <position> \
    -o output.xlsx \
    --sort-by count
```

Requires the intermediate `.json` files from `metagene_profiling.py`.

### 3. Weblogo

Generate weblogos from genome sequences after extracting regions with `extract_locus_data.py`.

**Helper script:** `merge_xlsx.py` merges multiple Excel sheets (useful for combining replicates and keeping best read counts per locus).

**Usage:**

```bash
bash weblogo/create_logo.sh \
    -x <input>.xlsx \
    -f genome.fa \
    -g annotation.gff \
    -o output_sequences.fa \
    -s output_logo.svg
```

Output formats: PNG (default), SVG (requires `pdf2svg`), or PDF (use `-t` parameter).

## Demo

Synthetic demo data is provided in `demo_data/` (no biological significance).

### Run metagene analysis:

```bash
python3 metagene_analysis/metagene_profiling.py -c demo_data/config.yaml
# or depending on your setup (this applies to all the following commands)
uv run python3 metagene_analysis/metagene_profiling.py -c demo_data/config.yaml
```

Creates plots with high peaks at +12nt for 5' mapping.

### Extract loci from peaks:

```bash
python3 locus_extraction/extract_locus_data.py \
    -i demo_data/metagene_output/readcounts_start_fiveprime.json \
    -g demo_data/annotation.gff \
    -s TIS-test-2 \
    -r 32 \
    -p 12 \
    -o demo_data/extraction/demo_extraction.xlsx \
    --sort-by count \
    --gff-output demo_data/extraction/extraction_annotated.gff
```

### Generate weblogo:

```bash
bash weblogo/create_logo.sh \
    -x demo_data/extraction/demo_extraction.xlsx \
    -f demo_data/genome.fa \
    -g demo_data/annotation.gff \
    -o demo_data/weblogo/output_sequences.fa \
    -s demo_data/weblogo/output_logo.png
```


## Tests

Unittests for the metagene_analysis code can be run using `pytest`.

```
pytest metagene_analysis
```

## License

GPL-3.0 license


## Citation

If you use this code, please cite:

```
[Add citation here once DOI is available]
```