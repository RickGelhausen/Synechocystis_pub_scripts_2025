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

> **⚠️ Platform Note:** This tool was developed and tested on a Linux system (Ubuntu 22.04.5 LTS, Ubuntu 20.04 LTS). It should not contain Linux-specific commands, but it was never tested on Windows or macOS.

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
## Reproduction Instructions

This section provides a detailed, step-by-step guide for reproducing all data used in the study. While some steps could be simplified, we describe the exact process as it was originally carried out.

### Prerequisites

Several steps in this guide rely on alignment files (`.bam`) as input. This includes the metagene profiling scripts described in this repository, as well as peak calling using ORFBounder. Generating these alignment files is therefore the first step.

> **Note:** All external tools used in this guide (HRIBO, ORFBounder, SRA Tools, etc.) have their own detailed documentation and usage instructions on their respective GitHub repositories. This guide focuses on how we applied these tools to reproduce our specific results, not on how to install or configure them in general.

### Step 1: Generate Alignment Files, Coverage Files, and Ribo-Seq-Based Predictions

Quality control, alignment, and all other analyses described in this subsection were performed using the [HRIBO pipeline](https://github.com/RickGelhausen/HRIBO).

#### 1.1 Download Raw Data from SRA

All raw data are available via [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246336) and can be downloaded using [SRA Tools](https://github.com/ncbi/sra-tools).

Use the `download.sh` script to download the data in the correct format for further processing.

> **Note:** Due to the timeline of the original analysis, TIS and TTS data were processed in two separate HRIBO runs.

#### 1.2 Download HRIBO

Download the HRIBO pipeline into each subfolder (`TIS_run` and `TTS_run`):

```bash
wget https://github.com/RickGelhausen/HRIBO/archive/refs/tags/1.8.1.tar.gz
tar -xzf 1.8.1.tar.gz; mv HRIBO-1.8.1 HRIBO; rm 1.8.1.tar.gz;
```

#### 1.3 Prepare Sample and Config Files

The sample sheets and config files used in the analysis are available in the [`reproduction_data` folder](https://github.com/RickGelhausen/Synechocystis_pub_scripts_2025/tree/main/reproduction_data/HRIBO_data) for both TIS and TTS. Copy them into the respective HRIBO folders you just created.

#### 1.4 Run HRIBO

> **Note:** This workflow involves read mapping. While the bacterial genome is relatively small, you may still encounter memory issues on machines with limited resources. A minimum of **16 GB RAM** is recommended.

```bash
snakemake --use-conda --use-singularity --singularity-args " -c " -s HRIBO/Snakefile --configfile config.yaml -j 20 --latency-wait 50
```

> **Note:** The `--singularity-args " -c "` flag was needed for the Snakemake version used during our analysis and may not be required for newer versions. Consult the [Snakemake documentation](https://snakemake.readthedocs.io/) if you encounter compatibility issues.

#### 1.5 Inspect Results

This generates the alignment files (`.bam`) required for downstream analyses, as well as an overview result file (`.xlsx`) with quality control metrics and prediction results.

---

### Step 2: Metagene Profiling and Weblogos

Metagene profiling was conducted using the scripts available in this repository. It requires the `.bam` files generated in Step 1.

#### 2.1 Prepare Input

The metagene profiling script requires a config file with the appropriate setup. The config files used for the study are available in the [`reproduction_data` folder](https://github.com/RickGelhausen/Synechocystis_pub_scripts_2025/tree/main/reproduction_data/metagene_data).

#### 2.2 Adjust the Config File

Update the config file to match your local file system. Make sure the paths to the annotation, genome, and `.bam` files are set correctly.

#### 2.3 Run the Metagene Profiling Script

```bash
uv run python3 metagene_profiling.py -c config_TIS.yaml
uv run python3 metagene_profiling.py -c config_TTS.yaml
```

This generates output tables and figures. Additional publication-quality figures were generated from these output tables.

#### 2.4 Extract Metagene Data for Weblogos

The previous step produces `.json` files containing detailed per-gene read count information. To generate weblogos, we identified genes contributing to positions of interest in the metagene plots and extracted them using the locus extraction script.

> **Note:** The extraction scripts are designed to be run from within their respective directories, as each step constituted an independent analysis.

---

### Step 3: ORFBounder

Running ORFBounder requires several input files: annotation, genome, alignment files, offset, read length, and mapped read counts.

- **Annotation, genome, and alignment files** can be taken directly from the HRIBO output (Step 1).
- **All other configuration files** (including the `mapped_reads` file, which is also available in the HRIBO `readcounting` folder) are provided in the [`reproduction_data` folder](https://github.com/RickGelhausen/Synechocystis_pub_scripts_2025/tree/main/reproduction_data/ORFBounder_data).

#### 3.1 Configure ORFBounder

Set the paths to all input files in `config_experiments.tsv`. This file allows you to specify and run multiple ORFBounder experiments sequentially.

#### 3.2 Run ORFBounder

```bash
uv run call_ORFBounder.py -c config_experiments.tsv -r orfbounder_output_folder
```

This creates an output folder for each experiment, containing Excel sheets and GFF files.

---

### Step 4: Filtering and Manual Inspection

All remaining analysis steps were performed manually using Excel/LibreOffice and a genome browser. Details on this process are described in the manuscript.

## License

GPL-3.0 license


## Citation

If you use this code, please cite:

```
[Add citation here once DOI is available]
```
