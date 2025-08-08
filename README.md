# Single-Cell RNA-seq Perturbation Pipeline

A Snakemake pipeline for processing Parse Biosciences single-cell RNA-seq data with CRISPR perturbations, including kallisto-bustools alignment, cell calling, quality control, and analysis-ready file generation.

## Overview

This pipeline processes Parse Bio scRNA-seq data with guide RNA (gRNA) readouts to analyze CRISPR perturbation experiments. It handles:
- FASTQ alignment using kallisto-bustools
- Multiple cell calling methods (EmptyDrops, knee/inflection point, etc.)
- Barcode recovery from undetermined reads
- Quality control and visualization
- Generation of analysis-ready files for downstream analysis

The documentation below uses the K562 CRISPR screen configuration (`config.yann_k562.yaml`) as a concrete example.

## Table of Contents
1. [Quick Start](#quick-start)
2. [Required Input Files](#required-input-files)
3. [Configuration File Structure](#configuration-file-structure)
4. [File Format Specifications](#file-format-specifications)
5. [Running the Pipeline](#running-the-pipeline)
6. [Output Structure](#output-structure)
7. [Example: K562 CRISPR Screen](#example-k562-crispr-screen)

---

## Quick Start

```bash
# 1. Clone the repository
git clone [repository_url]
cd perturb_pipeline

# 2. Activate the existing kb conda environment
conda activate /oak/stanford/groups/engreitz/Users/tonyzeng/miniconda3/envs/kb
# Note: The kb environment is pre-installed with all necessary dependencies
# If you need to set up your own environment, contact the maintainer for the environment.yml file

# 3. Copy and modify the example config
cp gw_analysis/config.yann_k562.yaml gw_analysis/config.my_experiment.yaml
# Edit gw_analysis/config.my_experiment.yaml with your paths and parameters

# 4. Prepare your input files (see Required Input Files section)

# 5. Run the pipeline
cd gw_analysis
CONFIG=config.my_experiment.yaml ./submit.sh
```

---

## Required Input Files

Before running the pipeline, you need to prepare the following files:

### 1. Sample Information Excel File (`sample_info.xlsx`)
**Purpose:** Defines all samples, their properties, and relationships  
**Location:** Referenced in config as `sample_info_file`

**Required columns:**
- `pool`: Pool identifier (e.g., "1k_pilot", "pool1", "pool2")
- `sample`: Sample name within the pool (e.g., "sample1_gex", "sample1_grna")
- `sample_id`: Unique identifier combining pool and sample (e.g., "1k_pilot:sample1_gex")
- `sample_type`: Either "gex" (gene expression) or "guide" (guide RNA)
- `expected_cells`: Expected number of cells for this library (e.g., 5000)
- `min_umi_threshold`: Minimum UMI count threshold for cell calling (e.g., 100)
- `fastq_dir`: Directory containing FASTQ files for this sample
- `sample_to_well_mapping`: Name of the plate mapping to use (e.g., "plate1")
- `i7_barcode`: i7 index sequence (8bp, e.g., "AAGGCTAT")
- `i5_barcode`: i5 index sequence (8bp, e.g., "AAACATCG")
- `notes`: Optional notes about the sample
- `read1_length`: Length of Read 1 (typically 150)
- `read2_length`: Length of Read 2 (typically 150)
- `paired_guide_pool`: Pool of the paired guide library (for GEX samples)
- `paired_guide_sample`: Sample name of paired guide library (for GEX samples)
- `paired_guide_sample_id`: Full ID of paired guide library (for GEX samples)

**Example (from K562 screen):**
```
pool        sample          sample_id                sample_type  expected_cells  fastq_dir
1k_pilot    sample1_gex     1k_pilot:sample1_gex     gex         5000           /data/fastqs/
1k_pilot    sample1_grna    1k_pilot:sample1_grna    guide       5000           /data/fastqs/
```

💡 **Tip:** See `config.yann_k562.yaml` for a complete working example configuration.

### 2. Plate Mapping Excel File (`plate_map.xlsx`)
**Purpose:** Maps well positions to biological samples  
**Location:** Referenced in config as `plate_maps_file`

**Structure:**
- Each sheet represents a different plate (sheet name matches `sample_to_well_mapping` in sample_info)
- Required columns per sheet:
  - `Well Position`: Well identifier (A1-H12 for 96-well plate)
  - `Sample`: Biological sample name
  - `Rep`: Replicate identifier

**Example (sheet "plate1"):**
```
Well Position    Sample    Rep
A1              rep1      rep1
A2              rep4      rep4
A3              rep1      rep1
...
```

### 3. Guide Library File (`guides.txt`)
**Purpose:** Defines all guide RNA sequences for kallisto indexing  
**Location:** Referenced in config as `guide_file`

**Format:** Tab-separated text file with columns:
- Column 1: Guide ID
- Column 2: Guide sequence (typically 20bp)

**Example:**
```
AAVS1_guide1    ACCCCACAGTGGGGCCACTA
AAVS1_guide2    GGGGCCACTAGGGACAGGAT
TP53_guide1     GAAATTTGCGTGTGGAGTAT
```

### 4. Guide Reference File (`guides_qc_reference.txt`)
**Purpose:** Maps guide IDs to target genes for QC analysis  
**Location:** Referenced in config as `guide_reference`

**Format:** Tab-separated with columns:
- `guide_id`: Must match IDs in guides.txt
- `target_gene`: Gene symbol
- `guide_type`: "targeting" or "non-targeting"
- `sequence`: Guide sequence (optional, for validation)

**Example:**
```
guide_id        target_gene    guide_type      sequence
AAVS1_guide1    AAVS1         targeting       ACCCCACAGTGGGGCCACTA
NTC_guide1      None          non-targeting   GCTAGCTAGCTAGCTAGCTA
```

### 5. Parse Bio Barcode Files
**Purpose:** Define the barcode sequences for demultiplexing  
**Locations:**
- `replace_file`: Replace barcode sequences (96 sequences for 96-well)
- `barcodes_file`: Cell barcode sequences (96 sequences for 96-well)

**Format:** Plain text, one barcode per line (must have exactly 96 lines for 96-well kit)

### 6. Reference Files
**Purpose:** Gene annotations and genome indices  
**Required files:**
- `ribosomal_genes`: List of ribosomal gene names (one per line)
- `cell_cycle_genes`: List of cell cycle genes (Regev lab format)
- `gene_database`: Gene annotation database (gencode format)
- Kallisto index files in `reference_base/nascent_genome/`

---

## Configuration File Structure

The config file has two main sections:

### Run-Specific Configuration
These settings change for each experiment:
```yaml
# Analysis identification
analysis_name: "yann_k562_analysis"  # Unique name for this analysis

# Sample information files
sample_info_file: "/path/to/sample_info.xlsx"
plate_maps_file: "/path/to/plate_map.xlsx"

# Guide library file
guide_file: "/path/to/guides.txt"

# Output directories
output_paths:
  results_base: "/path/to/results"
  logs_base: "/path/to/logs"

# Sample selection
sample_ids: []  # Empty = process all samples

# Analysis parameters
analysis:
  combinations:
    - ["main", "raw"]  # Source and processing type
  threads: 24
  strand: "forward"  # or "unstranded"
```

### Stable Infrastructure
These settings typically remain constant:
```yaml
# Parse Bio kit configuration
input_paths:
  replace_file: "/path/to/replace.96.txt"
  barcodes_file: "/path/to/barcodes.96.txt"
  guide_reference: "/path/to/guides_qc_reference.txt"
  
# Cell calling parameters
cell_calling:
  methods_to_run:
    - Expected_Cells
    - BarcodeRanks_Knee
  default_method: "Expected_Cells"
  
# Resource allocation
resources:
  alignment:
    mem_mb: 128000  # 128GB
    threads: 24
```

---

## File Format Specifications

### FASTQ File Naming Convention
The pipeline expects FASTQ files with specific naming patterns:

**Standard Illumina format:**
```
{sample_name}_S{number}_L{lane}_R{read}_001.fastq.gz
```

**Example:**
```
sample1_gex_S1_L001_R1_001.fastq.gz
sample1_gex_S1_L001_R2_001.fastq.gz
```

**For Parse Bio data:**
- R1: Contains cell barcode and UMI (150bp)
- R2: Contains transcript/guide sequence (150bp)

### Directory Structure
Organize your data as follows:
```
perturb_pipeline/
├── README.md
├── gw_analysis/
│   ├── Snakefile
│   ├── submit.sh
│   ├── config.yann_k562.yaml
│   ├── scripts/
│   └── references_yann_k562/    # Experiment-specific references
│       ├── sample_info.k562_yann.xlsx
│       ├── plate_map.xlsx
│       ├── k562_yann_guides.txt
│       └── k562_yann_guides_qc_reference.txt
└── results_yann_k562/            # Output directory (created by pipeline)
    └── logs_yann_k562/           # Log files (created by pipeline)
```

---

## Running the Pipeline

### 1. Activate Conda Environment
```bash
conda activate /oak/stanford/groups/engreitz/Users/tonyzeng/miniconda3/envs/kb
```

### 2. Navigate to Pipeline Directory
```bash
cd gw_analysis
```

### 3. Run the Pipeline

The pipeline is run using the `submit.sh` script, which handles SLURM submission and resource allocation. You must set the CONFIG environment variable to specify your configuration file.

**Basic usage:**
```bash
# Run the full pipeline
CONFIG=config.yann_k562.yaml ./submit.sh

# Dry run to check what will be executed
CONFIG=config.yann_k562.yaml ./submit.sh --dry-run

# Run with specific number of jobs
CONFIG=config.yann_k562.yaml ./submit.sh --jobs 50

# Force rerun based on file modification times
CONFIG=config.yann_k562.yaml ./submit.sh --rerun-triggers mtime

# Run specific targets only
CONFIG=config.yann_k562.yaml ./submit.sh kallisto_gex

# Run with custom Snakemake arguments
CONFIG=config.yann_k562.yaml ./submit.sh --dry-run --rerun-triggers mtime
```

**Note:** The `submit.sh` script:
- Automatically submits jobs to SLURM with appropriate resource requests
- Sets default resources (4 threads, 32GB RAM, 48-hour runtime)
- Handles retries and keeps incomplete files
- Runs up to 100 jobs in parallel by default

---

## Output Structure

The pipeline creates the following output structure:

```
results_base/
├── {pool}/
│   ├── {sample_id}/
│   │   ├── counts.txt                    # Read counts
│   │   ├── kb_all/                       # Kallisto-bustools output
│   │   │   └── counts_unfiltered/
│   │   │       └── adata.h5ad           # Expression matrix
│   │   ├── cell_calling/                 # Cell calling results
│   │   │   ├── results.tsv              # Summary of all methods
│   │   │   └── {sample}_{method}_cell_barcodes.txt
│   │   └── saturation/                   # Saturation analysis
│   │       └── saturation_curves.png
│   └── {pool}:Undetermined/              # Undetermined reads
├── qc_report/
│   ├── data/                             # QC metrics data
│   └── plots/                            # QC visualization
└── analysis_ready/                       # Final processed data
    ├── gex_all.h5ad                      # Combined GEX data
    └── guide_all.h5ad                    # Combined guide data
```

---

## Troubleshooting

### Common Issues

1. **Missing FASTQ files**
   - Check that `fastq_dir` in sample_info.xlsx points to correct location
   - Verify FASTQ naming matches expected pattern

2. **Cell calling fails**
   - Ensure `expected_cells` in sample_info.xlsx is reasonable
   - Check that kallisto alignment completed successfully

3. **Guide mapping issues**
   - Verify guide IDs match between guides.txt and guides_qc_reference.txt
   - Check guide sequences are correct length (typically 20bp)

4. **Memory errors**
   - Adjust resource allocations in config:
     ```yaml
     resources:
       alignment:
         mem_mb: 256000  # Increase memory
     ```

5. **Barcode issues**
   - Ensure using correct Parse Bio kit files (96-well vs 48-well)
   - Verify barcode files have correct number of sequences

### Validation Steps

Before running the full pipeline:

1. **Check sample info:**
   ```python
   import pandas as pd
   df = pd.read_excel('sample_info.xlsx')
   # Check for required columns
   required = ['pool', 'sample', 'sample_id', 'sample_type', 'fastq_dir']
   missing = [c for c in required if c not in df.columns]
   print(f"Missing columns: {missing}")
   ```

2. **Verify FASTQ files exist:**
   ```python
   import os
   for _, row in df.iterrows():
       fastq_dir = row['fastq_dir']
       sample = row['sample']
       r1 = f"{fastq_dir}/{sample}_R1.fastq.gz"
       if not os.path.exists(r1):
           print(f"Missing: {r1}")
   ```

3. **Check guide consistency:**
   ```bash
   # Count guides in each file
   wc -l guides.txt
   cut -f1 guides_qc_reference.txt | sort -u | wc -l
   ```

---

## Example: K562 CRISPR Screen

The `config.yann_k562.yaml` file provides a complete working example for a K562 CRISPR perturbation screen. This configuration:

### Dataset Details
- **Cell line:** K562 (human chronic myelogenous leukemia)
- **Perturbation:** CRISPR interference (CRISPRi)
- **Scale:** ~1,000 cells per sample (pilot experiment)
- **Libraries:** 2 GEX + 2 guide RNA libraries

### Key Configuration Settings
```yaml
# Analysis name identifies this experiment
analysis_name: "yann_k562_analysis"

# Cell calling uses multiple methods
cell_calling:
  methods_to_run:
    - Expected_Cells      # Top N cells by UMI
    - BarcodeRanks_Knee   # Knee point detection
  default_method: "Expected_Cells"

# Processing configuration
analysis:
  combinations:
    - ["main", "raw"]     # Process main reads only
  threads: 24
  strand: "forward"       # Parse Bio is forward-stranded
```

### Running the K562 Example
```bash
# Navigate to the pipeline directory
cd gw_analysis

# Use the K562 config directly
CONFIG=config.yann_k562.yaml ./submit.sh

# Or copy it as a template for your own experiment
cp config.yann_k562.yaml config.my_screen.yaml
# Edit paths and parameters for your experiment
CONFIG=config.my_screen.yaml ./submit.sh
```

### Expected Outputs for K562
The pipeline will generate:
- **Alignment results:** ~5,000 cells per library after filtering
- **Cell calling comparison:** Multiple methods compared in QC reports
- **Guide assignment:** Cells assigned to perturbations based on guide UMI counts
- **QC metrics:** Comprehensive quality control visualizations
- **Analysis-ready files:** Combined h5ad files for downstream analysis

---

## Troubleshooting

### Common Issues

1. **Missing FASTQ files**
   - Check that `fastq_dir` in sample_info.xlsx points to correct location
   - Verify FASTQ naming matches expected pattern
   - For Parse Bio data, ensure you have both R1 and R2 files

2. **Cell calling fails**
   - Ensure `expected_cells` in sample_info.xlsx is reasonable for your experiment
   - Check that kallisto alignment completed successfully
   - Review alignment logs for low alignment rates

3. **Guide mapping issues**
   - Verify guide IDs match between guides.txt and guides_qc_reference.txt
   - Check guide sequences are correct length (typically 20bp for CRISPR)
   - Ensure guide FASTQ files are properly demultiplexed

4. **Memory errors**
   - Adjust resource allocations in config:
     ```yaml
     resources:
       alignment:
         mem_mb: 256000  # Increase to 256GB for large datasets
     ```

5. **Parse Bio barcode issues**
   - Ensure using correct kit files (96-well vs 48-well)
   - Verify barcode files have correct number of sequences
   - Check that replace and barcode files match your kit version

### Validation Steps

Before running the full pipeline:

1. **Check sample info:**
   ```python
   import pandas as pd
   df = pd.read_excel('sample_info.xlsx')
   # Check for required columns
   required = ['pool', 'sample', 'sample_id', 'sample_type', 'fastq_dir']
   missing = [c for c in required if c not in df.columns]
   print(f"Missing columns: {missing}")
   ```

2. **Verify FASTQ files exist:**
   ```python
   import os
   for _, row in df.iterrows():
       fastq_dir = row['fastq_dir']
       sample = row['sample']
       r1 = f"{fastq_dir}/{sample}_R1.fastq.gz"
       if not os.path.exists(r1):
           print(f"Missing: {r1}")
   ```

3. **Check guide consistency:**
   ```bash
   # Count guides in each file
   wc -l guides.txt
   cut -f1 guides_qc_reference.txt | sort -u | wc -l
   # These counts should match
   ```

---

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## License

[Add your license information here]

## Contact

For issues or questions:
1. Check the [Issues](https://github.com/your-repo/issues) page
2. Review log files in `logs_base/` for detailed error messages
3. Ensure all input files match the specifications above

## Acknowledgments

This pipeline uses:
- [kallisto](https://pachterlab.github.io/kallisto/) and [bustools](https://bustools.github.io/) for alignment
- [Snakemake](https://snakemake.readthedocs.io/) for workflow management
- [scanpy](https://scanpy.readthedocs.io/) for single-cell analysis
- [Parse Biosciences](https://www.parsebiosciences.com/) split-pool barcoding technology