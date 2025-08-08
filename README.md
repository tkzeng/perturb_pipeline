# Perturb-seq pipeline

The documentation below uses the K562 CRISPR screen configuration (`analysis/config.yann_k562.yaml`) as a concrete example.

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
cp analysis/config.yann_k562.yaml analysis/config.my_experiment.yaml
# Edit analysis/config.my_experiment.yaml with your paths and parameters

# 4. Prepare your input files (see Required Input Files section)

# 5. Run the pipeline
cd analysis
# Dry run to check what will be executed
CONFIG=config.yann_k562.yaml ./submit.sh --dry-run
# Submit using slurm
CONFIG=config.my_experiment.yaml ./submit.sh
```

## Required Input Files

Before running the pipeline, you need to prepare the following files:

### 1. Sample Information Excel File (`sample_info.xlsx`)
**Purpose:** Defines all samples, their properties, and relationships  
**Location:** Referenced in config as `sample_info_file`

**Required columns (must be present for all samples):**
- `pool`: Pool identifier (e.g., "1k_pilot", "pool1", "pool2")
- `sample`: Sample name within the pool (e.g., "sample1_gex", "sample1_grna")
- `sample_id`: Unique identifier combining pool and sample (e.g., "1k_pilot:sample1_gex")
- `sample_type`: Either "gex" (gene expression) or "guide" (guide RNA)
- `fastq_dir`: Directory containing FASTQ files for this sample
- `expected_cells`: Expected number of cells for this library (e.g., 5000)
- `min_umi_threshold`: Minimum UMI count threshold for cell calling (e.g., 100)
- `sample_to_well_mapping`: Name of the plate mapping sheet (e.g., "plate1")
- `paired_guide_sample_id`: Full ID of paired guide library (for GEX samples, leave empty for guide samples)

**Optional columns (required only for specific features):**

*For undetermined read recovery:*
- `i7_barcode`: i7 index sequence (8bp, e.g., "AAGGCTAT")
- `i5_barcode`: i5 index sequence (8bp, e.g., "AAACATCG")

**Columns not used by pipeline (can be included for documentation):**
- `notes`: Optional notes about the sample
- `read1_length`: Length of Read 1 (typically 150)
- `read2_length`: Length of Read 2 (typically 150)
- `paired_guide_pool`: Pool of paired guide (derived from paired_guide_sample_id)
- `paired_guide_sample`: Sample name of paired guide (not used)

### 2. Plate Mapping Excel File (`plate_map.xlsx`)
**Purpose:** Maps well positions to biological samples  
**Location:** Referenced in config as `plate_maps_file`

**Structure:**
- Each sheet represents a different plate (sheet name matches `sample_to_well_mapping` in sample_info)
- Required columns per sheet:
  - `Well Position`: Well identifier (A1-H12 for 96-well plate)
  - `Sample`: Biological sample name

**Optional columns:**
- `Rep`: Replicate identifier (stored as metadata but not used by pipeline)
- Any other columns will be stored as cell metadata in the AnnData object

### 3. Guide Library File (`guides.txt`)
**Purpose:** Defines all guide RNA sequences for kallisto indexing  
**Location:** Referenced in config as `guide_file`

**Format:** Tab-separated text file with columns:
- Column 1: Guide sequence (typically 20-21bp)
- Column 2: Guide ID/name

### 4. Guide Reference File (`guides_qc_reference.txt`)
**Purpose:** Maps guide IDs to target genes for QC analysis  
**Location:** Referenced in config as `guide_reference`

**Required columns (tab-separated):**
- `ID`: Guide identifier (must match IDs in guides.txt)
- `gene`: Target gene name(s) - use separator for multi-gene targets

**Optional columns:**
- `separator`: Character to split multi-gene targets (defaults to '_AND_')
- Any other columns are ignored

### 5. Fill out the `config.my_experiment.yaml` file
