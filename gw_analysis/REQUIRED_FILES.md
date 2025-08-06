# Required Files for GW_PERTURB Pipeline

This document lists all files required to run the GW_PERTURB pipeline that are NOT included in the git repository due to size, privacy, or site-specific nature.

## Directory Structure Overview

```
GW_PERTURB/
├── gw_analysis/          # Pipeline code (in git)
├── data/                 # Raw sequencing data (NOT in git)
├── references/           # Reference files and indices (NOT in git)
└── analysis_results/     # Pipeline outputs (NOT in git)
```

## 1. Sample Configuration Files

### `../references/sample_info.xlsx`
**Purpose**: Master spreadsheet defining all samples and their metadata  
**Required Columns**:
- `sample_id`: Unique identifier (format: `pool:sample_name`)
- `pool`: Pool name (e.g., `pool1`, `pool2`)
- `sample`: Sample name without pool prefix (e.g., `gex_1`, `guide_1`)
- `sample_type`: Either `gex` or `guide` (rows with other values or NA are ignored)
- `fastq_dir`: Directory containing FASTQ files for this sample
- `expected_cells`: Number of expected cells for this sample
- `paired_guide_sample_id`: For GEX samples, the paired guide sample_id (required for guide-GEX pairing)
- `sample_to_well_mapping`: Name of plate map sheet in plate_maps.xlsx
- `min_umi_threshold`: Minimum UMI count for cell calling (optional)

**NA/Missing Value Handling**:
The following values are automatically parsed as NA/missing by pandas:
- Empty cells
- "NA", "N/A", "n/a", "NaN", "null", "NULL", "None", "#N/A", "#NA", "-NaN", "-nan", "nan"
- Rows with NA or any value other than `gex`/`guide` in `sample_type` are excluded from processing
- Guide samples must have a paired GEX sample; unpaired guides will cause an error

**Example**:
```
sample_id         pool    sample    sample_type  fastq_dir                     expected_cells  paired_guide_sample_id
pool1:gex_1       pool1   gex_1     gex         /path/to/pool1/fastqs/        5000            pool1:guide_1
pool1:guide_1     pool1   guide_1   guide       /path/to/pool1/fastqs/        5000            NA
pool1:control_1   pool1   control_1 NA          /path/to/pool1/fastqs/        NA              NA
```
Note: The `control_1` row would be ignored due to NA in sample_type.

### `../references/plate_maps.xlsx`
**Purpose**: Maps Parse Bio well barcodes to biological samples  
**Format**: Excel file with multiple sheets, each sheet representing a plate  
**Required Columns in Each Sheet**:
- `well`: Well position (e.g., `A1`, `B2`)
- `sample`: Biological sample name

## 2. Parse Bio Barcode Files

### `../references/barcodes.96.txt`
**Purpose**: Whitelist of valid Parse Bio cell barcodes  
**Format**: Text file with one barcode per line  
**Source**: Provided by Parse Biosciences for your kit version

### `../references/replace.96.txt`
**Purpose**: Barcode replacement/correction file for error correction  
**Format**: Tab-delimited file mapping incorrect to correct barcodes

### `../references/ParseBiosciences-Pipeline.1.2.1/`
**Purpose**: Parse Bio pipeline reference files  
**Required Files**:
- `splitpipe/barcodes/bc_data_n198_v5.csv` - 96-well barcode mappings
- `splitpipe/barcodes/bc_data_v2.csv` - 48-well barcode mappings

## 3. Guide Library Files

### Guide sequence file (specified in config.yaml as `guide_file`)
**Purpose**: Guide sequences for kallisto indexing  
**Format**: Tab-delimited with columns:
- Column 1: Guide sequence
- Column 2: Guide ID/name

**Note**: Specify the full path to this file in config.yaml as `guide_file`

### `../references/all_guides.20240917.tsv`
**Purpose**: Master reference of all guide sequences  
**Format**: Tab-delimited file for guide alignment

## 4. Genome Index Files

### `../references/nascent_all/`
**Purpose**: Kallisto indices for gene expression quantification  
**Required Files**:
- `transcriptome.idx` - Kallisto index for mature RNA
- `intron.idx` - Kallisto index for nascent RNA (introns)
- `cdna.fa` - cDNA sequences
- `intron.fa` - Intron sequences  
- `tr2g.txt` - Transcript to gene mapping

**Building Indices** (if not present):
```bash
# Example command to build nascent-aware index
kb ref \
  -i transcriptome.idx \
  -g tr2g.txt \
  -c1 cdna.fa \
  -c2 intron.fa \
  --workflow nucleus \
  genome.fa \
  genes.gtf
```

## 5. Gene Annotation Files

### `../references/gencode.v46.annotation.db`
**Purpose**: SQLite database with GENCODE gene annotations  
**Source**: Generated from GENCODE v46 GTF file  
**Used For**: Adding gene symbols, biotypes, chromosome locations

### `../references/regev_lab_cell_cycle_genes.txt`
**Purpose**: List of cell cycle genes for QC  
**Format**: One gene symbol per line  
**Source**: Regev lab cell cycle gene sets

### `../references/ribosomal.txt`
**Purpose**: List of ribosomal genes for QC  
**Format**: One gene symbol per line  
**Used For**: Calculating ribosomal gene percentage in cells

## 6. Raw Sequencing Data

### `../data/[pool_pattern]/`
**Purpose**: Raw FASTQ files from sequencing  
**Directory Structure**:
```
The fastq_dir column in sample_info.xlsx specifies where to find FASTQ files for each pool.
For example:
/path/to/pool1/
│   ├── gex_sample_R1.fastq.gz
│   ├── gex_sample_R2.fastq.gz
│   ├── guide_sample_R1.fastq.gz
│   ├── guide_sample_R2.fastq.gz
│   └── Undetermined_R1.fastq.gz  # For barcode recovery
```

**Naming Convention**: 
- R1: Contains cell barcode and UMI
- R2: Contains transcript/guide sequence
- File names must match `sample_id` in sample_info.xlsx

## 7. Optional Analysis Files

### `../references/split_lp_primer.xlsx`
**Purpose**: Maps i5/i7 index sequences to sample names  
**When Used**: ONLY needed for undetermined read recovery workflow  
**Not Required If**: You're only processing main reads (processing="raw")  
**Format**: Excel file with columns:
- Column B: Primer/sample name
- Column D: Index sequence
**Used By**: `create_undetermined_fastq` rule to assign recovered undetermined reads back to samples

## Setting Up for a New System

To run this pipeline on a new system, you need to:

1. **Clone the repository** (gets the code)
2. **Copy or link reference files** to `../references/`
3. **Copy or link raw data** to `../data/`
4. **Create sample_info.xlsx** with your experimental design
5. **Update config.yaml** with correct paths and parameters
6. **Set $SCRATCH** environment variable for temporary files

## File Size Estimates

- Genome indices: ~2-5 GB total
- Raw FASTQ files: 10-50 GB per sample
- Reference files: <100 MB
- Pipeline outputs: 5-20 GB per sample

## Data Privacy Notes

These files are not in git because they may contain:
- Proprietary guide sequences
- Patient/sample identifiers
- Large binary files inappropriate for version control
- Site-specific experimental data

## Verification Checklist

Run this command to verify all required files exist:
```bash
# From the gw_analysis directory
python -c "
import os
import sys

required_files = [
    '../references/sample_info.xlsx',
    '../references/plate_maps.xlsx',
    '../references/barcodes.96.txt',
    '../references/nascent_all/transcriptome.idx',
    '../references/gencode.v46.annotation.db',
]

missing = []
for f in required_files:
    if not os.path.exists(f):
        missing.append(f)
        print(f'✗ Missing: {f}')
    else:
        print(f'✓ Found: {f}')

if missing:
    print(f'\nMissing {len(missing)} required files')
    sys.exit(1)
else:
    print('\n✓ All required files present')
"
```

## Contact

For questions about obtaining these files:
- Parse Bio files: Contact Parse Biosciences support
- Genome indices: Build from GENCODE reference or contact bioinformatics core
- Guide libraries: Contact the lab that designed the CRISPR library