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
- `sample_type`: Either `gex` or `guide`
- `expected_cells`: Number of expected cells for this sample
- `sample_to_well_mapping`: Name of plate map sheet in plate_maps.xlsx
- `biological_sample`: Biological replicate name
- `min_umi_threshold`: Minimum UMI count for cell calling (optional)
- `emptydrops_lower`: Lower bound for EmptyDrops algorithm (optional)

**Example**:
```
sample_id         pool    sample_type  expected_cells  sample_to_well_mapping
pool1:gex_rep1    pool1   gex         5000            plate1_mapping
pool1:guide_rep1  pool1   guide       5000            plate1_mapping
```

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

### `../references/guides/[guide_set].txt`
**Purpose**: Guide sequences and metadata for your CRISPR library  
**Format**: Tab-delimited with columns:
- `id`: Unique guide identifier
- `seq`: Guide sequence (20-21bp)
- `gene`: Target gene name
- Other metadata columns as needed

**Note**: The `guide_set` parameter in config.yaml should match the filename (without .txt)

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
../data/
├── May2025_GW_pool1/
│   ├── pool1_gex_rep1_R1.fastq.gz
│   ├── pool1_gex_rep1_R2.fastq.gz
│   ├── pool1_guide_rep1_R1.fastq.gz
│   ├── pool1_guide_rep1_R2.fastq.gz
│   └── Undetermined_R1.fastq.gz  # For barcode recovery
└── May2025_GW_pool2/
    └── ... (similar structure)
```

**Naming Convention**: 
- R1: Contains cell barcode and UMI
- R2: Contains transcript/guide sequence
- File names must match `sample_id` in sample_info.xlsx

## 7. Optional Analysis Files

### `../references/split_lp_primer.xlsx`
**Purpose**: Primer sequences for demultiplexing (if applicable)  
**Format**: Excel file with primer information

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