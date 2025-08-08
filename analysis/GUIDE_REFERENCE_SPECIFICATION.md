# Guide Reference File Specification

## Overview
The guide reference file is a critical component of the GW_PERTURB pipeline that maps CRISPR guide RNA sequences to their target genes. This file is specified in `config.yaml` under `qc_analysis.guide_reference` and is used for quality control metrics and gene-level analysis.

## File Format
- **Type**: Tab-separated values (TSV)
- **Encoding**: UTF-8
- **Header**: Required (first row must contain column names)

## REQUIRED Columns

The pipeline **requires** exactly these two columns to function properly:

### 1. `ID` (Required)
- **Description**: Unique identifier for each guide RNA
- **Format**: String, typically following pattern `{SEQUENCE}_{GENE}_{SOURCE}`
- **Example**: `GAAAACAAACATGAGTGTGAA_CXCL11_crispick_guide4`
- **Notes**: 
  - Must be unique across all guides
  - Used as the primary key for guide-to-gene mapping
  - Must match the guide IDs in your kallisto/kb_python index

### 2. `gene` (Required)
- **Description**: Target gene name(s) for the guide
- **Format**: String, gene symbol(s)
- **Example**: `CXCL11` or `GENE1_AND_GENE2` for multi-gene targets
- **Special values**:
  - Empty string or NaN: Guide will be skipped in gene-level analysis
  - Multi-gene targets: Use the separator configured in `config.yaml` (see below)
- **Notes**: 
  - Gene names should match standard gene nomenclature (e.g., HGNC symbols)
  - Case-sensitive

## Multi-Gene Target Configuration

The separator for multi-gene targets is now specified in the guide reference file itself via a 'separator' column. If no separator column is present, the pipeline defaults to '_AND_' for backwards compatibility.

### File format with separator column:
```tsv
ID	gene	separator
GAAAACAAACATGAGTGTGAA_CXCL11_guide4	CXCL11	_AND_
MULTI_TARGET_GUIDE_1	MYC_AND_JUN_AND_FOS	_AND_
COMMA_SEPARATED_GUIDE	MYC,JUN,FOS	,
```

### Examples of different separator formats:
- **Default (`_AND_`)**: `MYC_AND_JUN_AND_FOS`
- **Comma (`,`)**: `MYC,JUN,FOS`
- **Semicolon (`;`)**: `MYC;JUN;FOS`
- **Pipe (`|`)**: `MYC|JUN|FOS`
- **None/empty**: All gene values treated as single genes (no splitting)

## Optional Columns

These columns are present in the example file but are NOT required by the pipeline:

- `seq (including G)`: Guide RNA sequence
- `Chromosome`: Genomic chromosome
- `Start`: Genomic start position
- `End`: Genomic end position
- `hgnc`: HGNC gene ID
- `ensg`: Ensembl gene ID
- `gene_type`: Type of gene (e.g., "protein coding")
- `strand`: Genomic strand (+/-)
- `genomic seq`: Genomic sequence context
- `pam`: PAM sequence
- `length`: Guide length
- `source`: Source of the guide design
- `validated_by`: Validation information
- `tss`: TSS information
- Various scoring columns (Hsu2013, DoenchCFD_specificityscore, etc.)

## Example Minimal Valid File

```tsv
ID	gene
GAAAACAAACATGAGTGTGAA_CXCL11_crispick_guide4	CXCL11
GAAAACAAACTGGGTCTCTAG_ATP8B4_crispick_guide19	ATP8B4
GAAAACAAAGAAAGGACAGAA_ECM2_crispick_guide9	ECM2
GAAAACAACAGCGCAGAGCCC_SYT4_crispick_guide9	SYT4
control_guide_1	
negative_control_2	
dual_target_guide_1	MYC_AND_JUN
```

### Example with different separators:

**Using comma separator (`multi_gene_separator: ","`):**
```tsv
ID	gene
guide_001	MYC,JUN,FOS
guide_002	TP53
guide_003	EGFR,KRAS
```

**No separator (`multi_gene_separator: null`):**
```tsv
ID	gene
guide_001	MYC
guide_002	TP53
guide_003	EGFR
```

## Related Files

### guides.txt (for kallisto/kb_python indexing)
The pipeline also uses a simpler two-column format file for creating the kallisto index:
- Specified in config.yaml as `guide_file` (full path)
- Format: `SEQUENCE<TAB>ID`
- Example:
```
GAAAACAAACATGAGTGTGAA	GAAAACAAACATGAGTGTGAA_CXCL11_crispick_guide4
GAAAACAAACTGGGTCTCTAG	GAAAACAAACTGGGTCTCTAG_ATP8B4_crispick_guide19
```

The `ID` column in this file MUST match the `ID` column in the guide reference TSV file.

## Pipeline Usage

1. **QC Metrics Calculation** (`calculate_qc_metrics_by_biological_sample.py`):
   - Loads the guide reference TSV file
   - Creates mapping: `guide_id -> [gene1, gene2, ...]`
   - Calculates `total_targeted_genes_cutoff{N}` metrics
   - Reports number of unique genes targeted

2. **Gene-Level Analysis**:
   - Aggregates guide-level counts to gene-level
   - Handles multi-gene targets appropriately
   - Filters out guides without gene annotations

## Validation Checklist

Before using a guide reference file:
- [ ] File is tab-separated (not comma-separated)
- [ ] Contains required columns: `ID` and `gene`
- [ ] All guide IDs are unique
- [ ] Multi-gene targets use `_AND_` separator
- [ ] Guide IDs match those in your kallisto index
- [ ] File is accessible at the path specified in config.yaml

## Error Handling

Common errors and solutions:

1. **KeyError: 'ID'** or **KeyError: 'gene'**
   - Solution: Ensure your TSV file has these exact column names (case-sensitive)

2. **Empty guide_to_genes mapping**
   - Solution: Check that `gene` column contains valid gene names (not all empty/NaN)

3. **Mismatch between guide counts and reference**
   - Solution: Verify guide IDs in reference match those in kallisto index

## Notes

- The pipeline is tolerant of missing genes (empty or NaN values)
- Additional columns are ignored but preserved if present
- The file can be very large (example file is ~45MB with extensive annotations)
- Consider keeping a full annotated version for reference and a minimal version for the pipeline