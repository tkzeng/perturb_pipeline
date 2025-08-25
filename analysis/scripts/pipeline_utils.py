"""
Unified pipeline utilities for single-cell perturbation data processing.

This module contains all utility functions and pipeline processing steps
for the single-cell perturbation analysis pipeline, organized by function.

Adapted for main pipeline structure - uses sample_info.tsv instead of hardcoded configs.
"""

import pandas as pd
import numpy as np
import subprocess
import sys
import time
import logging
from pathlib import Path
from typing import Optional
import scanpy as sc
import anndata as ad
import scipy
import gffutils
import gc
import sys
import os

# Configure logging only when needed, not at import time
logger = logging.getLogger(__name__)

def _configure_logging():
    """Configure logging only when explicitly needed."""
    if not logger.handlers:  # Only configure if not already configured
        handler = logging.StreamHandler(sys.stderr)  # Use stderr instead of stdout
        handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)

# GPU imports - only loaded if needed
try:
    import cupy as cp
    import cudf
    import rapids_singlecell as rsc
    import rmm

    GPU_AVAILABLE = True
except ImportError:
    GPU_AVAILABLE = False


# =============================================================================
# 0. HELPER FUNCTIONS FROM CONFIG
# =============================================================================

# Removed load_barcode_files - just use pd.read_csv directly at call sites


# Removed load_plate_maps - just use pd.read_excel directly at call sites


def create_well_to_sample_mapping_from_plates(plate_name, plate_maps_file):
    """Create mapping from wells to biological samples and all metadata for a specific plate."""
    import pandas as pd
    
    plate_maps = pd.read_excel(plate_maps_file, sheet_name=None)
    
    if plate_name not in plate_maps:
        raise ValueError(f"Plate {plate_name} not found in plate_maps file")
    
    plate_df = plate_maps[plate_name]
    
    # Create well→sample mapping
    well_to_sample = {}
    # Create well→metadata mapping for all columns
    well_to_metadata = {}
    
    for _, row in plate_df.iterrows():
        well = row['Well Position']
        biological_sample = row['Sample']
        
        if pd.isna(biological_sample):
            continue
            
        if well in well_to_sample:
            raise RuntimeError(f"CRITICAL: Well {well} maps to multiple samples in {plate_name}")
        
        well_to_sample[well] = biological_sample
        
        # Store all metadata columns (except Well Position)
        metadata = {}
        for col in plate_df.columns:
            if col != 'Well Position' and not pd.isna(row[col]):
                metadata[col] = row[col]
        well_to_metadata[well] = metadata
    
    return well_to_sample, well_to_metadata


# Removed get_expected_cells - just use the value directly at call sites


# =============================================================================
# 1. LOGGING & INFRASTRUCTURE
# =============================================================================


def log_print(*args, **kwargs):
    """Log message with INFO level and immediate flush."""
    _configure_logging()  # Ensure logging is configured
    message = " ".join(str(arg) for arg in args)
    logger.info(message)
    sys.stderr.flush()  # Flush stderr since we're logging there now


def align_datasets_by_cells(adata1, adata2, name1, name2):
    """Align two datasets keeping ALL cells from adata1 (GEX is authoritative).

    ⚠️  CRITICAL: GEX IS AUTHORITATIVE FOR CELL CALLING ⚠️
    ⚠️  KEEP ALL GEX CELLS, ADD GUIDE DATA WHERE AVAILABLE ⚠️
    ⚠️  DO NOT FILTER OUT GEX CELLS MISSING FROM GUIDE DATA ⚠️
    """
    common_cells = adata1.obs_names.intersection(adata2.obs_names)
    log_print(f"📊 Found {len(common_cells)} common cells between {name1} and {name2}")
    log_print(f"📊 {name1} total cells: {adata1.shape[0]}")
    log_print(f"📊 {name2} total cells: {adata2.shape[0]}")

    if len(common_cells) == 0:
        raise RuntimeError(
            f"CRITICAL: No common cells found between {name1} and {name2}. This indicates a fundamental data alignment issue."
        )

    # Keep ALL cells from adata1 (GEX is authoritative)
    adata1_all = adata1.copy()

    # Only get cells from adata2 that match adata1
    adata2_common = adata2[common_cells].copy()

    log_print(f"📊 Keeping ALL {adata1_all.shape[0]} cells from {name1} (authoritative)")
    log_print(f"📊 Adding data for {adata2_common.shape[0]} matching cells from {name2}")

    return adata1_all, adata2_common


# =============================================================================
# 2. CORE DATA PROCESSING (adapted for main pipeline)
# =============================================================================


def get_sample_pool_and_mapping(sample_id, sample_info_file):
    """Get pool and sample-to-well mapping for a specific sample from sample_info.tsv

    Args:
        sample_id: Sample ID to look up (e.g., "pool1:gex_1")
        sample_info_file: Path to sample_info.tsv

    Returns:
        tuple: (pool_name, sample_to_well_mapping_name)
    """
    sample_df = pd.read_csv(sample_info_file, sep='\t')
    sample_row = sample_df[sample_df['sample_id'] == sample_id]

    if sample_row.empty:
        raise ValueError(f"Sample {sample_id} not found in {sample_info_file}")

    # Pool is already in the dataframe, no need to extract from sample_id
    pool_name = sample_row.iloc[0]['pool']
    mapping_name = sample_row.iloc[0]['sample_to_well_mapping']

    return pool_name, mapping_name


def load_sample_info(sample_info_file):
    """Load and validate sample info with guide-GEX pairings.

    Args:
        sample_info_file: Path to sample_info.tsv

    Returns:
        pd.DataFrame: Sample info dataframe with validated columns
    """
    if not os.path.exists(sample_info_file):
        raise FileNotFoundError(f"Sample info file not found: {sample_info_file}")

    df = pd.read_csv(sample_info_file, sep='\t')

    # Validate required columns exist
    required_cols = ['sample_id', 'sample_type', 'pool']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Required columns missing from sample info: {missing_cols}")

    # Check if pairing columns exist (they're optional but if one exists, both should)
    pairing_cols = ['paired_guide_sample_id', 'paired_guide_pool']
    has_pairing = any(col in df.columns for col in pairing_cols)
    if has_pairing:
        missing_pairing = [col for col in pairing_cols if col not in df.columns]
        if missing_pairing:
            log_print(f"⚠️  Warning: Incomplete pairing columns. Missing: {missing_pairing}")

    return df


def get_guide_gex_pairings(sample_info_file):
    """Get guide-to-GEX and GEX-to-guide mappings from sample sheet.

    Uses paired_guide_sample_id column for explicit pairing.
    Since sample_ids now include pool prefix (e.g., "pool1:gex_1"), they are globally unique.

    Args:
        sample_info_file: Path to sample_info.tsv

    Returns:
        tuple: (guide_to_gex dict, gex_to_guide dict)
    """
    df = load_sample_info(sample_info_file)
    guide_to_gex = {}
    gex_to_guide = {}

    # Check if pairing column exists
    if 'paired_guide_sample_id' not in df.columns:
        raise ValueError(
            "CRITICAL: 'paired_guide_sample_id' column is REQUIRED in sample_info.tsv for guide-GEX pairing."
        )

    log_print("Building guide-GEX pairings from sample sheet")

    # Build mappings from paired_guide column
    for _, row in df.iterrows():
        if row['sample_type'] == 'gex' and 'paired_guide_sample_id' in row and pd.notna(row['paired_guide_sample_id']):
            gex_id = row['sample_id']  # e.g., "pool1:gex_1"
            guide_id = row['paired_guide_sample_id']  # e.g., "pool1:guide_1"

            # Validate guide exists
            guide_rows = df[df['sample_id'] == guide_id]
            if guide_rows.empty:
                raise ValueError(f"Paired guide {guide_id} not found for GEX sample {gex_id}")

            if guide_rows.iloc[0]['sample_type'] != 'guide':
                raise ValueError(f"Paired sample {guide_id} is not a guide sample (type: {guide_rows.iloc[0]['sample_type']})")

            # Simple mapping since sample_ids are now unique
            guide_to_gex[guide_id] = gex_id
            gex_to_guide[gex_id] = guide_id

    # Validate all guides have pairings
    all_guides = df[df['sample_type'] == 'guide']['sample_id'].tolist()
    unpaired_guides = [g for g in all_guides if g not in guide_to_gex]
    if unpaired_guides:
        raise ValueError(f"CRITICAL: Found unpaired guide samples: {unpaired_guides}. All guides must have paired GEX samples.")

    # Validate all GEX samples have pairings
    all_gex = df[df['sample_type'] == 'gex']['sample_id'].tolist()
    unpaired_gex = [g for g in all_gex if g not in gex_to_guide]
    if unpaired_gex:
        raise ValueError(f"CRITICAL: Found unpaired GEX samples: {unpaired_gex}. All GEX samples must have paired guide samples.")

    log_print(f"✅ Validated {len(gex_to_guide)} GEX-guide pairs")

    return guide_to_gex, gex_to_guide


def get_paired_sample(sample_id, sample_type, sample_info_file):
    """Get the paired sample for a given sample.

    Args:
        sample_id: Sample to find pair for
        sample_type: 'gex' or 'guide' (type of the input sample)
        sample_info_file: Path to sample_info.tsv

    Returns:
        str: Paired sample ID
    """
    guide_to_gex, gex_to_guide = get_guide_gex_pairings(sample_info_file)

    if sample_type == 'guide':
        if sample_id not in guide_to_gex:
            raise ValueError(f"No paired GEX sample found for guide {sample_id}")
        return guide_to_gex[sample_id]
    elif sample_type == 'gex':
        if sample_id not in gex_to_guide:
            raise ValueError(f"No paired guide sample found for GEX {sample_id}")
        return gex_to_guide[sample_id]
    else:
        raise ValueError(f"Invalid sample_type: {sample_type}. Must be 'gex' or 'guide'")


# REMOVED: process_single_library and process_guide_library functions
# These functions had complex path construction with hardcoded fallbacks.
# For combining sublibraries, use combine_sublibraries.py instead.


# Functions process_single_library and process_guide_library have been removed
# due to complex path construction and hardcoded fallbacks.



def add_guide_data(adata_gex, adata_guide):
    """Add separate guide data to GEX data in obsm format.

    CRITICAL UNDERSTANDING:
    - GEX and guide data come from the SAME physical cells (same cell prep)
    - Same barcodes appear in both datasets (same cell, different molecules)
    - After standardization, both have identical cell indices (e.g., "sample_1_BARCODE")
    - We merge by finding cells with identical obs.index (same barcode = same physical cell)

    Args:
        adata_gex: GEX expression data (genes × cells)
        adata_guide: Guide count data (guides × cells) - SAME CELLS as GEX

    Returns:
        Combined AnnData with:
        - X: GEX expression matrix
        - obsm["guide_counts"]: Guide count matrix for same cells
        - obs: Combined metadata from both datasets (no overwrites)
    """
    log_print("🧬 Adding guide data to GEX dataset...")

    # Align cells between datasets (GEX is authoritative)
    adata_gex_all, adata_guide_common = align_datasets_by_cells(adata_gex, adata_guide, "GEX dataset", "guide dataset")

    if not adata_gex_all:
        raise ValueError("Failed to align GEX and guide datasets")

    # ⚠️  CRITICAL: GEX IS AUTHORITATIVE FOR CELL CALLING ⚠️
    # ⚠️  KEEP ALL GEX CELLS, ADD GUIDE DATA WHERE AVAILABLE ⚠️
    log_print("🔄 Creating guide matrix for ALL GEX cells (GEX is authoritative)...")
    cells_gex_list = list(adata_gex_all.obs.index)
    cells_guide_set = set(adata_guide_common.obs.index)
    # Vectorized intersection preserving GEX order
    gex_index = pd.Index(cells_gex_list)
    common_mask = gex_index.isin(cells_guide_set)
    common_cells_ordered = gex_index[common_mask].tolist()

    log_print(f"📊 GEX cells (authoritative): {len(cells_gex_list):,}")
    log_print(f"📊 Guide cells available: {adata_guide_common.shape[0]:,}")
    log_print(f"📊 Cells with guide data: {len(common_cells_ordered):,}")
    log_print(f"📊 Cells missing guide data: {len(cells_gex_list) - len(common_cells_ordered):,}")

    # Create guide count matrix for ALL GEX cells (with zeros for missing cells)
    import scipy.sparse as sp

    n_gex_cells = adata_gex_all.shape[0]
    n_guides = adata_guide_common.shape[1]

    # Get indices of cells that have guide data using vectorized operations
    gex_index = pd.Index(adata_gex_all.obs.index)
    guide_index = pd.Index(adata_guide_common.obs.index)
    # Use get_indexer to vectorize the index lookup
    guide_indices = gex_index.get_indexer(guide_index)
    # Filter out any -1 values (cells not found)
    valid_mask = guide_indices != -1
    if not valid_mask.all():
        missing_count = (~valid_mask).sum()
        raise ValueError(f"Found {missing_count} guide cells not in GEX dataset")
    guide_indices = guide_indices[valid_mask].tolist()

    log_print(f"📊 Creating guide matrix for {n_gex_cells} GEX cells ({len(guide_indices)} with data)")

    # Get guide data and convert to coo format for efficient construction
    guide_data = adata_guide_common.X
    if not sp.issparse(guide_data):
        guide_data = sp.csr_matrix(guide_data)
    guide_coo = guide_data.tocoo()

    # Map cell indices from guide to GEX dataset coordinates
    row_mapping = np.array(guide_indices)
    new_rows = row_mapping[guide_coo.row]

    # Create sparse matrix directly with correct coordinates
    guide_counts_full = sp.csr_matrix((guide_coo.data, (new_rows, guide_coo.col)), shape=(n_gex_cells, n_guides))

    # Add guide data to obsm
    adata_gex_all.obsm["guide_counts"] = guide_counts_full
    adata_gex_all.uns["guide_names"] = adata_guide_common.var_names.tolist()
    log_print(f"✅ Added guide_counts matrix (zeros for {n_gex_cells - len(guide_indices)} cells)")

    # Transfer guide-specific obs columns (fill with NaN for missing cells)
    guide_specific_columns = ["guide_library_name", "guide_umi_counts"]
    for col in guide_specific_columns:
        if col in adata_guide_common.obs.columns:
            # Vectorized transfer: use reindex to align guide data to GEX cells
            guide_values = adata_guide_common.obs[col].reindex(adata_gex_all.obs.index)
            adata_gex_all.obs[col] = guide_values
            log_print(f"📋 Transferred guide column: {col}")

    cells_with_guides = len(adata_guide_common)
    total_gex_cells = len(adata_gex_all)
    log_print(f"📊 Added guide data: {n_guides} guides for {cells_with_guides}/{total_gex_cells} GEX cells")
    log_print(f"📋 Guide data stored in obsm['guide_counts'] (zeros for cells without guide data)")

    return adata_gex_all


# =============================================================================
# 2.5. BARCODE MAPPING FUNCTIONS (copied from statin_perturb)
# =============================================================================

def extract_and_validate_barcodes(adata, parsebio_barcodes):
    """Extract Round1 barcodes from cell names and validate.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    parsebio_barcodes : DataFrame
        Parse Bio barcode mapping CSV (maps sequences to wells)
    """
    log_print("🔍 Extracting Round1 barcodes from cell names...")

    # Extract round1 barcodes (last 8 characters) using vectorized string operations
    # CRITICAL: Use the obs_names directly as a pandas Index to preserve barcode names
    obs_index_series = pd.Series(adata.obs_names.values, index=adata.obs_names)
    round1_barcodes = obs_index_series.str[-8:]
    
    # Validate barcode extraction using vectorized length check
    barcode_lengths = round1_barcodes.str.len()
    invalid_mask = barcode_lengths != 8
    if invalid_mask.any():
        invalid_count = invalid_mask.sum()
        invalid_examples = round1_barcodes[invalid_mask].head(5).tolist()
        raise RuntimeError(f"Found {invalid_count} invalid barcodes (not 8 chars): {invalid_examples}")
    
    adata.obs["round1"] = round1_barcodes.values

    log_print(f"📊 Extracted barcodes for {len(adata.obs_names)} cells")

    # Validate barcodes against known Round1 barcodes
    if parsebio_barcodes is not None:
        all_valid_barcodes = set(parsebio_barcodes["sequence"])

        valid_mask = adata.obs["round1"].isin(all_valid_barcodes)
        invalid_count = (~valid_mask).sum()

        log_print(f"📊 {valid_mask.sum()}/{len(adata.obs_names)} cells have valid Round1 barcodes")

        if invalid_count > 0:
            invalid_barcodes = adata.obs.loc[~valid_mask, "round1"].unique()[:10]
            raise RuntimeError(
                f"CRITICAL: {invalid_count} cells have invalid Round1 barcodes. Examples: {list(invalid_barcodes)}"
            )

        log_print("✅ All Round1 barcodes are valid")
    else:
        raise RuntimeError("CRITICAL: Could not load Round1 barcode file for validation")


def map_cells_to_samples_with_plate(adata, plate_name, plate_maps_file, parsebio_barcodes):
    """Map cells to biological samples using barcode mappings and plate maps.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    plate_name : str
        Name of the plate for mapping
    plate_maps_file : str
        Path to plate maps Excel file
    parsebio_barcodes : DataFrame
        Parse Bio barcode mapping CSV (maps sequences to wells)
    """
    log_print(f"🗺️ Mapping cells to biological samples using plate {plate_name}...")

    # Create reverse mapping for this specific plate
    well_to_sample, well_to_metadata = create_well_to_sample_mapping_from_plates(plate_name, plate_maps_file)

    # Extract and validate barcodes
    extract_and_validate_barcodes(adata, parsebio_barcodes)

    # Create barcode lookup dictionary
    barcodes_dict = parsebio_barcodes.set_index("sequence")["well"].to_dict()

    # Map barcodes to wells
    adata.obs["well"] = adata.obs["round1"].map(barcodes_dict)

    # Map wells to biological samples
    adata.obs["biological_sample"] = adata.obs["well"].map(well_to_sample)

    # Add all metadata columns from plate maps using vectorized operations
    # Convert well_to_metadata dict to DataFrame for efficient merging
    metadata_records = []
    for well, metadata in well_to_metadata.items():
        record = {'well': well}
        record.update(metadata)
        metadata_records.append(record)
    
    if metadata_records:
        metadata_df = pd.DataFrame(metadata_records)
        # Remove 'Sample' column as it's already stored as 'biological_sample'
        if 'Sample' in metadata_df.columns:
            metadata_df = metadata_df.drop(columns=['Sample'])
        
        # Merge metadata with adata.obs based on well column
        # This is a vectorized left join operation
        # CRITICAL: Preserve the original index (cell barcodes) after merge
        original_index = adata.obs.index
        adata.obs = adata.obs.merge(metadata_df, on='well', how='left', suffixes=('', '_meta'))
        adata.obs.index = original_index  # Restore barcode names (left merge preserves order)
        
        # Handle any duplicate columns from merge (shouldn't happen but being safe)
        duplicate_cols = [col for col in adata.obs.columns if col.endswith('_meta')]
        if duplicate_cols:
            adata.obs = adata.obs.drop(columns=duplicate_cols)

    # Validate mapping results
    unmapped_wells = adata.obs[adata.obs["well"].isna()]
    unmapped_samples = adata.obs[adata.obs["biological_sample"].isna()]

    if len(unmapped_wells) > 0:
        raise ValueError(f"Failed to map {len(unmapped_wells)} cells to wells")

    if len(unmapped_samples) > 0:
        raise ValueError(f"Failed to map {len(unmapped_samples)} cells to biological samples")

    mapped_cells = len(adata.obs) - len(unmapped_samples)
    log_print(f"📊 Mapped {mapped_cells}/{len(adata.obs)} cells to biological samples")

    # Add plate information
    adata.obs["plate"] = plate_name

    return True




# =============================================================================
# 3. CONFIGURATION HELPERS (adapted for main pipeline)
# =============================================================================


# get_expected_cells function is imported from config.py


# =============================================================================
# 4. QC METRICS & FILTERING
# =============================================================================


def add_mitochondrial_metrics(adata):
    """Add mitochondrial gene metrics."""
    log_print("🧬 Adding mitochondrial gene metrics...")

    # Identify mitochondrial genes using gene names
    if 'gene' not in adata.var.columns:
        raise RuntimeError(
            "CRITICAL: Gene symbols required for mitochondrial gene detection. Run gene annotation first."
        )

    mt_genes = adata.var['gene'].str.startswith("MT-", na=False)
    log_print(f"📊 Using gene symbols: found {mt_genes.sum()} MT genes")

    if mt_genes.sum() == 0:
        raise RuntimeError("CRITICAL: No mitochondrial genes found. MT genes are required for quality control.")

    # Calculate mitochondrial metrics
    adata.var["mt"] = mt_genes
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    log_print(f"📊 Added mitochondrial metrics:")
    log_print(f"    Mean MT%: {adata.obs['pct_counts_mt'].mean():.2f}%")
    log_print(f"    Median MT%: {adata.obs['pct_counts_mt'].median():.2f}%")


def apply_cell_filters(adata, config=None):
    """Apply standard cell filtering based on QC metrics from config."""
    log_print("🧹 Applying cell quality filters...")

    # Make a copy to avoid view modification warnings
    adata = adata.copy()

    initial_cells = adata.shape[0]

    # Apply filters from config
    if config and "cell_filtering" in config:
        thresholds = config["cell_filtering"]
    else:
        # Fallback to default thresholds
        thresholds = {
            "max_total_counts": 40000,
            "max_mt_percent": 20,
            "min_genes": 1000,
            "default_min_counts": 2000
        }

    if "max_total_counts" in thresholds:
        adata = adata[adata.obs["total_counts"] < thresholds["max_total_counts"]]
    if "max_mt_percent" in thresholds and "pct_counts_mt" in adata.obs.columns:
        adata = adata[adata.obs["pct_counts_mt"] < thresholds["max_mt_percent"]]
    if "min_genes" in thresholds and "n_genes_by_counts" in adata.obs.columns:
        adata = adata[adata.obs["n_genes_by_counts"] > thresholds["min_genes"]]

    # Sample-specific minimum count thresholds
    if "default_min_counts" in thresholds:
        adata.obs["cutoff"] = thresholds["default_min_counts"]
        if "special_samples" in thresholds:
            for sample, cutoff in thresholds["special_samples"].items():
                if "sample" in adata.obs.columns:
                    adata.obs.loc[adata.obs["sample"] == sample, "cutoff"] = cutoff
        adata = adata[adata.obs["total_counts"] > adata.obs["cutoff"]]

    final_cells = adata.shape[0]
    removed_cells = initial_cells - final_cells

    log_print(f"📊 Cell filtering results:")
    log_print(f"    Initial: {initial_cells} cells")
    log_print(f"    Final: {final_cells} cells")
    log_print(f"    Removed: {removed_cells} cells ({removed_cells/initial_cells*100:.1f}%)")

    return adata


# =============================================================================
# 5. GENE ANNOTATION & FILTERING
# =============================================================================


def add_comprehensive_gene_annotations_fast(adata, gene_annotation_table_path):
    """Add comprehensive gene annotations using pre-computed table.
    
    This is a FAST replacement for add_comprehensive_gene_annotations that uses
    a pre-computed TSV table instead of querying gffutils database.
    
    Parameters:
    -----------
    adata : AnnData
        The annotated data object
    gene_annotation_table_path : str
        Path to the pre-computed gene annotation table (TSV format)
    """
    log_print("🧬 Adding comprehensive gene annotations (fast version)...")
    
    # Load pre-computed annotations from TSV
    gene_annotations = pd.read_csv(gene_annotation_table_path, sep='\t', index_col='gene_id')
    
    log_print(f"  📋 Loaded annotations for {len(gene_annotations)} genes from table")
    
    # Initialize columns with defaults
    adata.var['gene_type'] = ''
    adata.var['hgnc_id'] = ''
    adata.var['gene_name'] = ''  # Gene symbols
    adata.var['chromosome'] = ''
    adata.var['ribosomal'] = False
    adata.var['cell_cycle'] = False
    adata.var['cell_cycle_phase'] = ''
    adata.var['pct_cells_expressed'] = 0.0
    
    # Calculate % cells expressed (using expression threshold > 0)
    if scipy.sparse.issparse(adata.X):
        cells_expressing = (adata.X > 0).sum(axis=0).A1
    else:
        cells_expressing = (adata.X > 0).sum(axis=0)
    
    adata.var['pct_cells_expressed'] = (cells_expressing / adata.n_obs) * 100
    
    # Perform vectorized merge to annotate all genes at once
    log_print("  🔄 Performing vectorized gene annotation merge...")
    
    # Keep track of original var columns
    original_cols = list(adata.var.columns)
    
    # Merge annotations - this is MUCH faster than iterating
    adata.var = adata.var.merge(
        gene_annotations[['gene_type', 'hgnc_id', 'gene_name', 'chromosome', 
                          'ribosomal', 'cell_cycle', 'cell_cycle_phase']],
        left_index=True, 
        right_index=True, 
        how='left',
        suffixes=('_old', '')
    )
    
    # Fill NaN values with defaults
    adata.var['gene_type'] = adata.var['gene_type'].fillna('')
    adata.var['hgnc_id'] = adata.var['hgnc_id'].fillna('')
    adata.var['gene_name'] = adata.var['gene_name'].fillna('')
    adata.var['chromosome'] = adata.var['chromosome'].fillna('')
    adata.var['ribosomal'] = adata.var['ribosomal'].fillna(False)
    adata.var['cell_cycle'] = adata.var['cell_cycle'].fillna(False)
    adata.var['cell_cycle_phase'] = adata.var['cell_cycle_phase'].fillna('')
    
    # Add 'gene' column as alias for gene_name (for compatibility)
    adata.var['gene'] = adata.var['gene_name']
    
    # Clean up any duplicate columns from merge
    for col in adata.var.columns:
        if col.endswith('_old'):
            adata.var.drop(columns=[col], inplace=True)
    
    # Count successful annotations
    annotated_count = (adata.var['gene_type'] != '').sum()
    success_rate = annotated_count / adata.n_vars
    
    log_print(f"  📊 Annotation results: {annotated_count}/{adata.n_vars} genes annotated ({success_rate:.1%})")
    
    # Validate annotation success
    if annotated_count == 0:
        raise RuntimeError(f"CRITICAL: Gene annotation completely failed. 0/{adata.n_vars} genes annotated.")
    
    if success_rate < 0.5:
        raise RuntimeError(f"CRITICAL: Gene annotation mostly failed. Only {success_rate:.1%} genes annotated.")
    
    # Calculate mitochondrial genes using annotated gene symbols
    log_print("  🧬 Calculating mitochondrial gene scores...")
    mt_genes = adata.var["gene"].str.startswith("MT-", na=False)
    log_print(f"  📊 Found {mt_genes.sum()} mitochondrial genes")
    
    if mt_genes.sum() == 0:
        raise RuntimeError("CRITICAL: No mitochondrial genes found. MT genes are required for quality control.")
    
    # Calculate ribosomal scores
    log_print("  🧬 Calculating ribosomal gene scores...")
    ribo_genes = adata.var['ribosomal'] == True
    log_print(f"  📊 Found {ribo_genes.sum()} ribosomal genes")
    
    if ribo_genes.sum() > 0:
        # Calculate ribosomal percentage
        adata.obs['pct_counts_ribo'] = (adata[:, ribo_genes].X.sum(axis=1) / adata.X.sum(axis=1)) * 100
        if hasattr(adata.obs['pct_counts_ribo'], 'A1'):
            adata.obs['pct_counts_ribo'] = adata.obs['pct_counts_ribo'].A1
        log_print(f"  ✅ Added ribosomal percentage using {ribo_genes.sum()} genes")
    else:
        raise RuntimeError("CRITICAL: No ribosomal genes found. Ribosomal gene metrics are required for quality control.")
    
    # Summary statistics
    log_print(f"  ✅ Annotated {annotated_count}/{adata.n_vars} genes")
    log_print(f"  📊 Gene types: {adata.var['gene_type'].value_counts().head().to_dict()}")
    log_print(f"  📊 Ribosomal: {adata.var['ribosomal'].sum()}")
    log_print(f"  📊 Cell cycle: {adata.var['cell_cycle'].sum()}")
    log_print(f"  📊 Chromosomes: {adata.var['chromosome'].nunique()} unique")
    log_print(f"  📊 Mean % cells expressed: {adata.var['pct_cells_expressed'].mean():.1f}%")


def add_comprehensive_gene_annotations(adata, ribosomal_genes_path, gene_database_path, cell_cycle_genes_path):
    """Add comprehensive gene annotations to var dataframe.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    ribosomal_genes_path : str
        Path to ribosomal genes file
    gene_database_path : str
        Path to gene database file
    cell_cycle_genes_path : str
        Path to cell cycle genes file
    """
    log_print("🧬 Adding comprehensive gene annotations...")

    # Create reference files dict for internal use
    REFERENCE_FILES = {
        "ribosomal_genes": ribosomal_genes_path,
        "gene_database": gene_database_path,
        "cell_cycle_genes": cell_cycle_genes_path
    }

    # Load ribosomal genes from tab-separated file - use "HGNC ID" column
    ribosomal_hgnc_ids = set()
    with open(REFERENCE_FILES["ribosomal_genes"]) as f:
        header = f.readline().strip().split('\t')
        hgnc_col = header.index("HGNC ID")
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) > hgnc_col:
                ribosomal_hgnc_ids.add(cols[hgnc_col])

    log_print(f"  📋 Loaded {len(ribosomal_hgnc_ids)} ribosomal HGNC IDs")
    log_print(f"  📋 Sample ribosomal HGNC IDs: {list(ribosomal_hgnc_ids)[:5]}")

    # Load gencode database
    db = gffutils.FeatureDB(REFERENCE_FILES["gene_database"])

    # Initialize annotation columns
    adata.var['gene_type'] = ''
    adata.var['hgnc_id'] = ''
    adata.var['gene'] = ''  # Gene symbols
    adata.var['chromosome'] = ''
    adata.var['ribosomal'] = False
    adata.var['cell_cycle'] = False
    adata.var['cell_cycle_phase'] = ''
    adata.var['pct_cells_expressed'] = 0.0

    # Calculate % cells expressed (using expression threshold > 0)
    log_print("  📊 Calculating % cells expressed...")
    cells_expressed = (adata.X > 0).sum(axis=0)
    if hasattr(cells_expressed, 'A1'):
        cells_expressed = cells_expressed.A1
    adata.var['pct_cells_expressed'] = (cells_expressed / adata.n_obs) * 100

    # Load cell cycle genes for annotation
    with open(REFERENCE_FILES["cell_cycle_genes"]) as f:
        cell_cycle_genes_list = [x.strip() for x in f]
    cell_cycle_genes = set(cell_cycle_genes_list)  # For fast lookup
    s_phase_genes = set(cell_cycle_genes_list[:43])  # First 43 are S phase

    # Annotate from gencode database first
    log_print("  🗃️  Loading gene annotations from gencode database...")
    log_print(f"  🔍 Sample gene IDs from adata.var.index: {list(adata.var.index[:5])}")

    annotated_count = 0
    missing_count = 0
    error_count = 0

    for i, gene_id in enumerate(adata.var.index):
        if i < 5:  # Debug first 5 genes
            log_print(f"  🔍 Trying to annotate gene {i+1}: {gene_id}")

        try:
            # Query database for this gene by ID
            gene = db[gene_id]
            if i < 5:
                log_print(f"    ✅ Found gene {gene_id} in database")

            # Extract annotations - handle missing attributes gracefully
            # NOTE: We use .get() here because not all genes have all attributes (e.g., some lack hgnc_id)
            # This is expected and OK - we should NOT raise errors for missing optional attributes
            gene_type = gene.attributes.get('gene_type', [''])[0]
            hgnc_id = gene.attributes.get('hgnc_id', [''])[0]
            gene_name = gene.attributes.get('gene_name', [''])[0]

            if i < 5:
                log_print(f"    📋 Gene info: type={gene_type}, hgnc={hgnc_id}, name={gene_name}, chr={gene.seqid}")

            # Update annotations
            adata.var.loc[gene_id, 'gene_type'] = gene_type
            adata.var.loc[gene_id, 'hgnc_id'] = hgnc_id
            adata.var.loc[gene_id, 'gene'] = gene_name  # Add gene symbols
            adata.var.loc[gene_id, 'chromosome'] = gene.seqid

            # Check ribosomal genes (by HGNC ID)
            if hgnc_id in ribosomal_hgnc_ids:
                adata.var.loc[gene_id, 'ribosomal'] = True

            # Check cell cycle genes (by gene name)
            if gene_name in cell_cycle_genes:
                adata.var.loc[gene_id, 'cell_cycle'] = True
                # Determine phase using pre-computed S phase set
                if gene_name in s_phase_genes:
                    adata.var.loc[gene_id, 'cell_cycle_phase'] = 'S'
                else:
                    adata.var.loc[gene_id, 'cell_cycle_phase'] = 'G2M'

            annotated_count += 1

        except Exception as e:
            raise RuntimeError(f"Failed to annotate gene {gene_id}: {e}") from e

    log_print(f"  📊 Annotation results: {annotated_count} success, {missing_count} missing, {error_count} errors")

    # CRITICAL: Fail if gene annotation completely failed
    if annotated_count == 0:
        raise RuntimeError(f"CRITICAL: Gene annotation completely failed. 0/{adata.n_vars} genes annotated. Database lookup is broken.")

    # CRITICAL: Fail if most genes failed annotation
    success_rate = annotated_count / adata.n_vars
    if success_rate < 0.5:  # Less than 50% success
        raise RuntimeError(f"CRITICAL: Gene annotation mostly failed. Only {success_rate:.1%} ({annotated_count}/{adata.n_vars}) genes annotated successfully.")

    # Calculate mitochondrial genes using annotated gene symbols
    log_print("  🧬 Calculating mitochondrial gene scores...")
    mt_genes = adata.var["gene"].str.startswith("MT-", na=False)
    log_print(f"  📊 Found {mt_genes.sum()} mitochondrial genes")

    if mt_genes.sum() > 0:
        mt_gene_names = adata.var_names[mt_genes].tolist()

        log_print(f"  ✅ Added MT percentage and score using {len(mt_gene_names)} genes")
    else:
        raise RuntimeError("CRITICAL: No mitochondrial genes found for scoring. MT genes are required for quality control.")

    # Calculate ribosomal scores using annotated genes
    log_print("  🧬 Calculating ribosomal gene scores...")
    ribo_genes = adata.var['ribosomal'] == True
    log_print(f"  📊 Found {ribo_genes.sum()} ribosomal genes")

    if ribo_genes.sum() > 0:
        ribo_gene_names = adata.var_names[ribo_genes].tolist()

        # Calculate ribosomal percentage
        adata.obs['pct_counts_ribo'] = (adata[:, ribo_genes].X.sum(axis=1) / adata.X.sum(axis=1)) * 100
        if hasattr(adata.obs['pct_counts_ribo'], 'A1'):
            adata.obs['pct_counts_ribo'] = adata.obs['pct_counts_ribo'].A1

        log_print(f"  ✅ Added ribosomal percentage and score using {len(ribo_gene_names)} genes")
    else:
        raise RuntimeError("CRITICAL: No ribosomal genes found for scoring. Ribosomal gene metrics are required for quality control.")

    # Summary statistics
    log_print(f"  ✅ Annotated {annotated_count}/{adata.n_vars} genes")
    log_print(f"  📊 Gene types: {adata.var['gene_type'].value_counts().head().to_dict()}")
    log_print(f"  📊 Ribosomal: {adata.var['ribosomal'].sum()}")
    log_print(f"  📊 Cell cycle: {adata.var['cell_cycle'].sum()}")
    log_print(f"  📊 Chromosomes: {adata.var['chromosome'].nunique()} unique")
    log_print(f"  📊 Mean % cells expressed: {adata.var['pct_cells_expressed'].mean():.1f}%")


def _score_gene_list(adata_main, gene_list, score_name, use_gpu=False):
    """Score a gene list on normalized, log-transformed counts (no scaling needed).

    Note: Scaling is NOT needed for gene scoring as recommended in:
    https://satijalab.org/seurat/articles/cell_cycle_vignette
    Gene scoring should be performed on normalized, log-transformed counts.
    """
    if len(gene_list) == 0:
        log_print(f"⚠️ No genes found for {score_name}")
        return

    # Use all genes as universe (no subset needed)
    log_print(f"📊 Scoring {len(gene_list)} genes for {score_name} on all {adata_main.n_vars} genes")

    # Score genes directly on normalized, log-transformed data (no scaling)
    log_print(f"🔍 DEBUG: adata_main.X.max() RIGHT before score_genes call: {adata_main.X.max()}")
    if use_gpu:
        rsc.tl.score_genes(adata_main, gene_list=gene_list, score_name=score_name)
    else:
        sc.tl.score_genes(adata_main, gene_list=gene_list, score_name=score_name)

    log_print(f"✅ {score_name}: scored {len(gene_list)} genes")


# =============================================================================
# 6. GUIDE PROCESSING
# =============================================================================


def filter_guides_by_reference(adata):
    """Filter guide data by reference guide list.

    This is a placeholder for future implementation when guide processing
    is needed in the main pipeline.
    """
    raise NotImplementedError("Guide filtering functionality not yet implemented for main pipeline")


# Guide analysis functions moved to scripts/guide_analysis.py


# =============================================================================
# 7. ADVANCED ANALYSIS (PCA, UMAP, clustering)
# =============================================================================


def normalize_and_preprocess(adata, use_gpu=False):
    """Normalize data and create HVG subset for analysis.

    Correct order: Raw UMIs -> Filter expressed genes (>0.1% cells) -> HVG detection (seurat_v3 on raw) -> Normalize -> Log1p -> Subset to HVG
    """
    log_print("🔄 Normalizing and preprocessing data...")

    # ⚠️  CRITICAL: NEVER USE adata.raw - IT CAUSES DE/SCORING FUNCTIONS TO USE WRONG DATA ⚠️
    # ⚠️  MANY SCANPY FUNCTIONS AUTOMATICALLY USE .raw WHEN PRESENT, CAUSING SILENT BUGS ⚠️
    # ⚠️  ORIGINAL adata IS KEPT THROUGHOUT - JUST TRANSFER RESULTS BACK AT THE END ⚠️

    # Check that gene annotations are available
    if "gene_type" not in adata.var.columns:
        raise RuntimeError("CRITICAL: gene_type annotation missing. Run gene annotation before preprocessing.")

    # Filter genes expressed in >0.1% of cells AND protein-coding genes only
    # NOTE: Protein-coding filter makes a significant difference in clustering results
    min_cells = max(1, int(len(adata) * 0.001))  # 0.1% of cells
    genes_expressed = (adata.X > 0).sum(axis=0).A1 >= min_cells
    protein_coding = adata.var["gene_type"] == "protein_coding"
    genes_expressed = genes_expressed & protein_coding

    log_print(f"🧹 Using {genes_expressed.sum()} protein-coding genes expressed in ≥{min_cells} cells (0.1%)")
    log_print(f"📊 Protein-coding filter removed {(~protein_coding).sum()} non-coding genes")

    # Create expressed genes dataset for HVG detection
    adata_expressed = adata[:, genes_expressed].copy()
    log_print(f"📊 Created expressed genes dataset: {adata_expressed.shape}")

    # Calculate HVGs on raw counts (Seurat v3 approach)
    # Note: hvg_input column not needed since adata_expressed contains exactly those genes
    # TODO: Make n_top_genes configurable (currently hardcoded to 2000)
    sc.pp.highly_variable_genes(adata_expressed, flavor="seurat_v3", n_top_genes=2000, subset=False)

    # Transfer HVG information to original dataset
    hvg_columns = [col for col in adata_expressed.var.columns if col.startswith(('highly_variable', 'means', 'variances', 'highly_variable_rank', 'highly_variable_nbatches', 'highly_variable_intersection'))]

    # Normalize all expressed genes after HVG detection
    # TODO: Make target_sum configurable (currently hardcoded to 1e4)
    sc.pp.normalize_total(adata_expressed, target_sum=1e4)
    sc.pp.log1p(adata_expressed)
    for col in hvg_columns:
        if col not in adata.var.columns:
            adata.var[col] = 0.0 if adata_expressed.var[col].dtype in ['float64', 'float32'] else False
        adata.var.loc[adata_expressed.var.index, col] = adata_expressed.var[col]
    log_print(f"📊 Normalized expressed genes dataset: {adata_expressed.shape}")

    # Then subset to HVG for dimensionality reduction
    adata_hvg = adata_expressed[:, adata_expressed.var.highly_variable].copy()
    log_print(f"📊 Created HVG subset: {adata_hvg.shape}")

    return adata_hvg, adata_expressed


def perform_dimensionality_reduction(adata, use_gpu=False, n_neighbors=15, n_pcs=15):
    """Perform PCA and UMAP on data."""
    log_print(f"📊 Performing dimensionality reduction (n_neighbors={n_neighbors}, n_pcs={n_pcs})...")

    # PCA
    # TODO: Make n_comps configurable (currently hardcoded to 50)
    if use_gpu and GPU_AVAILABLE:
        rsc.pp.pca(adata, n_comps=50, use_highly_variable=True)
    else:
        sc.pp.pca(adata, n_comps=50, use_highly_variable=True)

    # Compute neighborhood graph
    if use_gpu and GPU_AVAILABLE:
        rsc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    # UMAP
    if use_gpu and GPU_AVAILABLE:
        rsc.tl.umap(adata)
    else:
        sc.tl.umap(adata)

    log_print("✅ Dimensionality reduction completed")

    return adata


def perform_clustering(adata, use_gpu=False, target_clusters=3):
    """Perform Leiden clustering discovering up to 6 clusters, starting at 0.05 resolution."""
    log_print(f"🎯 Performing Leiden clustering (discovering up to 6 clusters)...")

    # TODO: Make clustering parameters configurable (resolution start, max_clusters, limit)
    # Start at 0.05 and increment until we find up to 6 clusters
    resolution = 0.05
    max_clusters = 6
    best_resolution = 0.05

    while resolution <= 2.0:  # Safety limit
        if use_gpu and GPU_AVAILABLE:
            rsc.tl.leiden(adata, resolution=resolution, key_added=f"leiden_{resolution}")
        else:
            sc.tl.leiden(adata, resolution=resolution, key_added=f"leiden_{resolution}")

        n_clusters = len(adata.obs[f"leiden_{resolution}"].unique())
        log_print(f"  Resolution {resolution}: {n_clusters} clusters")

        best_resolution = resolution

        # Stop if we've reached the maximum number of clusters
        if n_clusters >= max_clusters:
            log_print(f"🎯 Reached {n_clusters} clusters at resolution {resolution}, stopping")
            break

        # Increment resolution
        resolution = round(resolution + 0.05, 2)

    log_print(f"✅ Clustering complete: resolution {best_resolution} with {len(adata.obs[f'leiden_{best_resolution}'].unique())} clusters")

    # Rename leiden columns to clustersN_resX format and remove originals
    log_print("🔄 Renaming leiden columns to clustersN_resX format...")
    # Use pandas string operations for efficient filtering
    leiden_mask = adata.obs.columns.str.startswith('leiden_')
    leiden_cols = adata.obs.columns[leiden_mask].tolist()
    for col in leiden_cols:
        resolution_val = col.replace('leiden_', '')
        n_clusters = len(adata.obs[col].unique())
        new_col_name = f"clusters{n_clusters}_res{resolution_val}"
        adata.obs[new_col_name] = adata.obs[col]
        log_print(f"  Renamed {col} → {new_col_name}")
    
    # Now remove the original leiden columns
    adata.obs.drop(columns=leiden_cols, inplace=True)
    log_print(f"🧹 Removed {len(leiden_cols)} original leiden columns")

    return adata


def ensure_cpu_data(adata, operation_name="operation"):
    """Ensure data is on CPU for operations that require it."""
    if GPU_AVAILABLE and hasattr(adata.X, "get"):  # CuPy array
        log_print(f"🔄 Moving data to CPU for {operation_name}...")
        adata_cpu = adata.copy()
        adata_cpu.X = adata_cpu.X.get()  # Move from GPU to CPU
        return adata_cpu
    return adata


# =============================================================================
# 6. GPU UTILITIES
# =============================================================================


def get_gpu_memory_usage():
    """Get GPU memory usage for all available GPUs."""
    if not GPU_AVAILABLE:
        return {}

    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=index,memory.used,memory.total", "--format=csv,noheader,nounits"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode == 0:
            gpu_info = {}
            for line in result.stdout.strip().split("\n"):
                if line.strip():
                    parts = line.split(", ")
                    if len(parts) == 3:
                        gpu_id = int(parts[0])
                        used_mb = int(parts[1])
                        total_mb = int(parts[2])
                        gpu_info[gpu_id] = {
                            "used_gb": used_mb / 1024,
                            "total_gb": total_mb / 1024,
                            "free_gb": (total_mb - used_mb) / 1024,
                        }
            return gpu_info
    except Exception as e:
        raise RuntimeError(f"Failed to get GPU memory info: {e}") from e

    return {}


def find_free_gpu(threshold_gb=2.0):
    """Find a GPU with less than threshold_gb used memory."""
    gpu_info = get_gpu_memory_usage()

    for gpu_id, info in gpu_info.items():
        if info["used_gb"] < threshold_gb:
            log_print(f"🎮 Found free GPU {gpu_id}: {info['used_gb']:.1f}GB used, {info['free_gb']:.1f}GB free")
            return gpu_id

    return None


def setup_gpu(use_gpu=False, gpu_threshold_gb=2.0):
    """Setup GPU for computations if requested and available.

    ⚠️  CRITICAL: DO NOT CHANGE THIS FUNCTION TO FALL BACK TO CPU ⚠️
    ⚠️  WHEN GPU IS REQUESTED, THE PIPELINE MUST FAIL IF GPU UNAVAILABLE ⚠️
    ⚠️  NO CPU FALLBACK ALLOWED - RAISE ERRORS INSTEAD ⚠️
    """
    if not use_gpu:
        log_print("🖥️  CPU-only mode requested")
        return False

    # ⚠️  DO NOT CHANGE: Must raise error, not fall back to CPU ⚠️
    if not GPU_AVAILABLE:
        raise RuntimeError(
            "GPU libraries not available but GPU mode was requested. Install RAPIDS or disable GPU mode."
        )

    # Find available GPU
    gpu_id = find_free_gpu(gpu_threshold_gb)
    # ⚠️  DO NOT CHANGE: Must raise error, not fall back to CPU ⚠️
    if gpu_id is None:
        raise RuntimeError(
            f"No GPU with <{gpu_threshold_gb}GB usage found but GPU mode was requested. Free up GPU memory or disable GPU mode."
        )

    try:
        # Set GPU device
        cp.cuda.Device(gpu_id).use()

        # Initialize memory pool
        mempool = cp.get_default_memory_pool()
        mempool.set_limit(size=None)  # No limit

        log_print(f"🎮 GPU {gpu_id} initialized successfully")
        return True

    except Exception as e:
        # ⚠️  DO NOT CHANGE: Must raise error, not fall back to CPU ⚠️
        raise RuntimeError(f"Failed to initialize GPU {gpu_id} but GPU mode was requested: {e}")


def cleanup_gpu_memory():
    """Clean up GPU memory."""
    if GPU_AVAILABLE:
        try:
            # Clear cupy cache
            cp.get_default_memory_pool().free_all_blocks()
            log_print("🧹 GPU memory cleaned up")
        except Exception as e:
            raise RuntimeError(f"GPU cleanup failed: {e}") from e


# =============================================================================
# 7. SAMPLE MAPPING HELPERS (for future use)
# =============================================================================


def map_cells_to_samples_from_sample_info(adata, sample_info_file):
    """Map cells to samples using sample_info.tsv and sample-to-well mappings.

    This is a placeholder for future implementation when you want to add
    sample annotation functionality similar to statin_perturb.
    """
    raise NotImplementedError("Sample mapping functionality not yet implemented for main pipeline")


# =============================================================================
# 8. GUIDE THRESHOLD UTILITIES
# =============================================================================


# Guide cutoff utilities moved to scripts/guide_analysis.py


