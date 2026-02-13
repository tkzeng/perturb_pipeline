"""
Unified pipeline utilities for single-cell perturbation data processing.

This module contains all utility functions and pipeline processing steps
for the single-cell perturbation analysis pipeline, organized by function.

Adapted for main pipeline structure - uses sample_info.xlsx instead of hardcoded configs.
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
        raise ValueError(f"Plate {plate_name} not found in plate_maps.xlsx")
    
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
    """Get pool and sample-to-well mapping for a specific sample from sample_info.xlsx

    Args:
        sample_id: Sample ID to look up (e.g., "pool1:gex_1")
        sample_info_file: Path to sample_info.xlsx

    Returns:
        tuple: (pool_name, sample_to_well_mapping_name)
    """
    sample_df = pd.read_excel(sample_info_file)
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
        sample_info_file: Path to sample_info.xlsx

    Returns:
        pd.DataFrame: Sample info dataframe with validated columns
    """
    if not os.path.exists(sample_info_file):
        raise FileNotFoundError(f"Sample info file not found: {sample_info_file}")

    df = pd.read_excel(sample_info_file)

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
        sample_info_file: Path to sample_info.xlsx

    Returns:
        tuple: (guide_to_gex dict, gex_to_guide dict)
    """
    df = load_sample_info(sample_info_file)
    guide_to_gex = {}
    gex_to_guide = {}

    # Check if pairing column exists
    if 'paired_guide_sample_id' not in df.columns:
        raise ValueError(
            "CRITICAL: 'paired_guide_sample_id' column is REQUIRED in sample_info.xlsx for guide-GEX pairing."
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
        sample_info_file: Path to sample_info.xlsx

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
    common_cells_ordered = [cell for cell in cells_gex_list if cell in cells_guide_set]

    log_print(f"📊 GEX cells (authoritative): {len(cells_gex_list):,}")
    log_print(f"📊 Guide cells available: {adata_guide_common.shape[0]:,}")
    log_print(f"📊 Cells with guide data: {len(common_cells_ordered):,}")
    log_print(f"📊 Cells missing guide data: {len(cells_gex_list) - len(common_cells_ordered):,}")

    # Create guide count matrix for ALL GEX cells (with zeros for missing cells)
    import scipy.sparse as sp

    n_gex_cells = adata_gex_all.shape[0]
    n_guides = adata_guide_common.shape[1]

    # Get indices of cells that have guide data
    cell_to_idx = {cell: i for i, cell in enumerate(adata_gex_all.obs.index)}
    guide_indices = [cell_to_idx[cell] for cell in adata_guide_common.obs.index]

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

    # Extract round1 barcodes (last 8 characters)
    round1_barcodes = [x[-8:] for x in adata.obs_names]

    # Validate barcode extraction
    invalid_barcodes = [bc for bc in round1_barcodes if len(bc) != 8]
    if invalid_barcodes:
        raise RuntimeError(f"Found {len(invalid_barcodes)} invalid barcodes (not 8 chars): {invalid_barcodes[:5]}")

    adata.obs["round1"] = round1_barcodes

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

    # Add all metadata columns from plate maps
    for _, row_idx in adata.obs.iterrows():
        well = row_idx.get("well")
        if well and well in well_to_metadata:
            metadata = well_to_metadata[well]
            for col_name, value in metadata.items():
                # Skip 'Sample' as it's already stored as 'biological_sample'
                if col_name != 'Sample':
                    adata.obs.at[row_idx.name, col_name] = value

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
    # make sure both mouse or human can be identified
    mt_genes = adata.var["gene"].astype(str).str.upper().str.startswith("MT-", na=False)
    #mt_genes = adata.var['gene'].str.startswith("MT-", na=False)
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

def add_comprehensive_gene_annotations(adata, ribosomal_genes_path, gene_database_path, cell_cycle_genes_path, species = 'human'):
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
    if species not in ("human", "mouse"):
        raise ValueError(f"Unknown species '{species}'. Use 'human' or 'mouse'.")
    log_print(f"🧬 Target species: {species}")
    
    # Create reference files dict for internal use
    REFERENCE_FILES = {
        "ribosomal_genes": ribosomal_genes_path,
        "gene_database": gene_database_path,
        "cell_cycle_genes": cell_cycle_genes_path
    }
    # ---- init annotation columns ----
    adata.var["gene_type"] = ""
    adata.var["gene"] = ""  # gene symbols
    adata.var["chromosome"] = ""
    adata.var["ribosomal"] = False
    adata.var["cell_cycle"] = False
    adata.var["cell_cycle_phase"] = ""
    adata.var["pct_cells_expressed"] = 0.0

    if species == "human":
        adata.var['hgnc_id'] = ''
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

    else:
        with open(REFERENCE_FILES["ribosomal_genes"]) as f:
            ribosomal_genes_list = [x.strip() for x in f]
            ribosomal_gene_symbols = set(ribosomal_genes_list)  # For fast lookup
        log_print(f"  📋 Loaded {len(ribosomal_gene_symbols)} ribosomal gene")
        log_print(f"  📋 Sample ribosomal gene: {list(ribosomal_gene_symbols)[:5]}")

    # Load gencode database
    db = gffutils.FeatureDB(REFERENCE_FILES["gene_database"])

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
        if i < 5:
            log_print(f"  🔍 Trying to annotate gene {i+1}: {gene_id}")
        try:
            try:
                gene = db[gene_id]
            except gffutils.exceptions.FeatureNotFoundError:
                clean_id = gene_id.split(".", 1)[0]
                try:
                    gene = db[clean_id]
                    if i < 5:
                        log_print(f"    ✅ Found gene {gene_id} in database (using ID {clean_id})")
                except gffutils.exceptions.FeatureNotFoundError:
                    missing_count += 1
                    continue

            # Extract annotations - handle missing attributes gracefully
            # NOTE: We use .get() here because not all genes have all attributes (e.g., some lack hgnc_id)
            # This is expected and OK - we should NOT raise errors for missing optional attributes
            gene_type = gene.attributes.get('gene_biotype', [''])[0]
            gene_name = gene.attributes.get('gene_name', [''])[0]

            if i < 5:
                log_print(f"    📋 Gene info: type={gene_type}, name={gene_name}, chr={gene.seqid}")

            # Update annotations
            adata.var.loc[gene_id, 'gene_type'] = gene_type
            adata.var.loc[gene_id, 'gene'] = gene_name  # Add gene symbols
            adata.var.loc[gene_id, 'chromosome'] = gene.seqid

            # Check ribosomal
            if species == "human":
                hgnc_id = gene.attributes.get("hgnc_id", [""])[0]
                adata.var.at[gene_id, "hgnc_id"] = hgnc_id
                if hgnc_id and (hgnc_id in ribosomal_hgnc_ids):
                    adata.var.at[gene_id, "ribosomal"] = True
            else:
                if gene_name and (gene_name in ribosomal_gene_symbols):
                    adata.var.at[gene_id, "ribosomal"] = True

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
    mt_genes = adata.var["gene"].astype(str).str.upper().str.startswith("MT-", na=False)

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


def calculate_guides_per_cell(guide_counts, cutoff):
    """Calculate guides per cell at a given cutoff.

    Args:
        guide_counts: Sparse matrix of guide counts (cells x guides)
        cutoff: Minimum count threshold for guide detection

    Returns:
        numpy array of guides per cell
    """
    binary_guides = (guide_counts >= cutoff).astype(int)
    return np.array(binary_guides.sum(axis=1)).flatten()


def calculate_guide_metrics_for_cells(guide_adata, cell_barcodes, guide_cutoffs, guide_to_genes=None):
    """Calculate guide metrics for a specific set of cells, including those with no guides.

    CRITICAL: This function includes ALL provided cells in calculations, even if they
    don't exist in the guide data. Missing cells contribute 0 UMIs to averages.

    Args:
        guide_adata: Guide count AnnData object (any subset)
        cell_barcodes: List/array of cell barcodes to calculate metrics for
        guide_cutoffs: List of UMI cutoffs for guide detection
        guide_to_genes: Optional dict mapping guide IDs to gene lists

    Returns:
        dict: Metrics including:
            - n_cells_total: Total cells provided
            - n_cells_in_guide_data: Cells found in guide data
            - guide_umis_per_cell: Mean guide UMIs (including 0s)
            - For each cutoff:
                - guides_per_cell_cutoff{N}: Mean guides per cell
                - fraction_cells_with_guides_cutoff{N}: Fraction with ≥1 guide
                - umis_per_guide_per_cell_cutoff{N}: Mean UMIs per guide (cells with guides only)
                - total_guides_detected_cutoff{N}: Unique guides across all cells
                - total_guide_detections_cutoff{N}: Total guide detections across all cells (sum)
                - total_targeted_genes_cutoff{N}: Unique genes targeted (if guide_to_genes provided)
    """
    metrics = {}
    n_cells_total = len(cell_barcodes)
    metrics['n_cells_total'] = n_cells_total

    if n_cells_total == 0:
        raise ValueError("No cell barcodes provided")

    # Find which cells exist in guide data using vectorized operations
    barcode_series = pd.Series(range(guide_adata.n_obs), index=guide_adata.obs_names)
    cell_indices_series = barcode_series.reindex(cell_barcodes)

    # Separate found and missing cells
    found_mask = ~cell_indices_series.isna()
    n_cells_found = found_mask.sum()
    metrics['n_cells_in_guide_data'] = int(n_cells_found)

    log_print(f"Found {n_cells_found:,} of {n_cells_total:,} cells in guide data ({n_cells_found/n_cells_total*100:.1f}%)")

    # Initialize arrays for ALL cells (including missing ones)
    all_umis_per_cell = np.zeros(n_cells_total)

    # If we have cells in guide data, get their counts
    if n_cells_found > 0:
        cell_indices = cell_indices_series.dropna().astype(int).values
        guide_counts_found = guide_adata.X[cell_indices]

        # Calculate UMIs for found cells
        umis_per_found_cell = np.array(guide_counts_found.sum(axis=1)).flatten()

        # Place these values in the correct positions in the full array
        all_umis_per_cell[found_mask] = umis_per_found_cell

        # Store guide counts for cutoff calculations
        guide_counts_sparse = scipy.sparse.csr_matrix(guide_counts_found)
    else:
        guide_counts_sparse = None

    # Calculate mean UMIs per cell (including cells with 0)
    metrics['guide_umis_per_cell'] = float(np.mean(all_umis_per_cell))

    # Calculate metrics for each cutoff
    for cutoff in guide_cutoffs:
        # Initialize arrays for ALL cells
        all_guides_per_cell = np.zeros(n_cells_total)

        if n_cells_found > 0 and guide_counts_sparse is not None:
            # Calculate guides per cell for found cells
            guides_per_found_cell = calculate_guides_per_cell(guide_counts_sparse, cutoff)
            all_guides_per_cell[found_mask] = guides_per_found_cell

            # Binary guide matrix for diversity calculations
            binary_guides = (guide_counts_sparse >= cutoff).astype(int)
            guides_detected = np.array(binary_guides.sum(axis=0)).flatten() > 0
            n_guides_detected = int(np.sum(guides_detected))
            
            # Total guide detections across all cells (sum of binary matrix)
            total_guide_detections = int(binary_guides.sum())
        else:
            n_guides_detected = 0
            total_guide_detections = 0

        # Mean guides per cell (including cells with 0)
        metrics[f'guides_per_cell_cutoff{cutoff}'] = float(np.mean(all_guides_per_cell))

        # Fraction of cells with guides (out of ALL cells)
        cells_with_guides = all_guides_per_cell > 0
        metrics[f'fraction_cells_with_guides_cutoff{cutoff}'] = float(np.sum(cells_with_guides) / n_cells_total)

        # UMIs per guide per cell (only for cells with guides)
        if np.any(cells_with_guides):
            # Get guide counts for cells with guides
            cells_with_guides_indices = np.where(found_mask & cells_with_guides)[0]
            if len(cells_with_guides_indices) > 0:
                # Map back to guide data indices
                guide_data_indices = cell_indices_series.iloc[cells_with_guides_indices].astype(int).values
                guide_counts_subset = guide_adata.X[guide_data_indices]

                # Only count UMIs from guides above cutoff
                binary_subset = (guide_counts_subset >= cutoff)
                filtered_counts = guide_counts_subset.multiply(binary_subset)
                umis_per_cell_with_guides = np.array(filtered_counts.sum(axis=1)).flatten()
                guides_per_cell_with_guides = all_guides_per_cell[cells_with_guides]

                umis_per_guide = umis_per_cell_with_guides / guides_per_cell_with_guides
                metrics[f'umis_per_guide_per_cell_cutoff{cutoff}'] = float(np.mean(umis_per_guide))
            else:
                metrics[f'umis_per_guide_per_cell_cutoff{cutoff}'] = np.nan
        else:
            metrics[f'umis_per_guide_per_cell_cutoff{cutoff}'] = np.nan

        # Total guides detected
        metrics[f'total_guides_detected_cutoff{cutoff}'] = n_guides_detected
        
        # Total guide detections across all cells
        metrics[f'total_guide_detections_cutoff{cutoff}'] = total_guide_detections

        # Gene-level metrics if mapping provided
        if guide_to_genes is not None and n_cells_found > 0:
            if 'guide_names' not in guide_adata.uns:
                log_print("Warning: 'guide_names' not found in guide adata.uns, skipping gene metrics")
                metrics[f'total_targeted_genes_cutoff{cutoff}'] = 0
            else:
                guide_names = guide_adata.uns['guide_names']
                unique_genes = set()
                for i, guide_name in enumerate(guide_names):
                    if guides_detected[i] and guide_name in guide_to_genes:
                        unique_genes.update(guide_to_genes[guide_name])
                metrics[f'total_targeted_genes_cutoff{cutoff}'] = len(unique_genes)
        else:
            metrics[f'total_targeted_genes_cutoff{cutoff}'] = 0

    return metrics


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

    # Mark genes used for analysis
    adata.var["hvg_input"] = genes_expressed

    # Calculate HVGs on raw counts (Seurat v3 approach)
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


def perform_dimensionality_reduction(adata, use_gpu=False):
    """Perform PCA and UMAP on data."""
    log_print("📊 Performing dimensionality reduction...")

    # PCA
    # TODO: Make n_comps configurable (currently hardcoded to 50)
    if use_gpu and GPU_AVAILABLE:
        rsc.pp.pca(adata, n_comps=50, use_highly_variable=True)
    else:
        sc.pp.pca(adata, n_comps=50, use_highly_variable=True)

    # Compute neighborhood graph
    # TODO: Make n_neighbors and n_pcs configurable (currently hardcoded to 15)
    if use_gpu and GPU_AVAILABLE:
        rsc.pp.neighbors(adata, n_neighbors=15, n_pcs=15)
    else:
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=15)

    # UMAP
    if use_gpu and GPU_AVAILABLE:
        rsc.tl.umap(adata)
    else:
        sc.tl.umap(adata)

    log_print("✅ Dimensionality reduction completed")

    return adata


def perform_clustering(adata, use_gpu=False, target_clusters=3):
    """Perform Leiden clustering discovering up to 10 clusters, starting at 0.05 resolution."""
    log_print(f"🎯 Performing Leiden clustering (discovering up to 10 clusters)...")

    # TODO: Make clustering parameters configurable (resolution start, max_clusters, limit)
    # Start at 0.05 and increment until we find up to 10 clusters
    resolution = 0.05
    max_clusters = 10
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

    # Assert that leiden_0.05 gives exactly 2 clusters
    if 'leiden_0.05' in adata.obs.columns:
        n_clusters_0_05 = len(adata.obs['leiden_0.05'].unique())
        log_print(f"🔍 Validating leiden_0.05 clustering: found {n_clusters_0_05} clusters")
        if n_clusters_0_05 != 2:
            raise ValueError(f"ASSERTION FAILED: leiden_0.05 clustering produced {n_clusters_0_05} clusters, expected exactly 2 clusters. "
                           f"This indicates the data structure has changed and pseudobulk analysis may not work correctly.")
        log_print("✅ leiden_0.05 validation passed: exactly 2 clusters found")

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
    """Map cells to samples using sample_info.xlsx and sample-to-well mappings.

    This is a placeholder for future implementation when you want to add
    sample annotation functionality similar to statin_perturb.
    """
    raise NotImplementedError("Sample mapping functionality not yet implemented for main pipeline")


