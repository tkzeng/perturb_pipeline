"""
Create analysis-ready GEX and guide files with cell calling filtering.
This script combines sublibrary data using specified cell calling methods.
No nascent RNA support - focuses on total RNA with advanced cell calling.

WARNING: THIS SCRIPT IS IN DEVELOPMENT AND NOT READY FOR PRODUCTION USE
This functionality is still being tested and may not work correctly.
Do not use for production analysis yet.
"""

import anndata as ad
import argparse
import numpy as np
import scanpy as sc
import pandas as pd
import gc
import time
import yaml
from pathlib import Path
import sys
import os

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.pipeline_utils import log_print, process_single_library, process_guide_library


def get_gex_libraries(config, pool_names=None):
    """Get all GEX libraries from sample_info.xlsx
    
    Args:
        config: Pipeline configuration dictionary
        pool_names: Optional list of pool names to filter (None means all pools)
        
    Returns:
        List of GEX sample IDs
    """
    # Load sample information
    sample_info_file = config['sample_info_file']
    sample_df = pd.read_excel(sample_info_file)
    
    # Filter for GEX samples
    gex_samples = sample_df[sample_df['sample_type'] == 'gex']
    
    # Filter by pools if specified
    if pool_names:
        gex_samples = gex_samples[gex_samples['pool'].isin(pool_names)]
    
    if gex_samples.empty:
        pool_msg = f" for pools {pool_names}" if pool_names else ""
        raise ValueError(f"No GEX samples found{pool_msg}")
    
    gex_libs = gex_samples['sample_id'].tolist()
    log_print(f"Found {len(gex_libs)} GEX libraries: {gex_libs}")
    return gex_libs


def get_guide_libraries(config, pool_names=None):
    """Get all guide libraries from sample_info.xlsx
    
    Args:
        config: Pipeline configuration dictionary
        pool_names: Optional list of pool names to filter (None means all pools)
        
    Returns:
        List of guide sample IDs
    """
    # Load sample information
    sample_info_file = config['sample_info_file']
    sample_df = pd.read_excel(sample_info_file)
    
    # Filter for guide samples
    guide_samples = sample_df[sample_df['sample_type'] == 'guide']
    
    # Filter by pools if specified
    if pool_names:
        guide_samples = guide_samples[guide_samples['pool'].isin(pool_names)]
    
    if guide_samples.empty:
        pool_msg = f" for pools {pool_names}" if pool_names else ""
        raise ValueError(f"No guide samples found{pool_msg}")
    
    guide_libs = guide_samples['sample_id'].tolist()
    log_print(f"Found {len(guide_libs)} guide libraries: {guide_libs}")
    return guide_libs


def create_unified_gex_file(config, output_file="gex_all.h5ad", cell_calling_method='Expected_Cells', pool_names=None):
    """Create single combined_gex.h5ad with total counts in X, using configurable cell calling
    
    Args:
        config: Pipeline configuration dictionary
        output_file: Output file path
        cell_calling_method: Which cell calling method to use for filtering
        pool_names: Optional list of pool names to filter (None means all pools)
    """
    
    log_print(f"📦 Creating unified GEX file with cell calling method: {cell_calling_method}")
    
    # Get GEX libraries for specified pools
    gex_libs = get_gex_libraries(config, pool_names)
    
    # Process libraries with cell calling
    log_print("📊 === Processing GEX libraries with cell calling ===")
    adata_list = []
    
    for i, lib in enumerate(gex_libs):
        log_print(f"[{i+1}/{len(gex_libs)}] Processing {lib}...")
        try:
            adata = process_single_library(lib, config, cell_calling_method)
            adata_list.append(adata)
            log_print(f"  ✅ Successfully processed {lib}: {adata.shape[0]} cells")
        except Exception as e:
            log_print(f"  ❌ Failed to process {lib}: {e}")
            raise
    
    if not adata_list:
        raise RuntimeError("No libraries processed successfully")
    
    # Concatenate all libraries
    log_print("🔗 Concatenating GEX libraries...")
    adata_combined = ad.concat(adata_list, join='outer')
    log_print(f"📊 Combined GEX shape: {adata_combined.shape}")
    
    # Log biological sample distribution if available
    if 'biological_sample' in adata_combined.obs.columns:
        sample_counts = adata_combined.obs['biological_sample'].value_counts()
        log_print(f"📊 Combined biological samples: {len(sample_counts)} unique samples")
        for sample, count in sample_counts.head(10).items():
            log_print(f"  {sample}: {count} cells")
        if len(sample_counts) > 10:
            log_print(f"  ... and {len(sample_counts) - 10} more samples")
    
    # Free memory
    del adata_list
    gc.collect()
    
    # Add metadata
    adata_combined.obs['cell_calling_method'] = cell_calling_method
    
    log_print(f"\n📊 Final unified GEX dataset:")
    log_print(f"  Shape: {adata_combined.shape}")
    log_print(f"  Cell calling method: {cell_calling_method}")
    
    log_print(f"\n💾 Writing to {output_file}...")
    write_start = time.time()
    adata_combined.write_h5ad(output_file)
    write_time = time.time() - write_start
    log_print(f"✅ File written in {write_time:.2f} seconds")
    
    log_print("Unified GEX file created successfully!")
    return adata_combined


def create_unified_guide_file(config, gex_file, output_file="guide_all.h5ad", pool_names=None):
    """Create guide file aligned to GEX cells
    
    Args:
        config: Pipeline configuration dictionary
        gex_file: Path to GEX file (for cell filtering)
        output_file: Output guide file path
        pool_names: Optional list of pool names to filter (None means all pools)
    """
    
    log_print("📦 Creating unified guide file aligned to GEX cells")
    
    # Load GEX data to get filtered cells
    log_print(f"Loading GEX data from {gex_file}...")
    adata_gex = sc.read_h5ad(gex_file)
    gex_cells = set(adata_gex.obs.index)
    log_print(f"Found {len(gex_cells)} cells in GEX data for filtering")
    
    # Free GEX memory - we only needed the cell list
    del adata_gex
    gc.collect()
    
    # Get guide libraries
    guide_libs = get_guide_libraries(config, pool_names)
    
    # Process guide libraries
    log_print("📊 === Processing guide libraries ===")
    adata_list = []
    
    for i, lib in enumerate(guide_libs):
        log_print(f"[{i+1}/{len(guide_libs)}] Processing {lib}...")
        try:
            adata = process_guide_library(lib, config, gex_cells)
            adata_list.append(adata)
            log_print(f"  ✅ Successfully processed {lib}: {adata.shape[0]} cells")
        except Exception as e:
            log_print(f"  ❌ Failed to process {lib}: {e}")
            raise
    
    if not adata_list:
        raise RuntimeError("No guide libraries processed successfully")
    
    # Concatenate all libraries
    log_print("🔗 Concatenating guide libraries...")
    adata_combined = ad.concat(adata_list, join='outer')
    log_print(f"📊 Combined guide shape: {adata_combined.shape}")
    
    # Mark all features as guides
    adata_combined.var["type"] = "guide"
    
    # Free memory
    del adata_list
    gc.collect()
    
    log_print(f"\n📊 Final unified guide dataset:")
    log_print(f"  Shape: {adata_combined.shape}")
    log_print(f"  Guide features: {adata_combined.shape[1]} guides")
    
    log_print(f"\n💾 Writing to {output_file}...")
    write_start = time.time()
    adata_combined.write_h5ad(output_file)
    write_time = time.time() - write_start
    log_print(f"✅ File written in {write_time:.2f} seconds")
    
    # Write metadata
    metadata_file = output_file.replace('.h5ad', '.tsv')
    log_print(f"Writing guide metadata to {metadata_file}...")
    adata_combined.obs.to_csv(metadata_file, sep='\t')
    
    log_print("Successfully created guide data!")
    return adata_combined


def main():
    parser = argparse.ArgumentParser(description='Create analysis-ready GEX and guide files')
    parser.add_argument('--config', required=True, help='Path to config.yaml file')
    parser.add_argument('--gex-output', required=True, help='Output GEX file path')
    parser.add_argument('--guide-output', required=True, help='Output guide file path')
    parser.add_argument('--cell-calling-method', default='Expected_Cells',
                        help='Cell calling method to use (default: Expected_Cells)')
    parser.add_argument('--pools', nargs='*', default=None,
                        help='Pool names to process (empty means all pools)')
    args = parser.parse_args()
    
    # Load config
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    # Create GEX file with cell calling
    log_print("=" * 60)
    log_print("STEP 1: Creating unified GEX file")
    log_print("=" * 60)
    adata_gex = create_unified_gex_file(
        config=config,
        output_file=args.gex_output,
        cell_calling_method=args.cell_calling_method,
        pool_names=args.pools
    )
    
    # Create guide file aligned to GEX cells
    log_print("\n" + "=" * 60)
    log_print("STEP 2: Creating unified guide file")
    log_print("=" * 60)
    adata_guide = create_unified_guide_file(
        config=config,
        gex_file=args.gex_output,
        output_file=args.guide_output,
        pool_names=args.pools
    )
    
    log_print("\n" + "=" * 60)
    log_print("✅ Analysis-ready files created successfully!")
    log_print(f"  GEX: {args.gex_output}")
    log_print(f"  Guide: {args.guide_output}")
    log_print("=" * 60)


if __name__ == "__main__":
    main()