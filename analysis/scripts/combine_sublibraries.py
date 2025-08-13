"""
Combine cell-called sublibraries into a single analysis-ready file.
Clean implementation with explicit file paths and no hardcoded fallbacks.

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
from pathlib import Path
import sys
import os
import yaml
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.pipeline_utils import log_print


def log_memory_usage():
    """Log current memory usage if psutil is available"""
    if PSUTIL_AVAILABLE:
        process = psutil.Process(os.getpid())
        mem_gb = process.memory_info().rss / 1024 / 1024 / 1024
        log_print(f"  Current memory usage: {mem_gb:.2f} GB")
    else:
        log_print("  Memory monitoring not available (psutil not installed)")


def get_sample_filter_params(sample_id, default_min_umi, default_mito_cutoff, per_sample_filters=None):
    """Get filtering parameters for a specific sample
    
    Args:
        sample_id: Sample identifier
        default_min_umi: Default minimum UMI threshold
        default_mito_cutoff: Default mitochondrial cutoff
        per_sample_filters: Dict mapping sample_id to filter params (optional)
        
    Returns:
        Dict with 'min_umi_filter' and 'mito_cutoff' keys
    """
    if per_sample_filters and sample_id in per_sample_filters:
        sample_params = per_sample_filters[sample_id]
        return {
            'min_umi_filter': sample_params.get('min_umi_filter', default_min_umi),
            'mito_cutoff': sample_params.get('mito_cutoff', default_mito_cutoff)
        }
    else:
        return {
            'min_umi_filter': default_min_umi,
            'mito_cutoff': default_mito_cutoff
        }


def filter_files_by_pools(h5ad_files, barcode_files, pools):
    """Filter file lists by pool names if specified
    
    Args:
        h5ad_files: List of h5ad file paths
        barcode_files: List of barcode file paths
        pools: Space-separated string of pool names to keep (empty for all)
    
    Returns:
        Tuple of (filtered_h5ad_files, filtered_barcode_files)
    """
    if not pools or pools.strip() == "":
        return h5ad_files, barcode_files
    
    pool_list = pools.strip().split()
    log_print(f"Filtering to pools: {pool_list}")
    
    # Extract sample IDs and filter by pools that appear in paths
    filtered_h5ad = []
    filtered_barcode = []
    
    for h5ad_path, barcode_path in zip(h5ad_files, barcode_files):
        # Check if any pool name appears in the path
        if any(pool in str(h5ad_path) for pool in pool_list):
            filtered_h5ad.append(h5ad_path)
            filtered_barcode.append(barcode_path)
    
    log_print(f"Filtered to {len(filtered_h5ad)} files from specified pools")
    return filtered_h5ad, filtered_barcode


def load_and_filter_library(h5ad_path, barcode_path, cell_calling_method, filter_params):
    """Load an h5ad file and filter it to called cells
    
    Args:
        h5ad_path: Path to annotated h5ad file
        barcode_path: Path to cell barcode file
        cell_calling_method: Cell calling method used for filtering
        filter_params: Dict with 'min_umi_filter' and 'mito_cutoff' keys
        
    Returns:
        Filtered AnnData object with library-prefixed cell names
    """
    # Load annotated h5ad file
    adata = sc.read_h5ad(h5ad_path)
    original_cells = adata.shape[0]
    
    # Get sample ID from the h5ad metadata (added by sublibrary_annotation.py)
    sample_id = adata.obs['sample_id'].iloc[0] if 'sample_id' in adata.obs.columns else "unknown"
    
    log_print(f"Loading {sample_id}...")
    
    # Load called cell barcodes
    with open(barcode_path) as f:
        called_cells = set(line.strip() for line in f)
    
    # Filter to called cells
    adata = adata[adata.obs_names.isin(called_cells)]
    log_print(f"  Cell calling filter {sample_id}: {adata.shape[0]} cells (from {original_cells} total)")
    
    # Apply UMI threshold filter
    min_umi_filter = filter_params['min_umi_filter']
    if min_umi_filter > 0:
        pre_umi_cells = adata.shape[0]
        adata = adata[adata.obs['total_counts'] >= min_umi_filter]
        log_print(f"  UMI filter {sample_id}: {adata.shape[0]} cells (removed {pre_umi_cells - adata.shape[0]} with <{min_umi_filter} UMIs)")
    
    # Apply mitochondrial percentage filter
    mito_cutoff = filter_params['mito_cutoff']
    if mito_cutoff > 0 and 'pct_counts_mt' in adata.obs.columns:
        pre_mito_cells = adata.shape[0]
        adata = adata[adata.obs['pct_counts_mt'] < mito_cutoff]
        log_print(f"  Mito filter {sample_id}: {adata.shape[0]} cells (removed {pre_mito_cells - adata.shape[0]} with ≥{mito_cutoff}% mito)")
    
    # Add library prefix to cell names
    adata.obs_names = [f"{sample_id}_{barcode}" for barcode in adata.obs_names]
    
    # Add cell calling method metadata
    adata.obs["cell_calling_method"] = cell_calling_method
    
    # Note: Library metadata (sample_id, pool) already added by sublibrary_annotation.py
    
    return adata


# Removed extract_guide_data function - no longer needed
# Output is a single combined h5ad file with both GEX and guide data


def combine_sublibraries(h5ad_files, barcode_files, filter_files, cell_calling_method, output_file, pools=""):
    """Main function to combine libraries into a single file
    
    Args:
        h5ad_files: List of annotated h5ad file paths
        barcode_files: List of cell barcode file paths
        filter_files: List of per-sample filter parameter YAML files
        cell_calling_method: Cell calling method used
        output_file: Output path for combined h5ad file
        pools: Space-separated pool names to filter (empty for all)
    """
    
    log_print("🧬 Combining filtered sublibraries")
    log_print(f"Cell calling method: {cell_calling_method}")
    log_print(f"Input files: {len(h5ad_files)} h5ad files")
    log_print(f"Filter files: {len(filter_files)} per-sample YAML files")
    
    # Filter by pools if specified
    h5ad_files, barcode_files = filter_files_by_pools(h5ad_files, barcode_files, pools)
    
    # Step 1: Load, filter, and incrementally concatenate libraries
    log_print("\n📥 STEP 1: Loading and filtering libraries with memory optimization...")
    log_memory_usage()
    
    combined_gex = None
    
    for i, (h5ad_path, barcode_path, filter_file) in enumerate(zip(h5ad_files, barcode_files, filter_files)):
        # Load filtering parameters from YAML file
        with open(filter_file) as f:
            filter_params = yaml.safe_load(f)
        
        # Load and filter current library
        adata = load_and_filter_library(h5ad_path, barcode_path, cell_calling_method, filter_params)
        
        if combined_gex is None:
            # First library becomes the base
            combined_gex = adata
            log_print(f"  Initialized with library 1/{len(h5ad_files)}: {adata.shape[0]} cells")
        else:
            # Concatenate with accumulated result
            log_print(f"  Concatenating library {i+1}/{len(h5ad_files)}: {adata.shape[0]} cells")
            combined_gex = ad.concat([combined_gex, adata], join='outer')
            
            # Immediately free memory from individual library
            del adata
            gc.collect()
            
            log_print(f"  Running total: {combined_gex.shape[0]} cells")
        
        # Log memory usage every few libraries
        if (i + 1) % 3 == 0 or i == len(h5ad_files) - 1:
            log_memory_usage()
    
    log_print(f"\n🔗 STEP 2: Final combined shape: {combined_gex.shape}")
    
    # Log biological sample distribution if available
    if 'biological_sample' in combined_gex.obs.columns:
        sample_counts = combined_gex.obs['biological_sample'].value_counts()
        log_print(f"📊 Biological samples: {len(sample_counts)} unique samples")
        for sample, count in sample_counts.head(10).items():
            log_print(f"  {sample}: {count} cells")
        if len(sample_counts) > 10:
            log_print(f"  ... and {len(sample_counts) - 10} more samples")
    
    # Step 3: Save combined file
    log_print("\n💾 STEP 3: Saving combined output...")
    log_memory_usage()
    
    log_print(f"Writing combined data to {output_file}...")
    write_start = time.time()
    combined_gex.write_h5ad(output_file)
    write_time = time.time() - write_start
    log_print(f"✅ File written in {write_time:.2f} seconds")
    
    log_print("\n✅ Combined sublibrary file created successfully!")
    log_print(f"  Combined: {output_file} ({combined_gex.shape[0]} cells × {combined_gex.shape[1]} genes)")
    log_print(f"  Guide data: Available in obsm['guide_counts'] ({len(combined_gex.uns.get('guide_names', []))} guides)")


def main():
    parser = argparse.ArgumentParser(description='Combine sublibraries into a single file')
    parser.add_argument('--h5ad-files', required=True, nargs='+', help='List of annotated h5ad file paths')
    parser.add_argument('--barcode-files', required=True, nargs='+', help='List of cell barcode file paths')
    parser.add_argument('--cell-calling-method', required=True, help='Cell calling method used')
    parser.add_argument('--output', required=True, help='Output combined h5ad file path')
    parser.add_argument('--pools', default='', help='Space-separated pool names to process (empty means all pools)')
    parser.add_argument('--min-umi-filter', type=int, default=100, help='Minimum UMI count threshold')
    parser.add_argument('--mito-cutoff', type=float, default=20, help='Maximum mitochondrial percentage threshold')
    parser.add_argument('--per-sample-filters', help='YAML file with per-sample filter parameters')
    args = parser.parse_args()
    
    # Load per-sample filters if provided
    per_sample_filters = None
    if args.per_sample_filters:
        with open(args.per_sample_filters) as f:
            per_sample_filters = yaml.safe_load(f)
    
    # Run the combination
    combine_sublibraries(
        h5ad_files=args.h5ad_files,
        barcode_files=args.barcode_files,
        cell_calling_method=args.cell_calling_method,
        output_file=args.output,
        pools=args.pools,
        min_umi_filter=args.min_umi_filter,
        mito_cutoff=args.mito_cutoff,
        per_sample_filters=per_sample_filters
    )


if __name__ == "__main__":
    main()