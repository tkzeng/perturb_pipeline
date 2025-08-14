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




def load_and_filter_library(h5ad_path, barcode_path, cell_calling_method, cell_lists_path, qc_method='gmm_2d_posterior_75'):
    """Load an h5ad file and filter it to called cells that pass QC
    
    Args:
        h5ad_path: Path to annotated h5ad file
        barcode_path: Path to cell barcode file  
        cell_calling_method: Cell calling method used for filtering
        cell_lists_path: Path to TSV file with QC cell lists
        qc_method: Column name in cell lists to use for filtering
        
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
    
    # Apply QC filtering using cell lists
    if not os.path.exists(cell_lists_path):
        raise FileNotFoundError(f"Cell lists file not found: {cell_lists_path}")
    
    cell_lists_df = pd.read_csv(cell_lists_path, sep='\t')
    
    if qc_method not in cell_lists_df.columns:
        available_methods = [col for col in cell_lists_df.columns if col != 'barcode']
        raise ValueError(f"QC method '{qc_method}' not found in cell lists file. "
                        f"Available methods: {', '.join(available_methods)}")
    
    # Get cells that pass QC
    qc_pass = cell_lists_df[cell_lists_df[qc_method] == True]['barcode'].tolist()
    pre_qc_cells = adata.shape[0]
    adata = adata[adata.obs_names.isin(qc_pass)]
    log_print(f"  QC filter {sample_id}: {adata.shape[0]} cells (removed {pre_qc_cells - adata.shape[0]} using {qc_method})")
    
    # Add library prefix to cell names
    adata.obs_names = [f"{sample_id}_{barcode}" for barcode in adata.obs_names]
    
    # Note: Library metadata (sample_id, pool) already added by sublibrary_annotation.py
    
    return adata


# Removed extract_guide_data function - no longer needed
# Output is a single combined h5ad file with both GEX and guide data


def combine_sublibraries(h5ad_files, barcode_files, filter_files, cell_calling_method, output_file, qc_method='gmm_2d_posterior_75'):
    """Main function to combine libraries into a single file
    
    Args:
        h5ad_files: List of annotated h5ad file paths
        barcode_files: List of cell barcode file paths
        filter_files: List of per-sample QC cell list TSV files
        cell_calling_method: Cell calling method used
        output_file: Output path for combined h5ad file
        qc_method: QC method column to use for filtering
    """
    
    log_print("🧬 Combining filtered sublibraries")
    log_print(f"Cell calling method: {cell_calling_method}")
    log_print(f"QC filtering method: {qc_method}")
    log_print(f"Input files: {len(h5ad_files)} h5ad files")
    log_print(f"Filter files: {len(filter_files)} per-sample cell list files")
    
    
    # Step 1: Load, filter, and incrementally concatenate libraries
    log_print("\n📥 STEP 1: Loading and filtering libraries with memory optimization...")
    log_memory_usage()
    
    combined_gex = None
    
    for i, (h5ad_path, barcode_path, filter_file) in enumerate(zip(h5ad_files, barcode_files, filter_files)):
        # Load and filter current library
        adata = load_and_filter_library(h5ad_path, barcode_path, cell_calling_method, filter_file, qc_method)
        
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
    
    # Add filtering metadata to uns (dataset-level metadata)
    combined_gex.uns['filtering_metadata'] = {
        'cell_calling_method': cell_calling_method,
        'qc_filtering_method': qc_method
    }
    log_print(f"📝 Added filtering metadata: cell_calling={cell_calling_method}, qc_method={qc_method}")
    
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
    log_print(f"  Combined: {output_file} ({combined_gex.shape[0]} cells x {combined_gex.shape[1]} genes)")
    log_print(f"  Guide data: Available in obsm['guide_counts'] ({len(combined_gex.uns.get('guide_names', []))} guides)")


def main():
    parser = argparse.ArgumentParser(description='Combine sublibraries into a single file')
    parser.add_argument('--h5ad-files', required=True, nargs='+', help='List of annotated h5ad file paths')
    parser.add_argument('--barcode-files', required=True, nargs='+', help='List of cell barcode file paths')
    parser.add_argument('--filter-files', required=True, nargs='+', help='List of QC cutoffs YAML file paths')
    parser.add_argument('--cell-calling-method', required=True, help='Cell calling method used')
    parser.add_argument('--output', required=True, help='Output combined h5ad file path')
    parser.add_argument('--qc-method', default='gmm_2d_posterior_75', 
                       help='QC method column to use for filtering (default: gmm_2d_posterior_75)')
    args = parser.parse_args()
    
    
    # Run the combination
    combine_sublibraries(
        h5ad_files=args.h5ad_files,
        barcode_files=args.barcode_files,
        filter_files=args.filter_files,
        cell_calling_method=args.cell_calling_method,
        output_file=args.output,
        qc_method=args.qc_method
    )


if __name__ == "__main__":
    main()