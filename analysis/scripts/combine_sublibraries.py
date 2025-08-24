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




def load_and_filter_library(h5ad_path, barcode_path, cell_calling_method, cell_lists_path):
    """Load an h5ad file and filter it to called cells, adding GMM probabilities if available
    
    Args:
        h5ad_path: Path to annotated h5ad file
        barcode_path: Path to cell barcode file  
        cell_calling_method: Cell calling method used for filtering
        cell_lists_path: Path to TSV file with GMM posterior probabilities
        
    Returns:
        AnnData object with called cells and GMM probabilities (if available)
    """
    # Load annotated h5ad file
    adata = sc.read_h5ad(h5ad_path)
    original_cells = adata.shape[0]
    
    # Get sample ID from the h5ad metadata (added by sublibrary_annotation.py)
    sample_id = adata.obs['sample_id'].iloc[0] if 'sample_id' in adata.obs.columns else "unknown"
    
    log_print(f"Loading {sample_id}...")
    
    # Delete layers immediately - they're not used in downstream analysis
    if adata.layers:
        layer_names = list(adata.layers.keys())
        log_print(f"  Deleting {len(layer_names)} layers to save memory: {layer_names}")
        for layer in layer_names:
            del adata.layers[layer]
    
    # Load called cell barcodes
    with open(barcode_path) as f:
        called_cells = set(line.strip() for line in f)
    
    # Filter to called cells
    adata = adata[adata.obs_names.isin(called_cells)]
    log_print(f"  Cell calling filter {sample_id}: {adata.shape[0]} cells (from {original_cells} total)")
    
    # Add cell quality posterior probabilities (always present after generate_qc_cell_lists runs)
    if os.path.exists(cell_lists_path):
        cell_lists_df = pd.read_csv(cell_lists_path, sep='\t')
        
        # Match barcodes and add the probability column
        cell_lists_df = cell_lists_df.set_index('barcode')
        
        # Add cell quality probabilities for cells that exist in both datasets
        common_barcodes = adata.obs_names.intersection(cell_lists_df.index)
        
        # Initialize column with NaN for all cells
        adata.obs['posterior_prob_compromised'] = np.nan
        
        # Fill in values for cells that have cell quality probabilities
        adata.obs.loc[common_barcodes, 'posterior_prob_compromised'] = \
            cell_lists_df.loc[common_barcodes, 'posterior_prob_compromised'].values
        
        n_with_gmm = (~adata.obs['posterior_prob_compromised'].isna()).sum()
        log_print(f"  Added cell quality posterior probabilities for {n_with_gmm}/{adata.shape[0]} cells")
    else:
        log_print(f"  Warning: Cell lists file not found: {cell_lists_path}")
    
    # Add library prefix to cell names using vectorized string operations
    adata.obs_names = pd.Index(adata.obs_names).astype(str).map(lambda x: f"{sample_id}_{x}")
    
    # Note: Library metadata (sample_id, pool) already added by sublibrary_annotation.py
    
    return adata


# Removed extract_guide_data function - no longer needed
# Output is a single combined h5ad file with both GEX and guide data


def combine_sublibraries(h5ad_files, barcode_files, filter_files, threshold_tables, cell_calling_method, output_file, config):
    """Main function to combine libraries into a single file
    
    Args:
        h5ad_files: List of annotated h5ad file paths
        barcode_files: List of cell barcode file paths
        filter_files: List of per-sample QC cell list TSV files (with GMM probabilities)
        cell_calling_method: Cell calling method used
        output_file: Output path for combined h5ad file
        config: Configuration dictionary
    """
    
    log_print("🧬 Combining filtered sublibraries")
    log_print(f"Cell calling method: {cell_calling_method}")
    log_print(f"Input files: {len(h5ad_files)} h5ad files")
    log_print(f"GMM probability files: {len(filter_files)} per-sample cell list files")
    
    
    # Step 1: Load, filter, and incrementally concatenate libraries
    log_print("\n📥 STEP 1: Loading and filtering libraries with memory optimization...")
    log_memory_usage()
    
    combined_gex = None
    
    for i, (h5ad_path, barcode_path, filter_file) in enumerate(zip(h5ad_files, barcode_files, filter_files)):
        # Load and filter current library
        adata = load_and_filter_library(h5ad_path, barcode_path, cell_calling_method, filter_file)
        
        if combined_gex is None:
            # First library becomes the base
            combined_gex = adata
            log_print(f"  Initialized with library 1/{len(h5ad_files)}: {adata.shape[0]} cells")
        else:
            # Verify that the files have matching structure before concatenation
            assert combined_gex.n_vars == adata.n_vars, \
                f"Gene count mismatch: {combined_gex.n_vars} vs {adata.n_vars}"
            assert combined_gex.var_names.equals(adata.var_names), \
                "Gene names don't match between files"
            # Layers check removed - we delete layers to save memory
            assert set(combined_gex.obsm.keys()) == set(adata.obsm.keys()), \
                f"Obsm keys mismatch: {set(combined_gex.obsm.keys())} vs {set(adata.obsm.keys())}"
            
            # Verify guide data consistency
            if 'guide_counts' in combined_gex.obsm:
                assert combined_gex.obsm['guide_counts'].shape[1] == adata.obsm['guide_counts'].shape[1], \
                    f"Guide count dimensions mismatch: {combined_gex.obsm['guide_counts'].shape[1]} vs {adata.obsm['guide_counts'].shape[1]}"
            
            if 'guide_names' in combined_gex.uns:
                import numpy as np
                assert np.array_equal(combined_gex.uns['guide_names'], adata.uns['guide_names']), \
                    "Guide names don't match or are in different order between files"
            
            # Concatenate with accumulated result
            log_print(f"  Concatenating library {i+1}/{len(h5ad_files)}: {adata.shape[0]} cells")
            # Store expected dimensions before concat
            expected_n_cells = combined_gex.n_obs + adata.n_obs
            expected_n_guides = combined_gex.obsm['guide_counts'].shape[1] if 'guide_counts' in combined_gex.obsm else 0
            
            # Use join='inner' since all files have identical genes
            # Use merge='same' to only keep var annotations that are identical across files
            # Use uns_merge='same' to preserve uns data that's identical (guide_names)
            combined_gex = ad.concat([combined_gex, adata], join='inner', merge='same', uns_merge='same')
            
            # Verify critical data is preserved after concatenation
            assert combined_gex.n_obs == expected_n_cells, \
                f"Cell count mismatch after concat: {combined_gex.n_obs} vs expected {expected_n_cells}"
            assert combined_gex.n_vars == adata.n_vars, \
                f"Gene count changed after concat: {combined_gex.n_vars} vs {adata.n_vars}"
            assert 'gene_type' in combined_gex.var.columns, \
                "Critical annotation 'gene_type' lost during concatenation"
            assert 'guide_counts' in combined_gex.obsm, \
                "Guide counts lost from obsm during concatenation"
            assert combined_gex.obsm['guide_counts'].shape[1] == expected_n_guides, \
                f"Guide dimension changed: {combined_gex.obsm['guide_counts'].shape[1]} vs expected {expected_n_guides}"
            if 'guide_names' in adata.uns:
                assert 'guide_names' in combined_gex.uns, \
                    "Guide names lost from uns during concatenation"
                assert len(combined_gex.uns['guide_names']) == len(adata.uns['guide_names']), \
                    f"Guide names dimension changed: {len(combined_gex.uns['guide_names'])} vs {len(adata.uns['guide_names'])}"
            
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
        'gmm_probabilities_added': True
    }
    n_with_prob = (~combined_gex.obs['posterior_prob_compromised'].isna()).sum()
    log_print(f"📝 Added filtering metadata: cell_calling={cell_calling_method}, GMM probabilities for {n_with_prob}/{combined_gex.n_obs} cells")
    
    # Log biological sample distribution if available
    if 'biological_sample' in combined_gex.obs.columns:
        sample_counts = combined_gex.obs['biological_sample'].value_counts()
        log_print(f"📊 Biological samples: {len(sample_counts)} unique samples")
        for sample, count in sample_counts.head(10).items():
            log_print(f"  {sample}: {count} cells")
        if len(sample_counts) > 10:
            log_print(f"  ... and {len(sample_counts) - 10} more samples")
    
    # Step 3: Add GMM threshold tables if provided
    if threshold_tables:
        log_print("\n📊 STEP 3: Adding GMM threshold tables...")
        
        # Read and combine all threshold tables
        threshold_dfs = []
        for threshold_file in threshold_tables:
            if os.path.exists(threshold_file):
                log_print(f"  Reading threshold table: {threshold_file}")
                df = pd.read_csv(threshold_file, sep='\t', index_col=0)
                threshold_dfs.append(df)
            else:
                log_print(f"  WARNING: Threshold table not found: {threshold_file}")
        
        if threshold_dfs:
            # Combine all tables (they should have the same cells in the same order)
            combined_thresholds = pd.concat(threshold_dfs, axis=0)
            
            # Remove duplicate indices (if same cells appear in multiple samples)
            combined_thresholds = combined_thresholds[~combined_thresholds.index.duplicated(keep='first')]
            
            # Align with combined data
            common_cells = combined_gex.obs.index.intersection(combined_thresholds.index)
            log_print(f"  Found thresholds for {len(common_cells)}/{combined_gex.n_obs} cells")
            
            # Get config values - crash if not present
            default_guide_cutoff = config['qc_analysis']['default_guide_cutoff']
            default_granularity = config['qc_analysis']['default_granularity']
            
            # Build the expected column name
            if isinstance(default_guide_cutoff, str) and default_guide_cutoff.startswith('gmm_'):
                level = default_guide_cutoff[4:]  # Extract number from "gmm_50"
                target_col = f'gmm_{level}_{default_granularity}_{cell_calling_method}'
                
                # Column must exist - crash if not
                combined_gex.obs['guide_threshold'] = -1  # Default value for cells without data
                combined_gex.obs.loc[common_cells, 'guide_threshold'] = combined_thresholds.loc[common_cells, target_col]
                n_valid = (combined_gex.obs['guide_threshold'] > 0).sum()
                log_print(f"  Added guide_threshold column (from {target_col}): {n_valid} cells with valid thresholds")
                
                # Store the guide calling method information in uns
                combined_gex.uns['guide_calling_metadata'] = {
                    'method': default_guide_cutoff,
                    'granularity': default_granularity,
                    'source_column': target_col,
                    'n_cells_with_thresholds': int(n_valid)
                }
            else:
                # Fixed threshold - no column needed, store value in uns
                combined_gex.uns['guide_calling_metadata'] = {
                    'method': 'fixed',
                    'threshold': int(default_guide_cutoff),
                    'granularity': 'uniform'
                }
    
    # Step 4: Save combined file
    log_print("\n💾 STEP 4: Saving combined output...")
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
    parser.add_argument('--threshold-tables', required=False, nargs='+', help='List of GMM threshold table paths')
    parser.add_argument('--cell-calling-method', required=True, help='Cell calling method used')
    parser.add_argument('--config', required=True, help='Path to config file')
    parser.add_argument('--output', required=True, help='Output combined h5ad file path')
    args = parser.parse_args()
    
    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)
    
    # Run the combination
    combine_sublibraries(
        h5ad_files=args.h5ad_files,
        barcode_files=args.barcode_files,
        filter_files=args.filter_files,
        threshold_tables=args.threshold_tables,
        cell_calling_method=args.cell_calling_method,
        output_file=args.output,
        config=config
    )


if __name__ == "__main__":
    main()