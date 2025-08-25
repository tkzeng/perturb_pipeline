#!/usr/bin/env python3
"""
Cell calling analysis using multiple approaches:
1. Expected cell count (top N barcodes)
2. Simple UMI threshold
3. BarcodeRanks (via R script)

This script performs cell calling analysis and saves results.
Plotting is handled separately by cell_calling_plots.py.
"""

import argparse
import subprocess
import sys
import yaml
import pandas as pd
import numpy as np
from pathlib import Path
import scanpy as sc
import os

def load_h5ad_matrix(h5ad_file):
    """Load count matrix from kallisto h5ad output."""
    adata = sc.read_h5ad(h5ad_file)
    
    # Delete guide data from obsm to save memory (not needed for cell calling)
    if 'guide_counts' in adata.obsm:
        print(f"Removing guide_counts matrix from memory ({adata.obsm['guide_counts'].shape})")
        del adata.obsm['guide_counts']
    if 'guide_assignment' in adata.obsm:
        print(f"Removing guide_assignment matrix from memory")
        del adata.obsm['guide_assignment']
    
    # Sum layers to guarantee X is total counts
    adata.X = adata.layers['mature'] + adata.layers['nascent'] + adata.layers['ambiguous']
    
    # Calculate basic QC metrics if not already present
    if 'total_counts' not in adata.obs.columns:
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    return adata

def run_dropletutils_r(kb_dir, sample_id, output_dir, expected_cells=0, lower=100):
    """Run R script for DropletUtils BarcodeRanks analysis."""
    # Get the path to the R script relative to this Python script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path = os.path.join(script_dir, "dropletutils_r_minimal.R")
    
    cmd = [
        "Rscript", "--vanilla", r_script_path,
        "--kb_dir", str(kb_dir),
        "--sample_id", sample_id,
        "--output_dir", str(output_dir),
        "--expected_cells", str(expected_cells),
        "--lower", str(lower)
    ]
    
    result = subprocess.run(cmd, check=True)
    print(f"BarcodeRanks R script completed successfully")

def expected_cell_method(adata, expected_cells):
    """Cell calling using expected cell count - take top N barcodes."""
    if 'total_counts' in adata.obs.columns:
        total_counts = adata.obs['total_counts'].values
    else:
        total_counts = np.array(adata.X.sum(axis=1)).flatten()
    
    sorted_indices = np.argsort(total_counts)[::-1]
    
    # Take top N barcodes
    n_cells = min(expected_cells, len(sorted_indices))
    cell_indices = sorted_indices[:n_cells]
    
    # Create boolean mask
    is_cell = np.zeros(adata.n_obs, dtype=bool)
    is_cell[cell_indices] = True
    
    return is_cell, total_counts[sorted_indices[n_cells-1]] if n_cells > 0 else 0

def threshold_method(adata, min_umi_threshold):
    """Cell calling using simple UMI threshold."""
    if 'total_counts' in adata.obs.columns:
        total_counts = adata.obs['total_counts'].values
    else:
        total_counts = np.array(adata.X.sum(axis=1)).flatten()
    
    is_cell = total_counts >= min_umi_threshold
    return is_cell, min_umi_threshold

def load_dropletutils_results(output_dir, sample_id, adata):
    """Load BarcodeRanks results from minimal R script output."""
    # Load BarcodeRanks results
    results_file = Path(output_dir) / f"{sample_id}_barcoderanks.tsv"
    df = pd.read_csv(results_file, sep='\t')
    
    results = {}
    
    # Get knee and inflection thresholds from the results file
    knee_threshold = df['knee_threshold'].iloc[0]
    inflection_threshold = df['inflection_threshold'].iloc[0]
    
    # Create boolean arrays for cell calling
    is_cell_knee = np.zeros(adata.n_obs, dtype=bool)
    is_cell_inflection = np.zeros(adata.n_obs, dtype=bool)
    
    # Vectorized mapping of barcode results back to adata indices
    # Create a mapping from barcodes to their positions in adata
    adata_barcode_series = pd.Series(range(adata.n_obs), index=adata.obs_names)
    
    # Get indices for barcodes that exist in both df and adata
    common_barcodes = df['barcode'][df['barcode'].isin(adata.obs_names)]
    indices = adata_barcode_series[common_barcodes].values
    
    # Vectorized assignment using the indices
    is_cell_knee[indices] = df.loc[df['barcode'].isin(adata.obs_names), 'is_cell_knee'].values
    is_cell_inflection[indices] = df.loc[df['barcode'].isin(adata.obs_names), 'is_cell_inflection'].values
    
    results['BarcodeRanks_Knee'] = (is_cell_knee, knee_threshold)
    results['BarcodeRanks_Inflection'] = (is_cell_inflection, inflection_threshold)
    
    return results

def calculate_umis_in_cells_pct(adata, methods_results):
    """Calculate percentage of total UMIs that fall in called cells."""
    if 'total_counts' in adata.obs.columns:
        total_counts = adata.obs['total_counts'].values
    else:
        total_counts = np.array(adata.X.sum(axis=1)).flatten()
    
    total_umis = np.sum(total_counts)
    
    percentages = {}
    for method, (is_cell, threshold) in methods_results.items():
        umis_in_cells = np.sum(total_counts[is_cell])
        percentages[method] = (umis_in_cells / total_umis * 100) if total_umis > 0 else 0.0
    
    return percentages

def get_cell_calling_params(config, sample_id):
    """Get cell calling parameters for a specific sample."""
    
    # Load sample-specific parameters from sample_info file (REQUIRED)
    sample_info_file = config['sample_info_file']
    if not os.path.exists(sample_info_file):
        raise FileNotFoundError(f"Sample info file not found: {sample_info_file}")
    
    sample_df = pd.read_csv(sample_info_file, sep='\t')
    
    sample_row = sample_df[sample_df['sample_id'] == sample_id]
    if sample_row.empty:
        raise ValueError(f"Sample {sample_id} not found in {sample_info_file}")
    
    row = sample_row.iloc[0]
    
    # REQUIRED: expected_cells must be specified
    if 'expected_cells' not in sample_df.columns or pd.isna(row['expected_cells']):
        raise ValueError(f"Expected cells must be specified for sample {sample_id} in {sample_info_file}")
    
    params = {'expected_cells': int(row['expected_cells'])}
    
    # Get min_umi_threshold from sample_info if provided
    if pd.notna(row.get('min_umi_threshold')):
        params['min_umi_threshold'] = int(row['min_umi_threshold'])
    
    
    return params

def main():
    parser = argparse.ArgumentParser(description="Cell calling analysis")
    parser.add_argument("--h5ad_file", required=True, help="Kallisto h5ad output file")
    parser.add_argument("--kb_dir", required=True, help="Kallisto bus output directory")
    parser.add_argument("--sample-id", required=True, help="Full sample ID (format: pool:sample)")
    parser.add_argument("--config", required=True, help="Config YAML file")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--ncores", type=int, required=True, help="Number of cores for parallel processing")
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)
    
    # Get sample-specific parameters
    params = get_cell_calling_params(config, args.sample_id)
    
    # Load h5ad matrix
    print(f"Loading h5ad file: {args.h5ad_file}")
    adata = load_h5ad_matrix(args.h5ad_file)
    
    print(f"Loaded matrix: {adata.n_obs} barcodes x {adata.n_vars} genes")
    
    # Apply different cell calling methods
    methods_results = {}
    
    # Get list of methods to run from config
    methods_to_run = config['cell_calling']['methods_to_run']
    print(f"Running cell calling methods: {methods_to_run}")
    
    # Method 1: Expected cell count
    if 'Expected_Cells' in methods_to_run:
        # Require expected_cells in params
        if 'expected_cells' not in params:
            raise ValueError("Expected_Cells method requires 'expected_cells' in sample_info.tsv")
        print(f"Running expected cell method (n={params['expected_cells']})...")
        methods_results['Expected_Cells'] = expected_cell_method(adata, params['expected_cells'])
    
    # Method 2: Threshold method
    if 'UMI_Threshold' in methods_to_run:
        # Require min_umi_threshold in params
        if 'min_umi_threshold' not in params:
            raise ValueError("UMI_Threshold method requires 'min_umi_threshold' in sample_info.tsv")
        print(f"Running threshold method (min_umi={params['min_umi_threshold']})...")
        methods_results['UMI_Threshold'] = threshold_method(adata, params['min_umi_threshold'])
    
    # Method 3-4: BarcodeRanks methods (via R)
    barcoderanks_methods = ['BarcodeRanks_Knee', 'BarcodeRanks_Inflection']
    run_barcoderanks = any(method in methods_to_run for method in barcoderanks_methods)
    
    if run_barcoderanks:
        print("Running BarcodeRanks analysis on pre-filtered matrix...")
        print("  - BarcodeRanks enabled")
            
        run_dropletutils_r(args.kb_dir, args.sample_id, output_dir, 
                          params['expected_cells'], 
                          100)  # lower threshold for BarcodeRanks
        
        dropletutils_results = load_dropletutils_results(output_dir, args.sample_id, adata)
        # Only include methods that were requested
        for method in methods_to_run:
            if method in dropletutils_results:
                methods_results[method] = dropletutils_results[method]
    
    # Calculate percentage of UMIs in cells
    umi_percentages = calculate_umis_in_cells_pct(adata, methods_results)
    
    # Save results
    results_data = []
    for method, (is_cell, threshold) in methods_results.items():
        n_cells = np.sum(is_cell)
        umis_in_cells_pct = umi_percentages[method]
        
        results_data.append({
            'sample_id': args.sample_id,
            'method': method,
            'n_cells_called': n_cells,
            'threshold_used': threshold,
            'umis_in_cells_pct': umis_in_cells_pct,
            'expected_cells_param': params['expected_cells'],
            'min_umi_threshold_param': params['min_umi_threshold']
        })
    
    # Save to TSV
    results_df = pd.DataFrame(results_data)
    results_df.to_csv(output_dir / 'results.tsv', 
                      sep='\t', index=False)
    
    # Save cell barcodes for each method
    for method, (is_cell, threshold) in methods_results.items():
        cell_barcodes = adata.obs_names[is_cell]
        with open(output_dir / f'{args.sample_id}_{method}_cell_barcodes.txt', 'w') as f:
            for barcode in cell_barcodes:
                f.write(f"{barcode}\n")
    
    print(f"Cell calling analysis completed. Results saved to {output_dir}")
    
    # Print summary
    print("\nSummary:")
    for method, (is_cell, threshold) in methods_results.items():
        n_cells = np.sum(is_cell)
        umis_pct = umi_percentages[method]
        print(f"  {method}: {n_cells} cells (threshold: {threshold:.0f}, "
              f"UMIs in cells: {umis_pct:.1f}%)")
    
    # Clean up memory
    del adata
    import gc
    gc.collect()

if __name__ == "__main__":
    main()