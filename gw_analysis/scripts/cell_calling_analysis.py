#!/usr/bin/env python3
"""
Cell calling analysis using multiple approaches:
1. Expected cell count (top N barcodes)
2. Simple UMI threshold
3. EmptyDrops (via R script)

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

def run_emptydrops_r(kb_dir, sample_id, output_dir, emptydrops_lower=100, ncores=4):
    """Run R script for emptyDrops analysis."""
    cmd = [
        "Rscript", "scripts/emptydrops_r.R",
        "--kb_dir", str(kb_dir),
        "--sample_id", sample_id,
        "--output_dir", str(output_dir),
        "--emptydrops_lower", str(emptydrops_lower),
        "--ncores", str(ncores)
    ]
    
    result = subprocess.run(cmd, check=True)
    print(f"EmptyDrops R script completed successfully")

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
    """Load DropletUtils results from R script output (EmptyDrops + barcodeRanks)."""
    # Load DropletUtils results
    results_file = Path(output_dir) / "emptydrops_results.tsv"
    df = pd.read_csv(results_file, sep='\t')
    
    # Map back to adata
    barcode_to_idx = {bc: i for i, bc in enumerate(adata.obs_names)}
    
    # Load results for all methods
    results = {}
    
    # 1. EmptyDrops with different FDR thresholds
    for fdr_threshold, suffix in [(0.001, 'FDR0001'), (0.01, 'FDR001'), (0.05, 'FDR005')]:
        is_cell = np.zeros(adata.n_obs, dtype=bool)
        
        # Get cells below FDR threshold
        cells_df = df[df['fdr'] <= fdr_threshold]
        for barcode in cells_df['barcode']:
            if barcode in barcode_to_idx:
                is_cell[barcode_to_idx[barcode]] = True
        
        results[f'EmptyDrops_{suffix}'] = (is_cell, fdr_threshold)
    
    # 2. BarcodeRanks knee and inflection points
    # The R script saves a summary file with these values
    summary_file = Path(output_dir) / "emptydrops_summary.tsv"
    if summary_file.exists():
        summary_df = pd.read_csv(summary_file, sep='\t')
        
        # Get knee and inflection thresholds
        knee_row = summary_df[summary_df['method'] == 'BarcodeRanks_Knee']
        inflection_row = summary_df[summary_df['method'] == 'BarcodeRanks_Inflection']
        
        if not knee_row.empty:
            knee_threshold = knee_row.iloc[0]['threshold_used']
            # Get cells above knee threshold
            if 'total_counts' in adata.obs.columns:
                total_counts = adata.obs['total_counts'].values
            else:
                total_counts = np.array(adata.X.sum(axis=1)).flatten()
            
            is_cell_knee = total_counts >= knee_threshold
            results['BarcodeRanks_Knee'] = (is_cell_knee, knee_threshold)
        
        if not inflection_row.empty:
            inflection_threshold = inflection_row.iloc[0]['threshold_used']
            # Get cells above inflection threshold
            if 'total_counts' in adata.obs.columns:
                total_counts = adata.obs['total_counts'].values
            else:
                total_counts = np.array(adata.X.sum(axis=1)).flatten()
            
            is_cell_inflection = total_counts >= inflection_threshold
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
    import os
    
    # Load sample-specific parameters from sample_info file (REQUIRED)
    sample_info_file = config.get('sample_info_file', 'sample_info.xlsx')
    if not os.path.exists(sample_info_file):
        raise FileNotFoundError(f"Sample info file not found: {sample_info_file}")
    
    if sample_info_file.endswith('.xlsx'):
        sample_df = pd.read_excel(sample_info_file)
    elif sample_info_file.endswith('.tsv'):
        sample_df = pd.read_csv(sample_info_file, sep='\t')
    else:
        raise ValueError(f"Unsupported file format: {sample_info_file}. Only .xlsx and .tsv files are supported.")
    
    sample_row = sample_df[sample_df['sample_id'] == sample_id]
    if sample_row.empty:
        raise ValueError(f"Sample {sample_id} not found in {sample_info_file}")
    
    row = sample_row.iloc[0]
    
    # REQUIRED: expected_cells must be specified
    if 'expected_cells' not in sample_df.columns or pd.isna(row['expected_cells']):
        raise ValueError(f"Expected cells must be specified for sample {sample_id} in {sample_info_file}")
    
    params = {'expected_cells': int(row['expected_cells'])}
    
    # Optional parameters with config defaults
    defaults = config.get('cell_calling', {}).get('defaults', {})
    if pd.notna(row.get('min_umi_threshold')):
        params['min_umi_threshold'] = int(row['min_umi_threshold'])
    else:
        params['min_umi_threshold'] = defaults.get('min_umi_threshold', 100)
    
    params['emptydrops_lower'] = defaults.get('emptydrops_lower', 100)
    
    return params

def main():
    parser = argparse.ArgumentParser(description="Cell calling analysis")
    parser.add_argument("--h5ad_file", required=True, help="Kallisto h5ad output file")
    parser.add_argument("--kb_dir", required=True, help="Kallisto bus output directory")
    parser.add_argument("--sample-id", required=True, help="Full sample ID (format: pool:sample)")
    parser.add_argument("--config", required=True, help="Config YAML file")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--ncores", type=int, default=4, help="Number of cores for parallel processing")
    
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
    
    # Method 1: Expected cell count
    print(f"Running expected cell method (n={params['expected_cells']})...")
    methods_results['Expected_Cells'] = expected_cell_method(adata, params['expected_cells'])
    
    # Method 2: Threshold method
    print(f"Running threshold method (min_umi={params['min_umi_threshold']})...")
    methods_results['UMI_Threshold'] = threshold_method(adata, params['min_umi_threshold'])
    
    # Method 3-7: DropletUtils methods (via R)
    print("Running DropletUtils analysis (EmptyDrops + barcodeRanks) on pre-filtered matrix...")
    run_emptydrops_r(args.kb_dir, args.sample_id, output_dir, 
                     params.get('emptydrops_lower', 100), 
                     args.ncores)
    
    dropletutils_results = load_dropletutils_results(output_dir, args.sample_id, adata)
    methods_results.update(dropletutils_results)
    
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