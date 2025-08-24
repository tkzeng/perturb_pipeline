#!/usr/bin/env python3
"""
Pseudobulk expression with TMM normalization using scanpy and edgeR
Processes the X matrix (not layers) for pseudobulk aggregation
Modified from statin_perturb pipeline to work with standard downstream analysis
"""

import argparse
import os
import pandas as pd
import numpy as np
import scanpy as sc
from scipy import sparse
import subprocess


def process_pseudobulk(adata, groupby_cols, covariate_cols):
    """Process X matrix for pseudobulk aggregation"""
    
    # Ensure grouping columns are strings
    adata_copy = adata.copy()
    for col in groupby_cols:
        adata_copy.obs[col] = adata_copy.obs[col].astype(str)
    
    # Pseudobulk aggregation using X matrix
    pdata = sc.get.aggregate(adata_copy, by=groupby_cols, func="sum", axis=0)
    
    # Ensure counts are integers
    if sparse.issparse(pdata.layers["sum"]):
        pdata.layers["sum"].data = np.round(pdata.layers["sum"].data).astype(int)
    else:
        pdata.layers["sum"] = np.round(pdata.layers["sum"]).astype(int)
    
    # Calculate cell counts and library sizes
    cell_counts = adata.obs.groupby(groupby_cols).size()
    
    # Create the same index structure that scanpy aggregate uses
    if len(groupby_cols) == 1:
        # Single column: use values directly
        cell_counts.index = cell_counts.index.astype(str)
    else:
        # Multiple columns: join with underscore
        cell_counts.index = cell_counts.index.map(lambda x: '_'.join(map(str, x)))
    
    if sparse.issparse(pdata.layers["sum"]):
        lib_sizes = np.array(pdata.layers["sum"].sum(axis=1)).flatten()
    else:
        lib_sizes = pdata.layers["sum"].sum(axis=1)
    
    # Add cell counts and library sizes to pdata.obs
    pdata.obs['n_cells'] = cell_counts.reindex(pdata.obs.index)
    pdata.obs['manual_total_lib_size'] = lib_sizes
    
    # Aggregate covariates if provided
    if covariate_cols:
        # Only aggregate numeric covariates that exist
        numeric_covariates = []
        for col in covariate_cols:
            if col in adata.obs.columns:
                if pd.api.types.is_numeric_dtype(adata.obs[col]):
                    numeric_covariates.append(col)
                else:
                    print(f"Skipping non-numeric covariate: {col}")
            else:
                print(f"Warning: Covariate column '{col}' not found in adata.obs")
        
        if numeric_covariates:
            covariate_stats = adata.obs.groupby(groupby_cols)[numeric_covariates].agg(['mean', 'median', 'std'])
            
            # Create the same index structure for covariate stats
            if len(groupby_cols) == 1:
                covariate_stats.index = covariate_stats.index.astype(str)
            else:
                covariate_stats.index = covariate_stats.index.map(lambda x: '_'.join(map(str, x)))
            
            # Flatten MultiIndex columns properly
            covariate_stats.columns = [f'{stat}_{col}' for col, stat in covariate_stats.columns]
            
            # Error checking: detect constant values (std = 0)
            std_cols = [col for col in covariate_stats.columns if col.startswith('std_')]
            for std_col in std_cols:
                zero_std_groups = covariate_stats[covariate_stats[std_col] == 0.0]
                if len(zero_std_groups) > 0:
                    covariate_name = std_col.replace('std_', '')
                    print(f"Warning: Found {len(zero_std_groups)} groups with constant {covariate_name} values (std=0)")
                    print(f"Groups: {zero_std_groups.index.tolist()[:5]}...")  # Show first 5
            
            # Merge covariate stats  
            for col in covariate_stats.columns:
                pdata.obs[col] = covariate_stats[col].reindex(pdata.obs.index)
    
    # Filter zero-count samples
    valid_samples = pdata.obs['manual_total_lib_size'] > 0
    if not valid_samples.all():
        print(f"Removing {(~valid_samples).sum()} samples with zero counts")
        pdata = pdata[valid_samples]
    
    return pdata


def run_tmm_workflow(adata, groupby_cols, covariate_cols, output_prefix):
    """Run TMM workflow: process X matrix -> edgeR -> save results"""
    
    print(f"\n=== Processing X matrix for pseudobulk ===")
    pdata = process_pseudobulk(adata, groupby_cols, covariate_cols)
    
    # Create counts DataFrame
    counts_df = pd.DataFrame(
        pdata.layers["sum"].T if sparse.issparse(pdata.layers["sum"]) else pdata.layers["sum"].T,
        index=pdata.var.index,
        columns=pdata.obs.index
    )
    counts_df.insert(0, 'gene', pdata.var.get('gene', pdata.var.index))
    
    # Save raw counts and run edgeR
    counts_df.to_csv(f"{output_prefix}_raw_counts_all_genes.tsv", sep="\t")
    input_file = f"{output_prefix}_input_counts.tsv"
    counts_df.reset_index().to_csv(input_file, sep="\t", index=False)
    
    # Run edgeR
    cmd = ["Rscript", "scripts/run_edger_tmm.R", input_file, output_prefix]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"edgeR failed: {result.stderr}")
    
    # Load edgeR results
    norm_factors_df = pd.read_csv(f"{output_prefix}_norm_factors.tsv", sep="\t")
    gene_filter_df = pd.read_csv(f"{output_prefix}_gene_filter.tsv", sep="\t")
    
    # Create metadata tables
    save_metadata_tables(pdata, norm_factors_df, gene_filter_df, groupby_cols, output_prefix)
    
    # Cleanup
    os.remove(input_file)
    
    keep_genes = gene_filter_df.set_index('gene_id')['kept_by_edger'].reindex(pdata.var.index).fillna(False).values
    print(f"Completed: {keep_genes.sum()}/{len(keep_genes)} genes kept by edgeR filtering")
    
    return pdata, norm_factors_df, gene_filter_df


def save_metadata_tables(pdata, norm_factors_df, gene_filter_df, groupby_cols, output_prefix):
    """Save group and gene metadata tables"""
    
    # Group metadata table
    metadata_df = pd.DataFrame({
        'cell_count': pdata.obs['n_cells'],
        'manual_total_lib_size': pdata.obs['manual_total_lib_size'],
        'edger_lib_size': norm_factors_df['lib_size'].values,
        'tmm_norm_factor': norm_factors_df['norm_factors'].values,
    }, index=pdata.obs.index)
    
    # Add grouping columns and covariate stats
    for col in groupby_cols:
        if col in pdata.obs.columns:
            metadata_df[col] = pdata.obs[col]
    
    for col in [c for c in pdata.obs.columns if any(c.startswith(f"{stat}_") for stat in ['mean', 'median', 'std'])]:
        metadata_df[col] = pdata.obs[col]
    
    metadata_df['effective_library_size'] = metadata_df['edger_lib_size'] * metadata_df['tmm_norm_factor']
    
    # Set index name and save with proper formatting
    metadata_df.index.name = 'group_id'
    metadata_df.to_csv(f"{output_prefix}_group_metadata.tsv", sep="\t")
    
    # Gene metadata table
    keep_genes = gene_filter_df.set_index('gene_id')['kept_by_edger'].reindex(pdata.var.index).fillna(False).values
    gene_metadata_df = pd.DataFrame({
        'kept_by_edger': keep_genes,
        'gene': pdata.var.get('gene', pdata.var.index),
    }, index=pdata.var.index)
    gene_metadata_df.index.name = 'gene_id'
    gene_metadata_df.to_csv(f"{output_prefix}_gene_metadata.tsv", sep="\t")


def pseudobulk_with_tmm(adata, groupby_cols, covariate_cols, output_prefix):
    """Main pseudobulk analysis with TMM normalization"""
    
    # Validate inputs
    missing_cols = [col for col in groupby_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Missing groupby columns in adata.obs: {missing_cols}")
    
    # Process X matrix with TMM normalization
    print(f"Processing X matrix with shape: {adata.shape}")
    print(f"Grouping by: {groupby_cols}")
    if covariate_cols:
        print(f"Aggregating covariates: {covariate_cols}")
    
    pdata, norm_factors_df, gene_filter_df = run_tmm_workflow(
        adata, groupby_cols, covariate_cols, output_prefix
    )
    
    print(f"\nPseudobulk processing complete!")
    print(f"Created {pdata.n_obs} pseudobulk samples from {adata.n_obs} cells")
    
    return pdata


def main():
    parser = argparse.ArgumentParser(description="Pseudobulk TMM normalization for X matrix")
    parser.add_argument("--input-file", required=True, help="Input h5ad file (preprocessed with metadata)")
    parser.add_argument("--raw-file", required=True, help="Raw h5ad file with X matrix (combined file)")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--groupby-cols", nargs='+', required=True, help="Grouping columns")
    parser.add_argument("--covariate-cols", nargs='*', default=[], help="Covariate columns (optional)")
    parser.add_argument("--output-prefix", required=True, help="Output prefix")
    
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load metadata from preprocessed file
    print(f"Loading metadata from: {args.input_file}")
    adata = sc.read_h5ad(args.input_file)
    print(f"Loaded metadata: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Load raw file in backed mode to get X matrix efficiently
    print(f"Loading raw counts from: {args.raw_file} (backed mode)")
    adata_raw = sc.read_h5ad(args.raw_file, backed='r')
    
    # Extract only the X matrix for the cells and genes we need
    # Using [:] forces loading the subset into memory as an array
    print("Extracting raw counts for preprocessed cells and genes...")
    adata.X = adata_raw[adata.obs.index, adata.var.index].X[:]
    
    print(f"Replaced X matrix with raw counts: {adata.X.shape}")
    
    # Check if X contains raw counts or normalized data
    if 'counts' in adata.layers:
        print("Warning: Found 'counts' layer. Consider using raw counts for pseudobulking.")
    
    # Process with TMM
    output_prefix = os.path.join(args.output_dir, args.output_prefix)
    result = pseudobulk_with_tmm(
        adata, args.groupby_cols, args.covariate_cols or [], output_prefix
    )
    
    print(f"\nAll output files saved with prefix: {output_prefix}")


if __name__ == "__main__":
    main()