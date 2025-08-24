#!/usr/bin/env python3
"""
Guide analysis module containing all guide-related functions.

This module provides functions for:
- Fitting Poisson-Gaussian mixture models to guide UMI distributions
- Calculating posterior probability thresholds
- Guide metrics calculation
- Guide threshold application
"""

import numpy as np
import pandas as pd
import scipy.sparse
import scipy.sparse as sp
from pathlib import Path
import sys
import os
from functools import partial
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
from scipy.stats import poisson, norm
from multiprocessing import Pool, cpu_count
import time

# Set global matplotlib parameters
plt.rcParams['savefig.dpi'] = 100
from itertools import product
import argparse
import scanpy as sc
import yaml

# Logging function
def log_print(message):
    """Print a message to stdout for logging."""
    print(message, flush=True)


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


def get_cutoff_suffix(cutoff_spec):
    """Convert cutoff specification to file suffix.
    
    Args:
        cutoff_spec: Either int or string like 'gmm_75'
        
    Returns:
        String suffix like 'cutoff1' or 'gmm75', or None if invalid
    """
    if isinstance(cutoff_spec, int):
        return f"cutoff{cutoff_spec}"
    elif isinstance(cutoff_spec, str) and cutoff_spec.startswith('gmm_'):
        level = cutoff_spec[4:]
        return f"gmm{level}"
    elif isinstance(cutoff_spec, str) and cutoff_spec.startswith('posterior_'):
        level = cutoff_spec[10:]
        return f"gmm{level}"
    return None


def map_guide_counts_to_gex_barcodes(guide_adata, gex_barcodes):
    """Map guide counts to GEX barcodes, filling zeros for missing cells.
    
    Args:
        guide_adata: Guide AnnData (may only have barcodes with reads)
        gex_barcodes: Complete list of GEX barcodes
        
    Returns:
        Sparse matrix with rows for ALL gex_barcodes (zeros for missing)
    """
    n_gex = len(gex_barcodes)
    n_guides = guide_adata.shape[1]
    
    barcode_series = pd.Series(range(guide_adata.n_obs), index=guide_adata.obs_names)
    cell_indices_series = barcode_series.reindex(gex_barcodes)
    
    found_mask = ~cell_indices_series.isna()
    n_found = found_mask.sum()
    
    if n_found == 0:
        return scipy.sparse.csr_matrix((n_gex, n_guides))
    
    found_indices = cell_indices_series.dropna().astype(int).values
    guide_counts_found = guide_adata.X[found_indices]
    
    full_matrix = scipy.sparse.lil_matrix((n_gex, n_guides))
    gex_indices = np.where(found_mask)[0]
    full_matrix[gex_indices] = guide_counts_found
    
    return full_matrix.tocsr()


def calculate_guide_metrics_for_cells(guide_adata, cell_barcodes, guide_cutoffs, obs_data, method_name, stratify_by):
    """Calculate guide metrics for cells, including those with no guides.
    
    Args:
        guide_adata: Guide count AnnData object (may be missing cells)
        cell_barcodes: List/array of cell barcodes to calculate metrics for
        guide_cutoffs: List of cutoff specifications (int or 'gmm_XX' strings)
        obs_data: DataFrame with per-cell metadata including GMM threshold columns
        method_name: Cell calling method name (required for GMM thresholds)
        stratify_by: Stratification level ('sample' or 'biological_sample')

    Returns:
        dict: Metrics including:
            - n_cells_total: Total cells provided
            - n_cells_in_guide_data: Cells found in guide data
            - guide_umis_per_cell: Mean guide UMIs (including 0s)
            - For each cutoff with appropriate suffix:
                - mean_guides_per_cell_{suffix}: Mean guides per cell
                - median_guides_per_cell_{suffix}: Median guides per cell
                - fraction_cells_with_guides_{suffix}: Fraction with ≥1 guide
                - umis_per_guide_per_cell_{suffix}: Mean UMIs per guide (cells with guides only)
    """
    metrics = {}
    n_cells_total = len(cell_barcodes)
    metrics['n_cells_total'] = n_cells_total

    if n_cells_total == 0:
        raise ValueError("No cell barcodes provided")

    # Check if guide_adata already contains exactly the requested cells
    # This avoids redundant remapping when the data is already properly aligned
    if set(guide_adata.obs_names) == set(cell_barcodes) and len(guide_adata.obs_names) == len(cell_barcodes):
        # Data already contains exactly the requested cells - no remapping needed
        mapped_counts = guide_adata.X
        n_cells_found = n_cells_total
    else:
        # Need to remap because guide_adata has different cells than requested
        mapped_counts = map_guide_counts_to_gex_barcodes(guide_adata, cell_barcodes)
        n_cells_found = (guide_adata.obs_names.isin(cell_barcodes)).sum()
    
    metrics['n_cells_in_guide_data'] = int(n_cells_found)

    log_print(f"Found {n_cells_found:,} of {n_cells_total:,} cells in guide data ({n_cells_found/n_cells_total*100:.1f}%)")

    all_umis_per_cell = np.array(mapped_counts.sum(axis=1)).flatten()
    guide_counts_sparse = mapped_counts

    metrics['guide_umis_per_cell'] = float(np.mean(all_umis_per_cell))

    for cutoff_spec in guide_cutoffs:
        suffix = get_cutoff_suffix(cutoff_spec)
        if suffix is None:
            continue
        
        if guide_counts_sparse.nnz > 0:
            binary_guides = apply_guide_threshold(guide_counts_sparse, cutoff_spec, obs_data, stratify_by, method_name)
            all_guides_per_cell = np.array(binary_guides.sum(axis=1)).flatten()
        else:
            all_guides_per_cell = np.zeros(n_cells_total)

        metrics[f'mean_guides_per_cell_{suffix}'] = float(np.mean(all_guides_per_cell))
        metrics[f'median_guides_per_cell_{suffix}'] = float(np.median(all_guides_per_cell))

        cells_with_guides = all_guides_per_cell > 0
        metrics[f'fraction_cells_with_guides_{suffix}'] = float(np.sum(cells_with_guides) / n_cells_total)

        if np.any(cells_with_guides):
            guide_umis_with_guides = all_umis_per_cell[cells_with_guides]
            guides_per_cell_with_guides = all_guides_per_cell[cells_with_guides]
            umis_per_guide = guide_umis_with_guides / guides_per_cell_with_guides
            metrics[f'umis_per_guide_per_cell_{suffix}'] = float(np.mean(umis_per_guide))
        else:
            metrics[f'umis_per_guide_per_cell_{suffix}'] = 0.0

    return metrics


def parse_guide_cutoffs(guide_cutoffs_config, gmm_thresholds):
    """
    Parse mixed guide cutoff specifications from config.
    
    Args:
        guide_cutoffs_config: List from config containing mixed threshold specs
        gmm_thresholds: Dict mapping posterior levels to UMI thresholds
        
    Returns:
        List of dicts with keys: type, value, label, threshold_umi
    """
    if guide_cutoffs_config is None:
        raise ValueError("guide_cutoffs_config cannot be None")
    
    if gmm_thresholds is None:
        gmm_thresholds = {}
    
    parsed_cutoffs = []
    
    for cutoff_spec in guide_cutoffs_config:
        if isinstance(cutoff_spec, int):
            parsed_cutoffs.append({
                'type': 'fixed',
                'value': cutoff_spec,
                'label': f'≥{cutoff_spec} UMIs',
                'threshold_umi': cutoff_spec
            })
        elif isinstance(cutoff_spec, str):
            if cutoff_spec.startswith('gmm_') or cutoff_spec.startswith('posterior_'):
                # Extract posterior probability level
                if cutoff_spec.startswith('gmm_'):
                    prob_str = cutoff_spec[4:]  # Remove 'gmm_'
                else:
                    prob_str = cutoff_spec[10:]  # Remove 'posterior_'
                
                try:
                    prob_level = int(prob_str)
                    if prob_level not in gmm_thresholds:
                        log_print(f"WARNING: GMM threshold for {prob_level}% not available, skipping {cutoff_spec}")
                        continue
                    
                    threshold_umi = gmm_thresholds[prob_level]
                    parsed_cutoffs.append({
                        'type': 'gmm',
                        'value': prob_level,
                        'label': f'{prob_level}% confidence (≥{threshold_umi} UMIs)',
                        'threshold_umi': threshold_umi
                    })
                except ValueError:
                    raise ValueError(f"Invalid posterior probability specification: {cutoff_spec}")
            else:
                raise ValueError(f"Unknown cutoff specification: {cutoff_spec}")
        else:
            raise ValueError(f"Cutoff must be int or str, got {type(cutoff_spec)}: {cutoff_spec}")
    
    if len(parsed_cutoffs) == 0:
        raise ValueError("No valid cutoffs after parsing")
    
    return parsed_cutoffs



def load_guide_reference(guide_ref_path):
    """Load guide reference and create mapping from guide sequence to list of target genes.
    
    The multi-gene separator is read from a 'separator' column in the file if present,
    otherwise defaults to '_AND_' for backwards compatibility.
    
    Args:
        guide_ref_path: Path to the guide reference TSV file
    """
    
    log_print(f"Loading guide reference from: {guide_ref_path}")
    
    # Read the TSV file
    guide_df = pd.read_csv(guide_ref_path, sep='\t')
    
    # Get separator from file if present, otherwise use default
    if 'separator' in guide_df.columns:
        # Get the first non-null separator value
        separator_values = guide_df['separator'].dropna().unique()
        if len(separator_values) > 0:
            multi_gene_separator = str(separator_values[0])
            log_print(f"Using multi-gene separator from file: '{multi_gene_separator}'")
        else:
            multi_gene_separator = '_AND_'
            log_print(f"No separator specified in file, using default: '{multi_gene_separator}'")
    else:
        multi_gene_separator = '_AND_'
        log_print(f"No separator column in file, using default: '{multi_gene_separator}'")
    
    # Create mapping from guide sequence to list of genes - vectorized approach
    guide_to_genes = {}
    
    # Filter out NaN or empty gene entries
    valid_mask = guide_df['gene'].notna() & (guide_df['gene'] != '')
    valid_df = guide_df[valid_mask].copy()
    
    if len(valid_df) > 0:
        # Vectorized processing
        guide_ids = valid_df['ID'].values
        gene_strs = valid_df['gene'].values
        
        for guide_id, gene_str in zip(guide_ids, gene_strs):
            # Split multi-gene targets if separator is configured
            if multi_gene_separator and multi_gene_separator in gene_str:
                genes = gene_str.split(multi_gene_separator)
            else:
                genes = [gene_str]
            
            guide_to_genes[guide_id] = genes
    
    log_print(f"Loaded {len(guide_to_genes)} guide IDs mapping to {len(set(sum(guide_to_genes.values(), [])))} unique genes")
    
    return guide_to_genes


def calculate_guide_metrics_for_barcodes(adata, barcode_set, guide_to_genes=None, cutoffs=None, stratify_by='biological_sample', method_name=None):
    """Calculate guide metrics for a given set of barcodes at multiple cutoffs.
    
    Now uses per-cell thresholds from adata.obs for GMM-based cutoffs.
    
    Args:
        method_name: Cell calling method name (required for GMM thresholds)
    """
    if cutoffs is None:
        raise ValueError("Guide cutoffs must be provided")
    
    # Convert barcode_set to list for the unified function
    if barcode_set:
        barcode_list = list(barcode_set)
    else:
        # Empty set - no cells
        barcode_list = []
    
    if len(barcode_list) == 0:
        # Return empty metrics when no barcodes are available
        metrics = {
            'guide_umis_per_cell': 0
        }
        # Add metrics for each cutoff matching the expected format
        if cutoffs:
            for cutoff_spec in cutoffs:
                suffix = get_cutoff_suffix(cutoff_spec)
                if suffix:
                    metrics[f'guides_per_cell_{suffix}'] = 0
                    metrics[f'fraction_cells_with_guides_{suffix}'] = 0.0
                    metrics[f'umis_per_guide_per_cell_{suffix}'] = 0
                    metrics[f'total_guides_detected_{suffix}'] = 0
                    metrics[f'total_guide_detections_{suffix}'] = 0
                    metrics[f'total_targeted_genes_{suffix}'] = 0
        return metrics
    
    # Extract guide data from adata
    # Get subset of cells we're analyzing
    adata_subset = adata[adata.obs.index.isin(barcode_list)]
    guide_counts = adata_subset.obsm['guide_counts']
    
    # Calculate total guide UMIs per cell
    guide_umis_per_cell = np.array(guide_counts.sum(axis=1)).flatten()
    metrics = {
        'guide_umis_per_cell': np.mean(guide_umis_per_cell) if len(guide_umis_per_cell) > 0 else 0
    }
    
    # Process each cutoff
    if cutoffs and len(cutoffs) > 0:
        for cutoff_spec in cutoffs:
            # Get file suffix for this cutoff
            suffix = get_cutoff_suffix(cutoff_spec)
            if not suffix:
                continue
            
            # Apply threshold using helper function
            binary_guides = apply_guide_threshold(guide_counts, cutoff_spec, adata_subset.obs, stratify_by, method_name)
            
            # Calculate metrics with the binary guide matrix
            guides_per_cell = np.array(binary_guides.sum(axis=1)).flatten()
            cells_with_guides = guides_per_cell > 0
            
            # Basic metrics
            metrics[f'guides_per_cell_{suffix}'] = np.mean(guides_per_cell) if len(guides_per_cell) > 0 else 0
            metrics[f'fraction_cells_with_guides_{suffix}'] = np.mean(cells_with_guides) if len(cells_with_guides) > 0 else 0.0
            
            # UMIs per guide per cell (only for cells with guides)
            if np.any(cells_with_guides):
                # Only consider UMIs from guides that pass the cutoff
                guide_counts_filtered = guide_counts[cells_with_guides].multiply(binary_guides[cells_with_guides])
                umis_per_cell_with_guides = np.array(guide_counts_filtered.sum(axis=1)).flatten()
                guides_per_cell_with_guides = guides_per_cell[cells_with_guides]
                umis_per_guide = umis_per_cell_with_guides / guides_per_cell_with_guides
                metrics[f'umis_per_guide_per_cell_{suffix}'] = np.mean(umis_per_guide)
            else:
                metrics[f'umis_per_guide_per_cell_{suffix}'] = 0
            
            # Total guides detected across all cells (binary_guides is sparse)
            # Get column-wise sum efficiently with sparse matrix
            guides_per_column = np.array(binary_guides.sum(axis=0)).flatten()
            guides_detected_by_any_cell = guides_per_column > 0
            metrics[f'total_guides_detected_{suffix}'] = np.sum(guides_detected_by_any_cell)
            
            # Total guide detections (sum of all guide-cell pairs)
            metrics[f'total_guide_detections_{suffix}'] = np.sum(binary_guides)
            
            # Total targeted genes (if guide_to_genes mapping provided)
            if guide_to_genes and 'guide_names' in adata.uns:
                guide_names = adata.uns['guide_names']
                detected_guide_names = [guide_names[i] for i in np.where(guides_detected_by_any_cell)[0]]
                targeted_genes = set()
                for guide_name in detected_guide_names:
                    if guide_name in guide_to_genes:
                        # guide_to_genes[guide_name] is always a list of genes
                        targeted_genes.update(guide_to_genes[guide_name])
                metrics[f'total_targeted_genes_{suffix}'] = len(targeted_genes)
            else:
                metrics[f'total_targeted_genes_{suffix}'] = 0
    
    return metrics


def get_threshold_for_spec(cutoff_spec, parsed_cutoffs):
    """Get the actual UMI threshold for a cutoff specification.
    
    Args:
        cutoff_spec: The specification from config (int or string like 'gmm_75')
        parsed_cutoffs: List of parsed cutoff dicts from parse_guide_cutoffs()
        
    Returns:
        int: The UMI threshold, or None if not found
    """
    for parsed in parsed_cutoffs:
        if isinstance(cutoff_spec, int) and parsed['type'] == 'fixed' and parsed['value'] == cutoff_spec:
            return parsed['threshold_umi']
        elif isinstance(cutoff_spec, str) and parsed['type'] == 'gmm':
            if cutoff_spec.startswith('gmm_') and parsed['value'] == int(cutoff_spec[4:]):
                return parsed['threshold_umi']
            elif cutoff_spec.startswith('posterior_') and parsed['value'] == int(cutoff_spec[10:]):
                return parsed['threshold_umi']
    return None


def get_cutoff_labels(cutoff_spec):
    """Get file suffix and plot label for a cutoff specification.
    
    Args:
        cutoff_spec: The specification from config (int or string like 'gmm_75')
        
    Returns:
        tuple: (file_suffix, plot_label)
    """
    if isinstance(cutoff_spec, int):
        return f"cutoff{cutoff_spec}", f"≥{cutoff_spec} UMIs"
    elif isinstance(cutoff_spec, str) and (cutoff_spec.startswith('gmm_') or cutoff_spec.startswith('posterior_')):
        prob_level = int(cutoff_spec[4:] if cutoff_spec.startswith('gmm_') else cutoff_spec[10:])
        return f"gmm{prob_level}", f"{prob_level}% confidence"
    raise ValueError(f"Invalid cutoff specification: {cutoff_spec}")


def apply_guide_threshold(guide_counts, cutoff_spec, obs_data, stratify_by='biological_sample', method_name=None):
    """Apply guide threshold to create binary guide matrix.
    
    This helper function handles both fixed and GMM-based thresholds to avoid
    redundant code throughout the script.
    
    Args:
        guide_counts: Sparse matrix of guide counts (cells x guides)
        cutoff_spec: Either an int (fixed cutoff) or string like 'gmm_75'
        obs_data: DataFrame-like with per-cell metadata (must have GMM columns if using GMM cutoff)
        stratify_by: Stratification level ('sample' or 'biological_sample')
        method_name: Cell calling method name (required for GMM thresholds)
        
    Returns:
        numpy array: Binary guide matrix (1 if guide >= threshold, 0 otherwise)
    """
    if isinstance(cutoff_spec, int):
        # Fixed cutoff - apply same threshold to all cells
        return (guide_counts >= cutoff_spec).astype(int)
        
    elif isinstance(cutoff_spec, str) and cutoff_spec.startswith('gmm_'):
        # GMM cutoff - use per-cell thresholds from obs column, specific to method
        if method_name is None:
            raise ValueError("method_name is required for GMM thresholds")
        
        # Determine column name based on stratification level and method
        if stratify_by == 'sample':
            col_name = f"{cutoff_spec}_sample_{method_name}"
        else:  # biological_sample
            col_name = f"{cutoff_spec}_biosample_{method_name}"
        
        # Get per-cell thresholds
        thresholds = obs_data[col_name].values
        
        # Vectorized approach using sparse matrix operations
        # Convert to COO format for efficient element-wise operations
        guide_coo = guide_counts.tocoo()
        
        # Get threshold for each non-zero element
        element_thresholds = thresholds[guide_coo.row]
        
        # Apply thresholds element-wise
        binary_data = (guide_coo.data >= element_thresholds).astype(int)
        
        # Create new sparse matrix with binary values
        binary_guides = sp.coo_matrix(
            (binary_data, (guide_coo.row, guide_coo.col)),
            shape=guide_counts.shape,
            dtype=int
        ).tocsr()
        
        return binary_guides
    else:
        raise ValueError(f"Invalid cutoff specification: {cutoff_spec}")


def create_binary_guide_assignment(guide_counts, thresholds):
    """
    Create a binary guide assignment matrix from guide counts and per-cell thresholds.
    
    Args:
        guide_counts: Sparse matrix of guide counts (cells x guides)
        thresholds: Array of threshold values per cell or single threshold value
        
    Returns:
        Sparse binary matrix indicating guide assignments (cells x guides)
    """
    # Handle single threshold value
    if np.isscalar(thresholds):
        return (guide_counts >= thresholds).astype(int)
    
    # Handle per-cell thresholds
    thresholds = np.asarray(thresholds)
    
    # Convert guide_counts to COO format for efficient element-wise operations
    guide_coo = guide_counts.tocoo()
    
    # Get threshold for each non-zero element
    element_thresholds = thresholds[guide_coo.row]
    
    # Apply thresholds element-wise
    binary_data = (guide_coo.data >= element_thresholds).astype(int)
    
    # Create binary matrix
    binary_guides = sp.coo_matrix(
        (binary_data, (guide_coo.row, guide_coo.col)),
        shape=guide_counts.shape,
        dtype=int
    ).tocsr()
    
    return binary_guides


def plot_top_guide_distributions(adata, output_dir, method_name, stratify_by):
    """
    Plot UMI distributions for top guides for a single cell calling method.
    
    Args:
        adata: AnnData object with cells from a single method, with guide_counts and guide_assignment in obsm
        output_dir: Directory to save plots
        method_name: Name of the cell calling method
        stratify_by: 'sample' or 'biological_sample'
    """
    plot_dir = Path(output_dir)
    
    log_print(f"  Generating top guide UMI distribution plots for {method_name}...")
    
    guide_counts = adata.obsm['guide_counts']
    guide_assignment = adata.obsm['guide_assignment']
    guide_names = adata.uns['guide_names']
    
    # Use guide_assignment matrix to get cells per guide
    cells_per_guide = np.array(guide_assignment.sum(axis=0)).flatten()
    
    if stratify_by == 'sample':
        # For sample level, create one figure with top 8 guides
        top_n = 8
        top_guide_indices = np.argsort(cells_per_guide)[-top_n:][::-1]
        
        # Create figure with 2x4 grid
        fig, axes = plt.subplots(2, 4, figsize=(16, 8))
        guide_axes = axes.flatten()
        
        # Plot each guide
        for plot_idx, guide_idx in enumerate(top_guide_indices):
            ax = guide_axes[plot_idx]
            guide_name = guide_names[guide_idx]
            
            # Get UMI counts for this guide
            guide_umis = np.array(guide_counts[:, guide_idx].todense()).flatten()
            guide_umis_nonzero = guide_umis[guide_umis > 0]
            
            if len(guide_umis_nonzero) > 0:
                # Create histogram with KDE
                sns.histplot(data=guide_umis_nonzero, bins=30, 
                            alpha=0.6, kde=True, ax=ax, 
                            log_scale=(True, False), color='steelblue')
            
            ax.set_xlabel('UMIs (log)', fontsize=9)
            ax.set_ylabel('Cells', fontsize=9)
            n_cells_with_guide = np.sum(guide_umis > 0)
            ax.set_title(f'{guide_name}\n{n_cells_with_guide:,} cells', fontsize=8)
            ax.grid(True, alpha=0.3, linestyle='--')
            ax.tick_params(axis='both', labelsize=8)
        
        # Overall title
        fig.suptitle(f'Top {top_n} Guide UMI Distributions ({method_name})', 
                    fontsize=14, y=1.02)
        
        # Adjust layout and save
        plt.tight_layout()
        metric_name = "top_guides_umi_distribution"
        metric_dir = plot_dir / metric_name / method_name / stratify_by / "linear"
        metric_dir.mkdir(parents=True, exist_ok=True)
        filepath = metric_dir / "plot.png"
        plt.savefig(filepath, bbox_inches='tight')
        plt.close()
        log_print(f"    Saved top {top_n} guides panel")
        
    elif stratify_by == 'biological_sample':
        # For biological sample level, create separate figures per group
        bio_samples = adata.obs['biological_sample'].dropna().unique()
        
        for bio_sample in bio_samples:
            # Get data for this biological sample
            bio_mask = adata.obs['biological_sample'] == bio_sample
            group_data = adata[bio_mask]
            group_guide_counts = group_data.obsm['guide_counts']
            group_guide_assignment = group_data.obsm['guide_assignment']
            
            # Get cells per guide for this group
            group_cells_per_guide = np.array(group_guide_assignment.sum(axis=0)).flatten()
            
            # Get top 4 guides for this group
            top_n = 4
            top_guide_indices = np.argsort(group_cells_per_guide)[-top_n:][::-1]
            
            # Create 2x2 figure for this biological sample
            fig, axes = plt.subplots(2, 2, figsize=(10, 8))
            guide_axes = axes.flatten()
            
            # Plot each guide
            for plot_idx, guide_idx in enumerate(top_guide_indices):
                ax = guide_axes[plot_idx]
                guide_name = guide_names[guide_idx]
                
                # Get UMI counts for this guide in this group
                guide_umis = np.array(group_guide_counts[:, guide_idx].todense()).flatten()
                guide_umis_nonzero = guide_umis[guide_umis > 0]
                
                if len(guide_umis_nonzero) > 0:
                    sns.histplot(data=guide_umis_nonzero, bins=30,
                                alpha=0.6, kde=True, ax=ax,
                                log_scale=(True, False), color='steelblue')
                
                ax.set_xlabel('UMIs (log)', fontsize=9)
                ax.set_ylabel('Cells', fontsize=9)
                n_cells_with_guide = np.sum(guide_umis > 0)
                ax.set_title(f'{guide_name}\n{n_cells_with_guide:,} cells', fontsize=8)
                ax.grid(True, alpha=0.3, linestyle='--')
                ax.tick_params(axis='both', labelsize=8)
            
            # Overall title
            fig.suptitle(f'Top {top_n} Guide UMI Distributions\n{bio_sample} ({method_name})',
                        fontsize=12, y=1.02)
            
            # Adjust layout and save
            plt.tight_layout()
            metric_name = "top_guides_umi_distribution"
            metric_dir = plot_dir / metric_name / method_name / stratify_by / bio_sample / "linear"
            metric_dir.mkdir(parents=True, exist_ok=True)
            filepath = metric_dir / "plot.png"
            plt.savefig(filepath, bbox_inches='tight')
            plt.close()
            
            log_print(f"    Saved {bio_sample} top {top_n} guides panel")



def poisson_gaussian_mixture_ll_batch(data, param_batch, min_umi_threshold):
    """
    Calculate log-likelihood for multiple parameter sets at once.
    
    Args:
        data: (n_data,) array of UMI counts (filtered ≥min_umi_threshold)
        param_batch: (n_params, 4) array with columns [lambda, mu, sigma, w]
        min_umi_threshold: Minimum UMI count threshold (1 or 2)
        
    Returns:
        (n_params,) array of log-likelihoods
    """
    # Pre-compute data terms
    data_int = data.astype(int)
    log_data = np.log10(data)
    
    # Extract parameters with shape (n_params, 1) for broadcasting
    lambdas = param_batch[:, 0:1]
    mus = param_batch[:, 1:2]
    sigmas = param_batch[:, 2:3]
    ws = param_batch[:, 3:4]
    
    # Truncated Poisson: P(X=k|X≥threshold) = P(X=k) / P(X≥threshold)
    p_poisson_regular = poisson.pmf(data_int[None, :], lambdas)
    p_0 = np.exp(-lambdas)  # P(X=0)
    
    if min_umi_threshold == 1:
        # P(X=k|X≥1) = P(X=k) / (1 - P(X=0))
        p_truncated = p_poisson_regular / (1 - p_0)
    elif min_umi_threshold == 2:
        # P(X=k|X≥2) = P(X=k) / (1 - P(X=0) - P(X=1))
        p_1 = lambdas * np.exp(-lambdas)  # P(X=1)
        p_truncated = p_poisson_regular / (1 - p_0 - p_1)
    else:
        raise ValueError(f"min_umi_threshold must be 1 or 2, got {min_umi_threshold}")
    
    # Gaussian component with change of variables
    p_gaussian = norm.pdf(log_data[None, :], mus, sigmas) / (data[None, :] * np.log(10))
    
    # Mixture probabilities
    p_mixture = ws * p_truncated + (1-ws) * p_gaussian
    
    # Log-likelihood for each parameter set
    lls = np.sum(np.log(np.maximum(p_mixture, 1e-300)), axis=1)
    
    return lls


def poisson_gaussian_mixture_ll_batch_cached(unique_vals, counts, param_batch, min_umi_threshold):
    """
    Calculate log-likelihood using cached unique values for speed.
    
    Args:
        unique_vals: (n_unique,) array of unique UMI counts
        counts: (n_unique,) array of counts for each unique value
        param_batch: (n_params, 4) array with columns [lambda, mu, sigma, w]
        min_umi_threshold: Minimum UMI count threshold (1 or 2)
        
    Returns:
        (n_params,) array of log-likelihoods
    """
    # Pre-compute data terms for unique values only
    data_int = unique_vals.astype(int)
    log_data = np.log10(unique_vals)
    
    # Extract parameters with shape (n_params, 1) for broadcasting
    lambdas = param_batch[:, 0:1]
    mus = param_batch[:, 1:2]
    sigmas = param_batch[:, 2:3]
    ws = param_batch[:, 3:4]
    
    # Truncated Poisson for unique values
    p_poisson_regular = poisson.pmf(data_int[None, :], lambdas)
    p_0 = np.exp(-lambdas)
    
    if min_umi_threshold == 1:
        p_truncated = p_poisson_regular / (1 - p_0)
    elif min_umi_threshold == 2:
        p_1 = lambdas * np.exp(-lambdas)
        p_truncated = p_poisson_regular / (1 - p_0 - p_1)
    else:
        raise ValueError(f"min_umi_threshold must be 1 or 2, got {min_umi_threshold}")
    
    # Gaussian component for unique values
    p_gaussian = norm.pdf(log_data[None, :], mus, sigmas) / (unique_vals[None, :] * np.log(10))
    
    # Mixture probabilities
    p_mixture = ws * p_truncated + (1 - ws) * p_gaussian
    
    # Log-likelihood weighted by counts
    # Use counts to weight the contribution of each unique value
    log_p_mixture = np.log(p_mixture + 1e-300)
    lls = np.sum(log_p_mixture * counts[None, :], axis=1)
    
    return lls


def fit_poisson_gaussian_mixture(data, min_umi_threshold, subsample_size, n_cores, verbose):
    """
    Fit Poisson-Gaussian mixture model to guide UMI data.
    
    Args:
        data: Array of UMI counts (already filtered ≥min_umi_threshold)
        min_umi_threshold: Minimum UMI count threshold (1 or 2)
        subsample_size: Number of points to subsample for fitting
        n_cores: Number of CPU cores to use for parallel processing
        verbose: Whether to print progress messages
        
    Returns:
        dict with keys: lambda, mu, sigma, weight, log_likelihood
    """
    # Subsample data for speed
    if len(data) > subsample_size:
        data_subset = np.random.choice(data, subsample_size, replace=False)
        if verbose:
            print(f"  Subsampled from {len(data)} to {len(data_subset)} points")
    else:
        data_subset = data
    
    # Cache unique values to speed up likelihood calculation
    unique_vals, counts = np.unique(data_subset, return_counts=True)
    if verbose:
        print(f"  Found {len(unique_vals)} unique values (from {len(data_subset)} total points)")
    
    # Define parameter ranges - reduced grid for faster computation
    # Lambda range: log-spaced for better coverage across scale
    lambda_range = np.logspace(np.log10(0.1), np.log10(15.0), 30)  # 30 points from 0.1 to 15.0 in log space
    
    # Mu range: log10 space from ~2 UMIs to ~100 UMIs  
    mu_range = np.linspace(0.3, 2.0, 30)  # Increased to 30 points
    
    # Sigma range: typical spreads for log-normal distributions
    sigma_range = np.linspace(0.1, 1.0, 30)  # Increased to 30 points
    
    # Weight range: expect more Gaussian since we removed the spike
    weight_range = np.linspace(0.1, 0.8, 30)  # Increased to 30 points
    
    if verbose:
        print(f"  Parameter ranges:")
        print(f"    λ: [{lambda_range[0]:.2f}, {lambda_range[-1]:.2f}]")
        print(f"    μ: [{mu_range[0]:.2f}, {mu_range[-1]:.2f}] (log10 scale)")
        print(f"    σ: [{sigma_range[0]:.2f}, {sigma_range[-1]:.2f}]")
        print(f"    w: [{weight_range[0]:.2f}, {weight_range[-1]:.2f}]")
    
    # Create parameter grid
    param_grid = np.array(list(product(lambda_range, mu_range, sigma_range, weight_range)))
    total_combinations = len(param_grid)
    
    if verbose:
        print(f"  Testing {total_combinations:,} parameter combinations using {n_cores} cores...")
    
    # Process in chunks using multiprocessing
    chunk_size = 10000
    chunks = [param_grid[i:min(i + chunk_size, total_combinations)] 
              for i in range(0, total_combinations, chunk_size)]
    
    # Create partial function with unique values, counts, and threshold
    evaluate_chunk = partial(poisson_gaussian_mixture_ll_batch_cached, unique_vals, counts, min_umi_threshold=min_umi_threshold)
    
    best_ll = -np.inf
    best_params = None
    
    start_time = time.time()
    
    # Process chunks in parallel
    with Pool(n_cores) as pool:
        for i, (chunk, lls) in enumerate(zip(chunks, pool.map(evaluate_chunk, chunks))):
            # Find best in this chunk
            chunk_best_idx = np.argmax(lls)
            if lls[chunk_best_idx] > best_ll:
                best_ll = lls[chunk_best_idx]
                best_params = chunk[chunk_best_idx]
            
            # Progress update
            if verbose and i % 5 == 0:
                progress = min((i + 1) * chunk_size, total_combinations)
                elapsed = time.time() - start_time
                rate = progress / elapsed if elapsed > 0 else 1
                remaining = (total_combinations - progress) / rate if rate > 0 else 0
                print(f"    Progress: {progress:,}/{total_combinations:,} ({100*progress/total_combinations:.1f}%) - "
                      f"Est. {remaining:.0f}s remaining")
    
    lambda_mix, mu_mix, sigma_mix, weight_mix = best_params
    
    return {
        'lambda': lambda_mix,
        'mu': mu_mix,
        'sigma': sigma_mix,
        'weight': weight_mix,
        'log_likelihood': best_ll
    }


def calculate_posterior_probabilities(data, params, min_umi_threshold):
    """
    Calculate posterior probabilities for component assignment.
    
    Args:
        data: Array of UMI counts
        params: Dictionary with keys lambda, mu, sigma, weight
        min_umi_threshold: Minimum UMI count threshold (1 or 2)
        
    Returns:
        dict with keys: posterior_poisson, posterior_gaussian
    """
    # Extract parameters
    lambda_mix = params['lambda']
    mu_mix = params['mu']
    sigma_mix = params['sigma']
    weight_mix = params['weight']
    
    # Calculate likelihoods
    data_int = data.astype(int)
    data_log = np.log10(data)
    
    # Truncated Poisson likelihood
    p_poisson_regular = poisson.pmf(data_int, lambda_mix)
    p_0 = np.exp(-lambda_mix)
    
    if min_umi_threshold == 1:
        poisson_likes = p_poisson_regular / (1 - p_0)
    elif min_umi_threshold == 2:
        p_1 = lambda_mix * np.exp(-lambda_mix)
        poisson_likes = p_poisson_regular / (1 - p_0 - p_1)
    else:
        raise ValueError(f"min_umi_threshold must be 1 or 2, got {min_umi_threshold}")
    
    # Gaussian likelihood with change of variables
    gaussian_likes = norm.pdf(data_log, mu_mix, sigma_mix) / (data * np.log(10))
    
    # Posterior probabilities
    posterior_poisson = (weight_mix * poisson_likes) / \
                       (weight_mix * poisson_likes + (1-weight_mix) * gaussian_likes + 1e-300)
    posterior_gaussian = 1 - posterior_poisson
    
    return {
        'posterior_poisson': posterior_poisson,
        'posterior_gaussian': posterior_gaussian
    }


def perform_guide_mixture_analysis(guide_counts_data, min_umi_threshold, subsample_size, n_threads, posterior_levels):
    """
    Perform Poisson-Gaussian mixture analysis on guide UMI data.
    
    Args:
        guide_counts_data: Sparse matrix of guide counts for cells
        min_umi_threshold: Minimum UMI threshold for filtering
        subsample_size: Number of points to subsample for fitting
        n_threads: Number of threads for parallel processing
        posterior_levels: List of posterior probability levels (e.g., [50, 75, 90])
        
    Returns:
        dict: Analysis results including posterior thresholds and mixture parameters
    """
    import time
    
    # Pool guide data across cells for mixture analysis - optimized for sparse matrices
    start_pool = time.time()
    n_cells = guide_counts_data.shape[0]
    n_guides = guide_counts_data.shape[1]
    target_points = 50000
    
    # Efficient sparse matrix operations
    # Convert to COO format for efficient access to non-zero values
    guide_coo = guide_counts_data.tocoo()
    
    # Filter by minimum threshold - this is vectorized and fast
    mask = guide_coo.data >= min_umi_threshold
    filtered_data = guide_coo.data[mask]
    
    # If we have more data than needed, randomly sample
    if len(filtered_data) > target_points:
        indices = np.random.choice(len(filtered_data), target_points, replace=False)
        pooled_data = filtered_data[indices]
    else:
        pooled_data = filtered_data
    
    print(f"      Pooling guide data: {time.time() - start_pool:.2f}s ({len(pooled_data)} points from {n_cells}x{n_guides} sparse matrix)")
    
    if len(pooled_data) < 100:
        print(f"  Warning: Only {len(pooled_data)} data points for mixture analysis")
        return None
    
    # Fit mixture model
    start_fit = time.time()
    mixture_params = fit_poisson_gaussian_mixture(
        data=pooled_data, 
        min_umi_threshold=min_umi_threshold,
        subsample_size=subsample_size,
        n_cores=n_threads,
        verbose=False
    )
    print(f"      Fitting mixture model: {time.time() - start_fit:.2f}s")
    
    # Calculate posterior probabilities
    start_post = time.time()
    posteriors = calculate_posterior_probabilities(
        data=pooled_data,
        params=mixture_params,
        min_umi_threshold=min_umi_threshold
    )
    print(f"      Calculating posteriors: {time.time() - start_post:.2f}s")
    
    # Find UMI thresholds for requested posterior levels using fine grid
    start_thresh = time.time()
    posterior_thresholds = {}
    
    # Create fine grid for theoretical threshold calculation (0.01 resolution)
    test_range = np.arange(min_umi_threshold, 30, 0.01)
    test_log = np.log10(test_range)
    
    # Calculate posterior probabilities for test range
    p_0 = np.exp(-mixture_params['lambda'])
    if min_umi_threshold == 1:
        test_poisson_likes = poisson.pmf(np.floor(test_range), mixture_params['lambda']) / (1 - p_0)
    else:  # min_umi_threshold == 2
        p_1 = mixture_params['lambda'] * np.exp(-mixture_params['lambda'])
        test_poisson_likes = poisson.pmf(np.floor(test_range), mixture_params['lambda']) / (1 - p_0 - p_1)
    
    test_gaussian_likes = norm.pdf(test_log, mixture_params['mu'], mixture_params['sigma']) / (test_range * np.log(10))
    test_posterior_gaussian = ((1-mixture_params['weight']) * test_gaussian_likes) / \
                             (mixture_params['weight'] * test_poisson_likes + (1-mixture_params['weight']) * test_gaussian_likes + 1e-300)
    
    for level in posterior_levels:
        target_prob = level / 100.0
        
        # Find where posterior crosses the threshold
        idx = np.where(test_posterior_gaussian >= target_prob)[0]
        if len(idx) > 0:
            threshold_umi = test_range[idx[0]]
        else:
            threshold_umi = min_umi_threshold
        
        posterior_thresholds[level] = threshold_umi
    
    print(f"      Finding thresholds: {time.time() - start_thresh:.2f}s")
    
    # Calculate additional statistics
    start_stats = time.time()
    n_data = len(pooled_data)
    
    # Cells assigned to each component
    gaussian_mask = posteriors['posterior_gaussian'] > 0.5
    poisson_mask = ~gaussian_mask
    
    # Statistics for each component
    poisson_mean = np.mean(pooled_data[poisson_mask]) if np.any(poisson_mask) else 0
    poisson_median = np.median(pooled_data[poisson_mask]) if np.any(poisson_mask) else 0
    gaussian_mean = np.mean(pooled_data[gaussian_mask]) if np.any(gaussian_mask) else 0
    gaussian_median = np.median(pooled_data[gaussian_mask]) if np.any(gaussian_mask) else 0
    
    # Convert log parameters back to linear scale for interpretability
    gaussian_mean_umis = 10 ** mixture_params['mu']
    
    print(f"      Computing statistics: {time.time() - start_stats:.2f}s")
    
    return {
        'posterior_thresholds': posterior_thresholds,
        'mixture_params': mixture_params,
        'pooled_data': pooled_data,
        'lambda': mixture_params['lambda'],
        'mu': mixture_params['mu'],
        'sigma': mixture_params['sigma'],
        'weight': mixture_params['weight'],
        'log_likelihood': mixture_params['log_likelihood'],
        'n_data_points': n_data,
        'poisson_mean_umis': poisson_mean,
        'poisson_median_umis': poisson_median,
        'gaussian_mean_umis': gaussian_mean,
        'gaussian_median_umis': gaussian_median,
        'gaussian_mean_log10': mixture_params['mu'],
        'fraction_poisson': np.mean(poisson_mask),
        'fraction_gaussian': np.mean(gaussian_mask),
        'fraction_gaussian_50pct': np.mean(posteriors['posterior_gaussian'] > 0.5),
        'mean_log_likelihood': mixture_params['log_likelihood'] / n_data if n_data > 0 else 0
    }


def calculate_group_gmm_thresholds(adata_group, min_umi_threshold, posterior_levels, n_threads):
    """Calculate GMM thresholds for a specific group of cells.
    
    Args:
        adata_group: AnnData object for this group
        min_umi_threshold: Minimum UMI threshold for filtering
        posterior_levels: List of posterior probability levels (e.g., [50, 75, 90])
        n_threads: Number of threads for parallel processing
        
    Returns:
        dict: Mapping of posterior level to UMI threshold, or empty dict if failed
    """
    if 'guide_counts' not in adata_group.obsm:
        return {}
    
    analysis_result = perform_guide_mixture_analysis(
        guide_counts_data=adata_group.obsm['guide_counts'],
        min_umi_threshold=min_umi_threshold,
        subsample_size=5000,
        n_threads=n_threads,
        posterior_levels=posterior_levels
    )
    
    if analysis_result is not None:
        return analysis_result['posterior_thresholds']
    return {}


def generate_mixture_samples(n_samples, params, min_umi_threshold, random_seed):
    """
    Generate samples from the fitted mixture model.
    
    Args:
        n_samples: Number of samples to generate
        params: Dictionary with keys lambda, mu, sigma, weight
        min_umi_threshold: Minimum UMI count threshold (1 or 2)
        random_seed: Random seed for reproducibility
        
    Returns:
        dict with keys: poisson_samples, gaussian_samples, combined_samples
    """
    np.random.seed(random_seed)
    
    # Extract parameters
    lambda_mix = params['lambda']
    mu_mix = params['mu']
    sigma_mix = params['sigma']
    weight_mix = params['weight']
    
    # Sample from Poisson component (truncated)
    n_poisson_samples = int(n_samples * weight_mix)
    poisson_samples = []
    while len(poisson_samples) < n_poisson_samples:
        candidate = np.random.poisson(lambda_mix)
        if candidate >= min_umi_threshold:
            poisson_samples.append(candidate)
    poisson_samples = np.array(poisson_samples[:n_poisson_samples])
    
    # Sample from Gaussian component (in log space, then transform)
    n_gaussian_samples = n_samples - len(poisson_samples)
    gaussian_samples_log = np.random.normal(mu_mix, sigma_mix, n_gaussian_samples)
    gaussian_samples = 10**gaussian_samples_log
    # Filter to match threshold
    gaussian_samples = gaussian_samples[gaussian_samples >= min_umi_threshold]
    
    # Combine samples
    combined_samples = np.concatenate([poisson_samples, gaussian_samples])
    
    return {
        'poisson_samples': poisson_samples,
        'gaussian_samples': gaussian_samples,
        'combined_samples': combined_samples
    }


def plot_mixture_analysis(data, params, output_path, title_prefix, min_umi_threshold, posterior_levels, axes=None):
    """
    Generate mixture model analysis plots.
    
    Args:
        data: Array of UMI counts (filtered ≥min_umi_threshold)
        params: Dictionary with fitted parameters
        output_path: Path to save the plot (ignored if axes provided)
        title_prefix: Prefix for plot titles
        min_umi_threshold: Minimum UMI count threshold (1 or 2)
        posterior_levels: List of posterior probability levels to plot (e.g., [50, 75, 90])
        axes: Optional tuple of (ax1, ax2) to plot on existing axes
    """
    # Generate samples from the model
    samples_dict = generate_mixture_samples(
        n_samples=len(data),
        params=params,
        min_umi_threshold=min_umi_threshold,
        random_seed=42
    )
    model_samples = samples_dict['combined_samples']
    
    # Create figure or use provided axes
    if axes is None:
        fig, (ax_low_counts, ax_mixture) = plt.subplots(1, 2, figsize=(14, 6))
        standalone_figure = True
    else:
        ax_low_counts, ax_mixture = axes
        fig = None
        standalone_figure = False
    
    # Panel 1: Low counts histogram comparison (2-15 UMIs)
    low_counts_real = data[data <= 15]
    low_counts_model = model_samples[model_samples <= 15]
    
    # Create combined dataset for seaborn with proper binning
    hist_data = []
    for val in low_counts_real:
        hist_data.append({'UMIs': val, 'Type': 'Real data'})
    for val in low_counts_model:
        hist_data.append({'UMIs': val, 'Type': 'Model samples'})
    
    hist_df = pd.DataFrame(hist_data)
    
    # Plot histograms with seaborn for consistent binning
    sns.histplot(data=hist_df, x='UMIs', hue='Type', bins=14, alpha=0.7, 
                ax=ax_low_counts, stat='density', multiple='dodge')
    
    ax_low_counts.set_xlabel('UMIs', fontsize=10)
    ax_low_counts.set_ylabel('Density', fontsize=10)
    ax_low_counts.set_title(f'Low Counts: Real vs Model\n(Real: n={len(low_counts_real)}, Model: n={len(low_counts_model)})', fontsize=11)
    ax_low_counts.grid(True, alpha=0.3, linestyle='--')
    ax_low_counts.set_xlim(min_umi_threshold - 0.5, 15.5)
    
    # Panel 2: KDE comparison with posterior probability thresholds
    # Calculate thresholds using grid (0.1 resolution is sufficient for integer UMIs)
    test_range = np.arange(min_umi_threshold, 30, 0.1)
    test_log = np.log10(test_range)
    
    # Calculate posterior probabilities for test range
    p_0_test = np.exp(-params['lambda'])
    if min_umi_threshold == 1:
        test_poisson_likes = poisson.pmf(np.floor(test_range), params['lambda']) / (1 - p_0_test)
    else:  # min_umi_threshold == 2
        p_1_test = params['lambda'] * np.exp(-params['lambda'])
        test_poisson_likes = poisson.pmf(np.floor(test_range), params['lambda']) / (1 - p_0_test - p_1_test)
    
    test_gaussian_likes = norm.pdf(test_log, params['mu'], params['sigma']) / (test_range * np.log(10))
    test_posterior_gaussian = ((1-params['weight']) * test_gaussian_likes) / \
                             (params['weight'] * test_poisson_likes + (1-params['weight']) * test_gaussian_likes + 1e-300)
    
    # Find threshold UMIs for requested posterior probabilities
    thresholds = {}
    for level in posterior_levels:
        prob = level / 100.0
        idx = np.where(test_posterior_gaussian >= prob)[0]
        if len(idx) > 0:
            thresholds[prob] = test_range[idx[0]]
    
    # Plot KDEs
    sns.kdeplot(data=data, ax=ax_mixture,
               log_scale=(True, False), color='black', linewidth=3,
               label='Real data')
    if len(model_samples) > 0:
        sns.kdeplot(data=model_samples, ax=ax_mixture,
                   log_scale=(True, False), color='red', linewidth=2,
                   label='Model samples')
    
    # Add threshold lines with colors from matplotlib colormap
    colors = cm.get_cmap('tab10')
    
    for idx, (prob, umi_thresh) in enumerate(thresholds.items()):
        color = colors(idx % 10)  # Cycle through tab10 colors
        ax_mixture.axvline(umi_thresh, color=color, linestyle='--', 
                          linewidth=2, alpha=0.7, label=f'{int(prob*100)}%: {umi_thresh:.1f} UMIs')
    
    ax_mixture.set_xlabel('UMIs (log scale)', fontsize=10)
    ax_mixture.set_ylabel('Density', fontsize=10)
    ax_mixture.set_title('Data with Posterior Probability Thresholds', fontsize=11)
    ax_mixture.legend(fontsize=8)
    ax_mixture.grid(True, alpha=0.3, linestyle='--')
    
    # Only handle figure-level operations if standalone
    if standalone_figure:
        # Overall title
        fig.suptitle(f'{title_prefix} - Poisson-Gaussian Mixture Model (≥{min_umi_threshold} UMIs)', fontsize=14, y=1.02)
        
        # Adjust layout and save
        plt.tight_layout()
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()


def create_cell_threshold_table(adata, guide_cutoffs, posterior_levels, min_umi_threshold, n_threads, cell_barcodes_dict, plot_dir=None, guide_assignment_method=None, config=None):
    """Create a per-cell table with GMM thresholds at different stratification levels for all cell calling methods.
    
    This table contains the GMM threshold values for each cell based on the group it belongs to and the cell calling method.
    For example, all cells in a biological sample get the same biosample-level threshold for their method.
    
    Args:
        adata: AnnData object with guide counts in obsm['guide_counts']
        guide_cutoffs: List of guide cutoffs from config
        posterior_levels: List of posterior probability levels to calculate (e.g., [50, 75, 90])
        min_umi_threshold: Minimum UMI threshold for mixture model fitting
        n_threads: Number of threads for parallel processing
        cell_barcodes_dict: Dictionary mapping method names to sets of cell barcodes
        plot_dir: Optional directory to save mixture model plots (already includes per-sample path)
        guide_assignment_method: Method for determining guide assignments (for plotting)
        config: Configuration dictionary
        
    Returns:
        pd.DataFrame: Per-cell table with GMM thresholds for different stratification levels and methods
    """
    # Initialize threshold table
    threshold_df = pd.DataFrame(index=adata.obs.index)
    
    # Calculate sample-level GMM thresholds for each cell calling method
    print("Calculating sample-level GMM thresholds for each cell calling method...")
    
    for method_name, method_barcodes in cell_barcodes_dict.items():
        print(f"  Processing {method_name}...")
        
        # Filter to this method's cells
        method_mask = adata.obs.index.isin(method_barcodes)
        adata_method = adata[method_mask]
        
        if len(adata_method) == 0:
            print(f"    No cells found for {method_name}, skipping")
            continue
        
        # Perform full mixture analysis for this method
        sample_analysis = perform_guide_mixture_analysis(
            guide_counts_data=adata_method.obsm['guide_counts'],
            min_umi_threshold=min_umi_threshold,
            subsample_size=5000,
            n_threads=n_threads,
            posterior_levels=posterior_levels
        )
        
        if sample_analysis is not None:
            sample_gmm = sample_analysis['posterior_thresholds']
            
            # Generate plot if requested
            if plot_dir is not None:
                
                plot_dir = Path(plot_dir)
                # plot_dir already includes per_cell/{source}_{processing}/{sample_id}, so just append the rest
                metric_dir = plot_dir / "guide_mixture_analysis" / method_name / "sample" / "linear"
                metric_dir.mkdir(parents=True, exist_ok=True)
                
                plot_path = metric_dir / "plot.png"
                plot_mixture_analysis(
                    data=sample_analysis['pooled_data'],
                    params=sample_analysis['mixture_params'],
                    output_path=plot_path,
                    title_prefix=f"Sample Level ({method_name})",
                    min_umi_threshold=min_umi_threshold,
                    posterior_levels=posterior_levels
                )
        else:
            sample_gmm = {}
        
        # Add sample-level thresholds for this method's cells
        for level in posterior_levels:
            col_name = f'gmm_{level}_sample_{method_name}'
            # Initialize column with NaN for all cells (means "not applicable")
            if col_name not in threshold_df.columns:
                threshold_df[col_name] = np.nan
            if level in sample_gmm:
                threshold_df.loc[method_mask, col_name] = sample_gmm[level]
    
    # Calculate biological sample-level GMM thresholds if that column exists
    if 'biological_sample' in adata.obs.columns:
        print("\nCalculating biological sample-level GMM thresholds for each cell calling method...")
        bio_samples = adata.obs['biological_sample'].dropna().unique()
        
        for method_name, method_barcodes in cell_barcodes_dict.items():
            print(f"  Processing {method_name}...")
            
            # Filter to this method's cells
            method_mask = adata.obs.index.isin(method_barcodes)
            
            # Create multi-panel plot if requested
            if plot_dir is not None and len(bio_samples) > 1:
                
                n_groups = len(bio_samples)
                fig_width = min(20, 6 * min(n_groups, 3))
                fig_height = 6 * ((n_groups + 2) // 3)
                fig = plt.figure(figsize=(fig_width, fig_height))
            
            for i, bio_sample in enumerate(bio_samples):
                print(f"    Calculating for {bio_sample}...")
                bio_mask = (adata.obs['biological_sample'] == bio_sample) & method_mask
                
                if bio_mask.sum() == 0:
                    continue
                
                # Perform full analysis for this biological sample and method
                bio_analysis = perform_guide_mixture_analysis(
                    guide_counts_data=adata[bio_mask].obsm['guide_counts'],
                    min_umi_threshold=min_umi_threshold,
                    subsample_size=5000,
                    n_threads=n_threads,
                    posterior_levels=posterior_levels
                )
                
                if bio_analysis is not None:
                    bio_gmm = bio_analysis['posterior_thresholds']
                    
                    # Add to multi-panel plot if requested
                    if plot_dir is not None and len(bio_samples) > 1:
                        # Create subplot (2 panels: histogram and KDE)
                        row = i // 3
                        col = i % 3
                        
                        # Each group gets 2 subplots (histogram and KDE)
                        ax1 = plt.subplot2grid((((n_groups + 2) // 3), 6), (row, col*2), colspan=1)
                        ax2 = plt.subplot2grid((((n_groups + 2) // 3), 6), (row, col*2 + 1), colspan=1)
                        
                        # Use the plot_mixture_analysis function with provided axes
                        plot_mixture_analysis(
                            data=bio_analysis['pooled_data'],
                            params=bio_analysis['mixture_params'],
                            output_path=None,  # Not used when axes provided
                            title_prefix=bio_sample,
                            min_umi_threshold=min_umi_threshold,
                            posterior_levels=posterior_levels,
                            axes=(ax1, ax2)
                        )
                else:
                    bio_gmm = {}
                
                # Add biological sample-specific thresholds for this method
                for level in posterior_levels:
                    col_name = f'gmm_{level}_biosample_{method_name}'
                    if level in bio_gmm:
                        threshold_df.loc[bio_mask, col_name] = bio_gmm[level]
                    else:
                        # Initialize if not already present
                        if col_name not in threshold_df.columns:
                            threshold_df[col_name] = np.nan
            
            # Save multi-panel plot if created
            if plot_dir is not None and len(bio_samples) > 1:
                plt.suptitle(f'GMM Thresholds by Biological Sample ({method_name})', fontsize=14)
                plt.tight_layout()
                
                metric_dir = Path(plot_dir) / "guide_mixture_analysis" / method_name / "biological_sample" / "linear"
                metric_dir.mkdir(parents=True, exist_ok=True)
                
                plot_path = metric_dir / "plot.png"
                plt.savefig(plot_path, bbox_inches='tight')
                plt.close()
    
    # Generate individual guide distribution plots if plot_dir is provided
    if plot_dir is not None:
        print("\nGenerating individual guide distribution plots...")
        
        # Get default cell calling method from config
        default_method = config.get('cell_calling', {}).get('default_method', 'BarcodeRanks_Knee')
        
        # Get cells for default method
        default_method_barcodes = cell_barcodes_dict.get(default_method)
        if default_method_barcodes is None:
            print(f"  Warning: Default method {default_method} not found in cell barcodes, skipping plots")
            return threshold_df_filtered
        
        # Filter to default method cells
        method_mask = adata.obs.index.isin(default_method_barcodes)
        adata_default = adata[method_mask]
        
        # Get thresholds for default method cells
        guide_counts = adata_default.obsm['guide_counts']
        
        if isinstance(guide_assignment_method, int):
            # Fixed threshold
            thresholds = guide_assignment_method
        elif guide_assignment_method.startswith('gmm_'):
            # GMM-based threshold for default method
            level = int(guide_assignment_method[4:].split('_')[0])
            threshold_col = f'gmm_{level}_sample_{default_method}'
            thresholds = threshold_df.loc[method_mask, threshold_col].values
        
        # Create binary guide assignment matrix
        guide_assignment = create_binary_guide_assignment(guide_counts, thresholds)
        
        # Store guide_assignment in adata_default.obsm for the plotting function
        adata_default.obsm['guide_assignment'] = guide_assignment
        
        # Call plotting function with only default method cells
        plot_top_guide_distributions(adata_default, plot_dir, default_method, 'sample')
        if 'biological_sample' in adata_default.obs.columns:
            plot_top_guide_distributions(adata_default, plot_dir, default_method, 'biological_sample')
    
    # Filter out rows where ALL values are NaN (cells not called by any method)
    # This significantly reduces file size by removing uncalled cells
    mask = threshold_df.notna().any(axis=1)
    threshold_df_filtered = threshold_df[mask]
    
    print(f"Filtered threshold table from {len(threshold_df):,} to {len(threshold_df_filtered):,} cells (removed uncalled cells)")
    
    return threshold_df_filtered


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Calculate GMM thresholds for guide calling")
    parser.add_argument("--h5ad", required=True, help="Path to h5ad file")
    parser.add_argument("--cell-calling-dir", required=True, help="Path to cell calling results directory")
    parser.add_argument("--sample-id", required=True, help="Sample ID")
    parser.add_argument("--config", required=True, help="Path to config file")
    parser.add_argument("--output", required=True, help="Output path for threshold table")
    parser.add_argument("--plot-dir", help="Directory for plots")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    
    args = parser.parse_args()
    
    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)
    
    # Load adata
    print(f"Loading {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    
    # Load cell barcodes for all methods
    cell_calling_dir = Path(args.cell_calling_dir)
    cell_barcodes_dict = {}
    
    # Find all cell barcode files
    for barcode_file in cell_calling_dir.glob(f"{args.sample_id}_*_cell_barcodes.txt"):
        # Extract method name from filename
        method = barcode_file.stem.replace(f"{args.sample_id}_", "").replace("_cell_barcodes", "")
        
        with open(barcode_file) as f:
            barcodes = set(line.strip() for line in f)
            cell_barcodes_dict[method] = barcodes
            print(f"Loaded {len(barcodes)} cells for {method}")
    
    if not cell_barcodes_dict:
        raise ValueError(f"No cell barcode files found in {cell_calling_dir}")
    
    # Get guide cutoffs and posterior levels from config
    guide_cutoffs = config.get('qc_analysis', {}).get('guide_cutoffs', [])
    posterior_levels = [50]  # Always include 50%
    for cutoff in guide_cutoffs:
        if isinstance(cutoff, str) and cutoff.startswith('gmm_'):
            level = int(cutoff[4:])
            if level not in posterior_levels:
                posterior_levels.append(level)
    
    min_umi_threshold = config.get('guide_mixture_model', {}).get('min_umi_threshold', 2)
    
    # Get required config values - crash if not present
    default_cell_method = config['cell_calling']['default_method']
    default_guide_cutoff = config['qc_analysis']['default_guide_cutoff']
    default_granularity = config['qc_analysis']['default_granularity']
    
    # Build the guide assignment method string
    if isinstance(default_guide_cutoff, str) and default_guide_cutoff.startswith('gmm_'):
        level = default_guide_cutoff[4:]  # Extract the number (e.g., "50" from "gmm_50")
        guide_assignment_method = f"gmm_{level}_{default_granularity}_{default_cell_method}"
    else:
        guide_assignment_method = default_guide_cutoff
    
    # Calculate thresholds and generate plots for all methods
    threshold_table = create_cell_threshold_table(
        adata, guide_cutoffs, posterior_levels, min_umi_threshold, 
        args.threads, cell_barcodes_dict, plot_dir=args.plot_dir,
        guide_assignment_method=guide_assignment_method, config=config
    )
    
    # Save threshold table
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    threshold_table.to_csv(args.output, sep='\t')
    print(f"Saved threshold table to {args.output}")












