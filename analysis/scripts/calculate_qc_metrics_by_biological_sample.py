#!/usr/bin/env python3
"""
Calculate comprehensive QC metrics at three stratification levels.

TODO: Refactor this file - it's getting long (1500+ lines)
- Could split plotting functions into a separate qc_plotting_helpers.py module
- plot_per_cell_distributions() is 400+ lines and could be broken up

This script loads the annotated h5ad file and cell calling results to compute metrics
at different stratification levels based on the --stratify-by parameter.

All output files contain the following metrics for EACH cell calling method:
- n_cells: Number of cells identified
- median/mean_pct_mito_cells: Mitochondrial percentage in cells
- median/mean_umis_per_cell: UMI counts per cell  
- median/mean_genes_per_cell: Genes detected per cell
- total_umis_in_cells: Total UMIs in called cells
- reads_per_umi_in_cells: Sequencing saturation metric (only meaningful at sample level)
- guides_per_cell_cutoff[1-5]: Average guides per cell at different count thresholds
- Plus equivalent metrics for non-cells (background barcodes)

Stratification levels (--stratify-by parameter):
- sample: One row per sequencing sample
  * Includes all read statistics (total_reads, mapped_reads, pcr_duplication_rate_pct, etc.)
  * Includes guide read statistics (prefixed with "guide_")
  * This is the most comprehensive output
- biological_sample: One row per biological replicate (e.g., resting_rep1)
  * Only QC metrics (no read statistics)
  * Shows which wells contribute to each biological sample
- well: One row per plate well (e.g., A1, B2)
  * Only QC metrics (no read statistics)
  * Shows which biological sample is in each well
"""

import argparse
import os
import sys
import gc
from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import yaml
import scipy.sparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
from datetime import datetime
from scipy.stats import median_abs_deviation, poisson, norm, pearsonr

# Set global matplotlib parameters
plt.rcParams['savefig.dpi'] = 100
# Import shared guide utility function
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.guide_analysis import apply_guide_threshold, get_cutoff_suffix, load_guide_reference, calculate_guide_metrics_for_barcodes, get_threshold_for_spec, get_cutoff_labels


# Removed calculate_2d_cell_quality_statistics - moved to generate_qc_cell_lists.py


def plot_violin_with_mean(data, x, y, order, ax, log_scale, show_statistics, metric_type):
    """Helper function to create violin plot with mean markers.
    
    Args:
        data: DataFrame with data to plot
        x: Column name for x-axis grouping
        y: Column name for y-axis values
        order: Order of groups on x-axis
        ax: Matplotlib axis to plot on
        log_scale: Whether to use log scale
        show_statistics: Whether to add statistical lines (only for mitochondrial percentage)
        metric_type: Type of metric being plotted ('pct_mito' for mitochondrial percentage)
        
    Returns:
        dict: Fraction of non-zero values for each group
    """
    # Calculate fraction of non-zero values for each group
    nonzero_fractions = {}
    for group in order:
        group_data = data[data[x] == group][y]
        if len(group_data) > 0:
            nonzero_fractions[group] = (group_data > 0).sum() / len(group_data)
        else:
            nonzero_fractions[group] = 0
    
    # Filter to only non-zero values for plotting
    plot_data = data[data[y] >= 1].copy()
    
    if log_scale:
        ax.set_yscale('log')
    
    # Create violin plot with quartiles, no fill, and cut=0
    # Only plot if there's data for this group
    if len(plot_data) > 0:
        sns.violinplot(data=plot_data, x=x, y=y, order=order, ax=ax, split=True,
                       inner="quart", inner_kws=dict(linewidth=2, color='black'), 
                       fill=False, cut=0)
    
    # Add mean markers (calculated on non-zero data)
    means = plot_data.groupby(x)[y].mean()
    for i, group in enumerate(order):
        if group in means:
            ax.plot(i, means[group], marker='o', color='red', markersize=8, 
                   markeredgecolor='darkred', markeredgewidth=2)
    
    # Update x-tick labels to include non-zero fraction
    labels = []
    for group in order:
        frac = nonzero_fractions.get(group, 0)
        labels.append(f"{group}\n({frac:.0%} non-zero)")
    ax.set_xticklabels(labels)
    
    # Add statistical lines if requested and this is mitochondrial percentage
    if show_statistics and metric_type == 'pct_mito' and len(plot_data) > 0:
        # Calculate simple statistics on all data
        all_values = data[y].values
        median = np.median(all_values)
        mad = median_abs_deviation(all_values)
        
        # MAD-based (solid lines)
        ax.axhline(y=median + 2*mad, color='cyan', linestyle='-', linewidth=1,
                  label=f"Median+2×MAD ({median + 2*mad:.1f}%)")
        ax.axhline(y=median + 3*mad, color='navy', linestyle='-', linewidth=1,
                  label=f"Median+3×MAD ({median + 3*mad:.1f}%)")
        
        # Percentiles (dashed lines)
        p75 = np.percentile(all_values, 75)
        p95 = np.percentile(all_values, 95)
        p99 = np.percentile(all_values, 99)
        
        ax.axhline(y=p75, color='green', linestyle='--', linewidth=1,
                  label=f"75th percentile ({p75:.1f}%)")
        ax.axhline(y=p95, color='orange', linestyle='--', linewidth=1,
                  label=f"95th percentile ({p95:.1f}%)")
        ax.axhline(y=p99, color='red', linestyle='--', linewidth=1,
                  label=f"99th percentile ({p99:.1f}%)")
        
        # Add legend
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    
    return nonzero_fractions


# Removed save_qc_cell_lists - moved to generate_qc_cell_lists.py

def load_cell_barcodes(cell_calling_dir, sample_id):
    """Load cell barcodes for each cell calling method."""
    cell_calling_dir = Path(cell_calling_dir)
    
    # Find all cell barcode files using the full sample_id
    methods = []
    cell_barcodes = {}
    
    for barcode_file in cell_calling_dir.glob(f"{sample_id}_*_cell_barcodes.txt"):
        # Extract method name from filename
        # Format: {sample_id}_{method}_cell_barcodes.txt
        method = barcode_file.stem.replace(f"{sample_id}_", "").replace("_cell_barcodes", "")
        methods.append(method)
        
        with open(barcode_file, 'r') as f:
            barcodes = [line.strip() for line in f if line.strip()]
            cell_barcodes[method] = set(barcodes)
            print(f"Loaded {len(barcodes)} cell barcodes for {method}")
    
    return cell_barcodes


# load_guide_reference moved to scripts/guide_analysis.py


def calculate_metrics_for_barcodes(adata, barcode_set, is_cell=True):
    """Calculate QC metrics for a given set of barcodes."""
    # Create mask for the barcode set
    if barcode_set:
        mask = adata.obs.index.isin(barcode_set)
    else:
        # Empty set - no cells
        mask = np.zeros(adata.n_obs, dtype=bool)
    
    if not is_cell:
        # For non-cells, invert the mask
        mask = ~mask
    
    n_barcodes = np.sum(mask)
    
    if n_barcodes == 0:
        return {
            'n': 0,
            'median_pct_mito': np.nan,
            'mean_pct_mito': np.nan,
            'median_umis': np.nan,
            'mean_umis': np.nan,
            'median_genes': np.nan,
            'mean_genes': np.nan,
            'total_umis': 0
        }
    
    # Subset data
    subset = adata[mask]
    
    # Calculate metrics
    metrics = {
        'n': n_barcodes,
        'median_pct_mito': float(np.median(subset.obs['pct_counts_mt'])),
        'mean_pct_mito': float(np.mean(subset.obs['pct_counts_mt'])),
        'median_umis': float(np.median(subset.obs['total_counts'])),
        'mean_umis': float(np.mean(subset.obs['total_counts'])),
        'median_genes': float(np.median(subset.obs['n_genes_by_counts'])),
        'mean_genes': float(np.mean(subset.obs['n_genes_by_counts'])),
        'total_umis': int(subset.obs['total_counts'].sum())
    }
    
    return metrics


# calculate_guide_metrics_for_barcodes moved to scripts/guide_analysis.py


def calculate_cell_specific_guide_overlap(adata, cell_barcodes):
    """Calculate what percentage of total guide UMIs are captured in called cells.
    
    This is a cell-specific version of the overlap statistics, calculated AFTER
    cell calling to show the true guide capture efficiency in real cells only.
    
    Args:
        adata: Full AnnData object with guide_counts in obsm
        cell_barcodes: Set of cell barcodes from cell calling
        
    Returns:
        float: Percentage of total guide UMIs found in called cells
    """
    if 'guide_counts' not in adata.obsm:
        return np.nan
    
    # Get total guide UMIs across ALL barcodes
    guide_counts = adata.obsm['guide_counts']
    total_guide_umis = int(guide_counts.sum())
    
    if total_guide_umis == 0:
        return 0.0
    
    # Get guide UMIs in cells only
    cell_mask = adata.obs.index.isin(cell_barcodes)
    guide_umis_in_cells = int(guide_counts[cell_mask].sum())
    
    # Calculate percentage
    percentage = (guide_umis_in_cells / total_guide_umis) * 100
    
    return percentage


def calculate_metrics_for_group(adata, cell_barcodes, group_mask, total_reads, group_name, guide_to_genes=None, guide_cutoffs=None, guide_total_reads=None, mixture_metrics=None, qc_cell_lists_df=None, default_method=None, stratify_by='biological_sample'):
    """Calculate all metrics for a specific group (sample/biological sample/well).
    
    Args:
        adata: AnnData object
        cell_barcodes: Dict of cell barcodes by method
        group_mask: Boolean mask for this group
        total_reads: Total reads for the sample
        group_name: Name of the group (e.g., 'Overall', biological sample name, well name)
        guide_to_genes: Dictionary mapping guide names to gene names
        guide_cutoffs: List of guide count cutoffs
        guide_total_reads: Total guide reads (only for sample level)
        mixture_metrics: Dict of mixture model metrics by method
    """
    # Initialize results dictionary
    results = {}
    
    # Calculate metrics for each cell calling method
    for method, method_barcodes in cell_barcodes.items():
        # Get cells that are both in this method AND this group
        method_mask = adata.obs.index.isin(method_barcodes)
        combined_mask = method_mask & group_mask
        group_barcodes = set(adata.obs.index[combined_mask])
        
        # Calculate metrics for cells
        cell_metrics = calculate_metrics_for_barcodes(adata[group_mask], group_barcodes, is_cell=True)
        
        # Calculate metrics for non-cells
        noncell_metrics = calculate_metrics_for_barcodes(adata[group_mask], group_barcodes, is_cell=False)
        
        # Calculate guide metrics for cells
        guide_metrics = calculate_guide_metrics_for_barcodes(adata[group_mask], group_barcodes, guide_to_genes, guide_cutoffs, stratify_by, method)
        
        # Calculate reads per UMI for cells
        # Note: For biological samples, we can't calculate a meaningful reads_per_umi
        # because we don't know what fraction of total reads came from each bio sample
        # Only meaningful for Overall
        if group_name == 'Overall':
            if cell_metrics['total_umis'] == 0:
                reads_per_umi = 0.0  # No UMIs found, set to 0
            else:
                reads_per_umi = total_reads / cell_metrics['total_umis']
        
        # Add all metrics directly to results
        
        # Cell calling performance
        results[f'n_cells_{method}'] = cell_metrics['n']
        results[f'total_umis_in_cells_{method}'] = cell_metrics['total_umis']
        if group_name == 'Overall':
            results[f'reads_per_umi_in_cells_{method}'] = reads_per_umi
            
            # Calculate reads per cell metrics (only valid at sample level)
            n_cells = cell_metrics['n']
            if n_cells > 0:
                results[f'total_reads_per_cell_{method}'] = total_reads / n_cells
                # Add guide reads per cell if guide reads are available
                if guide_total_reads is not None:
                    results[f'guide_reads_per_cell_{method}'] = guide_total_reads / n_cells
            
            # Calculate percentage of UMIs in cells (only valid at sample level)
            # Need total UMIs across ALL barcodes (cells + background)
            total_umis_all_barcodes = adata[group_mask].X.sum()
            if total_umis_all_barcodes > 0:
                results[f'umis_in_cells_pct_{method}'] = cell_metrics['total_umis'] / total_umis_all_barcodes * 100
            else:
                results[f'umis_in_cells_pct_{method}'] = 0.0
        
        # Cell quality
        results[f'median_pct_mito_cells_{method}'] = cell_metrics['median_pct_mito']
        results[f'mean_pct_mito_cells_{method}'] = cell_metrics['mean_pct_mito']
        results[f'median_umis_per_cell_{method}'] = cell_metrics['median_umis']
        results[f'mean_umis_per_cell_{method}'] = cell_metrics['mean_umis']
        results[f'median_genes_per_cell_{method}'] = cell_metrics['median_genes']
        results[f'mean_genes_per_cell_{method}'] = cell_metrics['mean_genes']
        
        # Background quality
        results[f'median_pct_mito_noncells_{method}'] = noncell_metrics['median_pct_mito']
        results[f'mean_pct_mito_noncells_{method}'] = noncell_metrics['mean_pct_mito']
        results[f'median_umis_per_noncell_{method}'] = noncell_metrics['median_umis']
        results[f'mean_umis_per_noncell_{method}'] = noncell_metrics['mean_umis']
        results[f'median_genes_per_noncell_{method}'] = noncell_metrics['median_genes']
        results[f'mean_genes_per_noncell_{method}'] = noncell_metrics['mean_genes']
        
        # Guide capture
        results[f'guide_umis_per_cell_{method}'] = guide_metrics['guide_umis_per_cell']
        
        # Calculate cell-specific guide overlap (only for sample level where we have full data)
        if group_name == 'Overall':
            guide_overlap_pct = calculate_cell_specific_guide_overlap(adata[group_mask], group_barcodes)
            results[f'guide_umis_in_cells_pct_{method}'] = guide_overlap_pct
        
        for cutoff_spec in guide_cutoffs:
            suffix = get_cutoff_suffix(cutoff_spec)
            results[f'guides_per_cell_{suffix}_{method}'] = guide_metrics[f'guides_per_cell_{suffix}']
            results[f'fraction_cells_with_guides_{suffix}_{method}'] = guide_metrics[f'fraction_cells_with_guides_{suffix}']
            results[f'umis_per_guide_per_cell_{suffix}_{method}'] = guide_metrics[f'umis_per_guide_per_cell_{suffix}']
            
            # Guide diversity
            results[f'total_guides_detected_{suffix}_{method}'] = guide_metrics[f'total_guides_detected_{suffix}']
            results[f'total_guide_detections_{suffix}_{method}'] = guide_metrics[f'total_guide_detections_{suffix}']
            results[f'total_targeted_genes_{suffix}_{method}'] = guide_metrics[f'total_targeted_genes_{suffix}']
        
        # Add mixture model metrics if available
        if mixture_metrics is not None:
            metrics = None
            if method in mixture_metrics:
                metrics = mixture_metrics[method]
            elif f"{method}_{group_name}" in mixture_metrics:
                metrics = mixture_metrics[f"{method}_{group_name}"]
            
            if metrics:
                # Add posterior thresholds for all levels
                posterior_thresholds = metrics.get('posterior_thresholds', {})
                for level, threshold in posterior_thresholds.items():
                    results[f'guide_posterior_threshold_{level}pct_{method}'] = threshold
                
                # Add other mixture metrics
                results[f'guide_mixture_mean_log_likelihood_{method}'] = metrics.get('mean_log_likelihood', None)
                results[f'guide_mixture_lambda_{method}'] = metrics.get('lambda', None)
                results[f'guide_mixture_weight_poisson_{method}'] = metrics.get('weight_poisson', None)
                results[f'guide_mixture_gaussian_mean_{method}'] = metrics.get('gaussian_mean_umis', None)
                results[f'guide_mixture_gaussian_median_{method}'] = metrics.get('gaussian_median_umis', None)
                results[f'guide_mixture_fraction_gaussian_50pct_{method}'] = metrics.get('fraction_gaussian_50pct', None)
    
    # Add percentage of compromised cells for default method
    if qc_cell_lists_df is not None and default_method is not None:
        prob_compromised = qc_cell_lists_df['posterior_prob_compromised']
        pct_compromised = (prob_compromised > 0.50).mean() * 100
        results[f'pct_compromised_cells_50_{default_method}'] = pct_compromised
    
    # Add GMM threshold values for each method
    # The thresholds are the same for all cells in the group, so we just need to get one value
    for method, method_barcodes in cell_barcodes.items():
        # Get cells that are both in this method AND this group
        method_mask = adata.obs.index.isin(method_barcodes)
        combined_mask = method_mask & group_mask
        
        if np.sum(combined_mask) > 0:
            # Get the first cell's thresholds (all cells in group have same value)
            first_cell_idx = np.where(combined_mask)[0][0]
            
            # Check for GMM threshold columns
            for col in adata.obs.columns:
                if col.startswith('gmm_') and (col.endswith('_sample') or col.endswith('_biosample')):
                    # Get the threshold value for this group
                    threshold_value = adata.obs.iloc[first_cell_idx][col]
                    
                    # Only add if it's a valid threshold (not -1)
                    if threshold_value > 0:
                        # Add as a metric with the method suffix
                        # e.g., gmm_50_sample_threshold_BarcodeRanks_Knee
                        results[f'{col}_threshold_{method}'] = threshold_value
    
    return results


    if analysis_result is not None:
        return analysis_result['posterior_thresholds']
    return {}


# get_threshold_for_spec and get_cutoff_labels moved to scripts/guide_analysis.py


def plot_per_cell_distributions(adata, cell_barcodes, groups, stratify_by, output_dir, sample_id, source, processing, guide_cutoffs, pool, read_stats, plot_method_filter):
    """Generate violin plots for per-cell metrics as individual PNG files.
    
    Note: Assumes GMM threshold columns are already in adata.obs (e.g., gmm_50_sample, gmm_75_biosample)
    
    Returns:
        dict: Mixture model metrics for each method and group
    """
    # Initialize dictionary to store mixture metrics
    mixture_metrics = {}
    
    # Skip plotting for wells (too many groups)
    if stratify_by == 'well':
        print("Skipping per-cell distribution plots for well-level analysis (too many groups)")
        return mixture_metrics
    
    # Use provided plot directory
    plot_dir = Path(output_dir)
    plot_dir.mkdir(exist_ok=True)
    
    print(f"\nGenerating per-cell distribution plots in {plot_dir}")
    
    # Set plot style
    sns.set_style("whitegrid")
    plt.rcParams['figure.dpi'] = 100  # For display
    
    # Define metrics to plot
    metrics = [
        ('total_counts', 'UMIs per cell', 'umis_per_cell'),
        ('n_genes_by_counts', 'Genes per cell', 'genes_per_cell'),
        ('pct_counts_mt', 'Mitochondrial %', 'pct_mito')
    ]
    
    # Check if guide data exists
    has_guides = 'guide_counts' in adata.obsm
    # guide_cutoffs is always required, caller must handle this
    
    # Determine which methods to plot
    if plot_method_filter == "all":
        methods_to_plot = cell_barcodes.items()
        print(f"  Plotting all {len(cell_barcodes)} methods")
    elif plot_method_filter in cell_barcodes:
        methods_to_plot = [(plot_method_filter, cell_barcodes[plot_method_filter])]
        print(f"  Plotting only method: {plot_method_filter}")
    else:
        raise ValueError(f"Method '{plot_method_filter}' not found. Available: {', '.join(cell_barcodes.keys())}")
    
    # Process each cell calling method
    for method, method_barcodes in methods_to_plot:
        print(f"  Plotting distributions for {method}...")
        
        # Skip if no cells for this method
        if len(method_barcodes) == 0:
            print(f"    No cells found for {method}, skipping")
            continue
        
        # Collect data for all groups
        group_data = []
        group_labels = []
        
        for group in groups:
            # Create mask for this group
            if stratify_by == 'sample' or group == 'Overall':
                group_mask = np.ones(adata.n_obs, dtype=bool)
            else:
                stratify_column = 'biological_sample' if stratify_by == 'biological_sample' else 'well'
                group_mask = adata.obs[stratify_column] == group
            
            # Get cells that are both in this method AND this group
            method_mask = adata.obs.index.isin(method_barcodes)
            combined_mask = method_mask & group_mask
            
            if np.sum(combined_mask) > 0:
                group_data.append(adata[combined_mask])
                group_labels.append(group)
        
        # Skip if no data for any group
        if len(group_data) == 0:
            print(f"    No data found for {method} in any group, skipping")
            continue
        
        # Calculate mean values for sorting
        group_means = {}
        for i, (data, label) in enumerate(zip(group_data, group_labels)):
            group_means[label] = {
                'total_counts': np.mean(data.obs['total_counts']),
                'n_genes_by_counts': np.mean(data.obs['n_genes_by_counts']),
                'pct_counts_mt': np.mean(data.obs['pct_counts_mt'])
            }
        
        # Plot each metric
        for obs_col, metric_label, file_suffix in metrics:
            # Sort groups by mean value for this metric (highest to lowest)
            sorted_groups = sorted(group_labels, key=lambda x: group_means[x][obs_col], reverse=True)
            
            # Prepare data for plotting
            plot_data = []
            plot_labels = []
            
            for group in sorted_groups:
                idx = group_labels.index(group)
                values = group_data[idx].obs[obs_col].values
                plot_data.extend(values)
                plot_labels.extend([group] * len(values))
            
            # Create DataFrame for seaborn
            plot_df = pd.DataFrame({
                'value': plot_data,
                'group': plot_labels
            })
            
            # Apply percentile cutoffs (1% and 99%)
            percentile_1 = np.percentile(plot_df['value'], 1)
            percentile_99 = np.percentile(plot_df['value'], 99)
            plot_df['value'] = plot_df['value'].clip(percentile_1, percentile_99)
            
            # Create two plots: linear and log scale
            for scale_type in ['linear', 'log']:
                fig, ax = plt.subplots(figsize=(max(8, len(sorted_groups) * 0.8), 6))
                
                # Make violin plot
                # Use helper function for violin plot
                show_stats = (file_suffix == 'pct_mito')  # Only show statistics for mitochondrial percentage
                plot_violin_with_mean(plot_df, 'group', 'value', sorted_groups, ax, 
                                    log_scale=(scale_type == 'log'), 
                                    show_statistics=show_stats, 
                                    metric_type=file_suffix)
            
                # Customize plot
                ax.set_xlabel(stratify_by.replace('_', ' ').title())
                ax.set_ylabel(metric_label + (' (log scale)' if scale_type == 'log' else ''))
                
                # Add read count and cell count to title
                n_cells = len(method_barcodes)
                if read_stats and 'total_reads' in read_stats:
                    total_reads = read_stats['total_reads']
                    reads_per_cell = total_reads // n_cells if n_cells > 0 else 0
                    title = f'{metric_label} Distribution by {stratify_by.replace("_", " ").title()} ({method}) - {scale_type.capitalize()} Scale\nCells: {n_cells:,} | Reads/cell: {reads_per_cell:,} | Percentiles: 1%-99%'
                else:
                    title = f'{metric_label} Distribution by {stratify_by.replace("_", " ").title()} ({method}) - {scale_type.capitalize()} Scale\nCells: {n_cells:,} | Percentiles: 1%-99%'
                ax.set_title(title)
                
                # Rotate x labels if many groups
                if len(sorted_groups) > 5:
                    plt.xticks(rotation=45, ha='right')
                
                # Adjust layout
                plt.tight_layout()
                
                # Save plot with category subdirectory
                # Create clean directory structure
                # per_cell/{source}_{processing}/{pool}/{sample}/{file_suffix}/{method}/{stratify_by}/plot_{scale_type}.png
                if pool:
                    metric_dir = plot_dir / file_suffix / method / stratify_by / scale_type
                else:
                    metric_dir = plot_dir / file_suffix / method / stratify_by / scale_type
                metric_dir.mkdir(parents=True, exist_ok=True)
                
                filename = "plot.png"
                filepath = metric_dir / filename
                plt.savefig(filepath, bbox_inches='tight')
                plt.close()
                
                
                print(f"    Saved {filename} ({scale_type} scale)")
        
        # Plot guide metrics if available
        if has_guides:
            # Iterate through each cutoff specification from the config
            for cutoff_spec in guide_cutoffs:
                # Get file suffix and plot label
                cutoff_suffix, plot_title_suffix = get_cutoff_labels(cutoff_spec)
                
                # Calculate guides per cell for this cutoff specification
                guides_per_cell_data = []
                guides_per_cell_labels = []
                guides_group_means = {}
                
                for i, (data, label) in enumerate(zip(group_data, group_labels)):
                    # Get guide counts for these cells
                    guide_counts = data.obsm['guide_counts']
                    
                    # Apply threshold using helper function
                    binary_guides = apply_guide_threshold(guide_counts, cutoff_spec, data.obs, stratify_by, method)
                    
                    # Sum guides per cell
                    guides_per_cell = np.array(binary_guides.sum(axis=1)).flatten()
                    
                    guides_per_cell_data.extend(guides_per_cell)
                    guides_per_cell_labels.extend([label] * len(guides_per_cell))
                    guides_group_means[label] = np.mean(guides_per_cell)
                
                # Sort groups by mean guides per cell (highest to lowest)
                sorted_groups = sorted(group_labels, key=lambda x: guides_group_means[x], reverse=True)
                
                # Create DataFrame for seaborn
                plot_df = pd.DataFrame({
                    'value': guides_per_cell_data,
                    'group': guides_per_cell_labels
                })
                
                # Apply percentile cutoffs (1% and 99%)
                percentile_1 = np.percentile(plot_df['value'], 1)
                percentile_99 = np.percentile(plot_df['value'], 99)
                plot_df['value'] = plot_df['value'].clip(percentile_1, percentile_99)
                
                # Create two plots: linear and log scale
                for scale_type in ['linear', 'log']:
                    fig, ax = plt.subplots(figsize=(max(8, len(sorted_groups) * 0.8), 6))
                    
                    # Make violin plot using helper function for consistency
                    plot_violin_with_mean(plot_df, 'group', 'value', sorted_groups, ax, log_scale=(scale_type == 'log'), show_statistics=False, metric_type=None)
                
                    # Customize plot
                    ax.set_xlabel(stratify_by.replace('_', ' ').title())
                    ax.set_ylabel(f'Guides per cell ({plot_title_suffix})' + (' (log scale)' if scale_type == 'log' else ''))
                    
                    # Add guide read count and cell count to title
                    n_cells = len(method_barcodes)
                    if read_stats and 'guide_total_reads' in read_stats:
                        guide_total_reads = read_stats['guide_total_reads']
                        guide_reads_per_cell = guide_total_reads // n_cells if n_cells > 0 else 0
                        title = f'Guides per Cell Distribution ({plot_title_suffix}) by {stratify_by.replace("_", " ").title()} ({method}) - {scale_type.capitalize()} Scale\nCells: {n_cells:,} | Guide reads/cell: {guide_reads_per_cell:,} | Percentiles: 1%-99%'
                    else:
                        title = f'Guides per Cell Distribution ({plot_title_suffix}) by {stratify_by.replace("_", " ").title()} ({method}) - {scale_type.capitalize()} Scale\nCells: {n_cells:,} | Percentiles: 1%-99%'
                    ax.set_title(title)
                    
                    # Rotate x labels if many groups
                    if len(sorted_groups) > 5:
                        plt.xticks(rotation=45, ha='right')
                    
                    # Adjust layout
                    plt.tight_layout()
                    
                    # Save plot with clean directory structure
                    # per_cell/{source}_{processing}/{pool}/{sample}/guides_per_cell_{cutoff_suffix}/{method}/{stratify_by}/plot_{scale_type}.png
                    metric_name = f"guides_per_cell_{cutoff_suffix}"
                    metric_dir = plot_dir / metric_name / method / stratify_by / scale_type
                    metric_dir.mkdir(parents=True, exist_ok=True)
                    
                    filename = "plot.png"
                    filepath = metric_dir / filename
                    plt.savefig(filepath, bbox_inches='tight')
                    plt.close()
                    
                    
                    print(f"    Saved {filename} ({scale_type} scale)")
            
            # Plot guide UMIs per cell (cutoff-independent)
            guide_umis_data = []
            guide_umis_labels = []
            guide_umis_means = {}
            
            for i, (data, label) in enumerate(zip(group_data, group_labels)):
                # Get total guide UMIs per cell
                guide_counts = data.obsm['guide_counts']
                total_umis_per_cell = np.array(guide_counts.sum(axis=1)).flatten()
                
                guide_umis_data.extend(total_umis_per_cell)
                guide_umis_labels.extend([label] * len(total_umis_per_cell))
                guide_umis_means[label] = np.mean(total_umis_per_cell)
            
            # Sort groups by mean guide UMIs per cell (highest to lowest)
            sorted_groups = sorted(group_labels, key=lambda x: guide_umis_means[x], reverse=True)
            
            # Create DataFrame for seaborn
            plot_df = pd.DataFrame({
                'value': guide_umis_data,
                'group': guide_umis_labels
            })
            
            # Apply percentile cutoffs (1% and 99%)
            percentile_1 = np.percentile(plot_df['value'], 1)
            percentile_99 = np.percentile(plot_df['value'], 99)
            plot_df['value'] = plot_df['value'].clip(percentile_1, percentile_99)
            
            # Create two plots: linear and log scale
            for scale_type in ['linear', 'log']:
                fig, ax = plt.subplots(figsize=(max(8, len(sorted_groups) * 0.8), 6))
                
                # Make violin plot
                # Use helper function for violin plot
                plot_violin_with_mean(plot_df, 'group', 'value', sorted_groups, ax, log_scale=(scale_type == 'log'), show_statistics=False, metric_type=None)
            
                # Customize plot
                ax.set_xlabel(stratify_by.replace('_', ' ').title())
                ax.set_ylabel('Guide UMIs per cell' + (' (log scale)' if scale_type == 'log' else ''))
                
                # Add guide read count and cell count to title
                n_cells = len(method_barcodes)
                if read_stats and 'guide_total_reads' in read_stats:
                    guide_total_reads = read_stats['guide_total_reads']
                    guide_reads_per_cell = guide_total_reads // n_cells if n_cells > 0 else 0
                    title = f'Guide UMIs per Cell Distribution by {stratify_by.replace("_", " ").title()} ({method}) - {scale_type.capitalize()} Scale\nCells: {n_cells:,} | Guide reads/cell: {guide_reads_per_cell:,} | Percentiles: 1%-99%'
                else:
                    title = f'Guide UMIs per Cell Distribution by {stratify_by.replace("_", " ").title()} ({method}) - {scale_type.capitalize()} Scale\nCells: {n_cells:,} | Percentiles: 1%-99%'
                ax.set_title(title)
                
                # Rotate x labels if many groups
                if len(sorted_groups) > 5:
                    plt.xticks(rotation=45, ha='right')
                
                # Adjust layout
                plt.tight_layout()
                
                # Save plot with clean directory structure
                metric_name = "guide_umis_per_cell"
                metric_dir = plot_dir / metric_name / method / stratify_by / scale_type
                metric_dir.mkdir(parents=True, exist_ok=True)
                
                filename = "plot.png"
                filepath = metric_dir / filename
                plt.savefig(filepath, bbox_inches='tight')
                plt.close()
                
                
                print(f"    Saved {filename} ({scale_type} scale)")
            
            # Plot UMIs per guide per cell for each cutoff
            for cutoff_spec in guide_cutoffs:
                # Get file suffix and plot label
                cutoff_suffix, plot_title_suffix = get_cutoff_labels(cutoff_spec)
                
                umis_per_guide_data = []
                umis_per_guide_labels = []
                umis_per_guide_means = {}
                
                for i, (data, label) in enumerate(zip(group_data, group_labels)):
                    # Get guide counts for these cells
                    guide_counts = data.obsm['guide_counts']
                    
                    # Apply threshold using helper function
                    binary_guides = apply_guide_threshold(guide_counts, cutoff_spec, data.obs, stratify_by, method)
                    
                    # Sum guides per cell
                    guides_per_cell = np.array(binary_guides.sum(axis=1)).flatten()
                    
                    # Initialize array for ALL cells (including those without guides)
                    umis_per_guide_all_cells = np.zeros(len(guides_per_cell))
                    
                    # Only calculate for cells with at least one guide
                    cells_with_guides = guides_per_cell > 0
                    if np.any(cells_with_guides):
                        # Only consider UMIs from guides that pass the cutoff
                        guide_counts_filtered = guide_counts[cells_with_guides].multiply(binary_guides[cells_with_guides])
                        umis_per_cell_with_guides = np.array(guide_counts_filtered.sum(axis=1)).flatten()
                        guides_per_cell_with_guides = guides_per_cell[cells_with_guides]
                        # Calculate UMIs per guide for each cell
                        umis_per_guide = umis_per_cell_with_guides / guides_per_cell_with_guides
                        
                        # Place values in the full array
                        umis_per_guide_all_cells[cells_with_guides] = umis_per_guide
                        
                        # Calculate mean only for cells with guides (for sorting)
                        umis_per_guide_means[label] = np.mean(umis_per_guide)
                    
                    # Add ALL cells to the data (including zeros)
                    umis_per_guide_data.extend(umis_per_guide_all_cells)
                    umis_per_guide_labels.extend([label] * len(umis_per_guide_all_cells))
                
                # Only plot if we have data
                if len(umis_per_guide_data) > 0:
                    # Sort groups by mean UMIs per guide (highest to lowest)
                    sorted_groups_umis = sorted([g for g in group_labels if g in umis_per_guide_means], 
                                               key=lambda x: umis_per_guide_means[x], reverse=True)
                    
                    # Create DataFrame for seaborn
                    plot_df = pd.DataFrame({
                        'value': umis_per_guide_data,
                        'group': umis_per_guide_labels
                    })
                    
                    # Apply percentile cutoffs (1% and 99%)
                    percentile_1 = np.percentile(plot_df['value'], 1)
                    percentile_99 = np.percentile(plot_df['value'], 99)
                    plot_df['value'] = plot_df['value'].clip(percentile_1, percentile_99)
                    
                    # Create two plots: linear and log scale
                    for scale_type in ['linear', 'log']:
                        fig, ax = plt.subplots(figsize=(max(8, len(sorted_groups_umis) * 0.8), 6))
                        
                        # Make violin plot
                        # Use helper function for violin plot
                        plot_violin_with_mean(plot_df, 'group', 'value', sorted_groups_umis, ax, log_scale=(scale_type == 'log'), show_statistics=False, metric_type=None)
                    
                        # Customize plot
                        ax.set_xlabel(stratify_by.replace('_', ' ').title())
                        ax.set_ylabel(f'UMIs per guide per cell ({plot_title_suffix})' + (' (log scale)' if scale_type == 'log' else ''))
                        
                        # Add guide read count and cell count to title
                        n_cells = len(method_barcodes)
                        if read_stats and 'guide_total_reads' in read_stats:
                            guide_total_reads = read_stats['guide_total_reads']
                            guide_reads_per_cell = guide_total_reads // n_cells if n_cells > 0 else 0
                            title = f'UMIs per Guide per Cell Distribution ({plot_title_suffix}) by {stratify_by.replace("_", " ").title()} ({method}) - {scale_type.capitalize()} Scale\nCells: {n_cells:,} | Guide reads/cell: {guide_reads_per_cell:,} | Percentiles: 1%-99%'
                        else:
                            title = f'UMIs per Guide per Cell Distribution ({plot_title_suffix}) by {stratify_by.replace("_", " ").title()} ({method}) - {scale_type.capitalize()} Scale\nCells: {n_cells:,} | Percentiles: 1%-99%'
                        ax.set_title(title)
                        
                        # Rotate x labels if many groups
                        if len(sorted_groups_umis) > 5:
                            plt.xticks(rotation=45, ha='right')
                        
                        # Adjust layout
                        plt.tight_layout()
                        
                        # Save plot with category subdirectory
                        # Save plot with clean directory structure
                        metric_name = f"umis_per_guide_per_cell_{cutoff_suffix}"
                        metric_dir = plot_dir / metric_name / method / stratify_by / scale_type
                        metric_dir.mkdir(parents=True, exist_ok=True)
                        
                        filename = "plot.png"
                        filepath = metric_dir / filename
                        plt.savefig(filepath, bbox_inches='tight')
                        plt.close()
                        
                        
                        print(f"    Saved {filename} ({scale_type} scale)")
    
    print("Per-cell distribution plots completed")
    return mixture_metrics


def _plot_umi_scatter_generic(adata, cell_barcodes, sample_id, stratify_by, plot_dir, 
                              x_values_func, x_label,
                              y_values_func, y_label,
                              color_values_func, color_label, color_map, color_vmin, color_vmax,
                              metric_name, plot_title_suffix):
    """Generic function to create scatter plots with custom coloring.
    
    Args:
        adata: AnnData object
        cell_barcodes: Dictionary of cell barcodes by method
        sample_id: Sample identifier
        stratify_by: How to stratify the data
        plot_dir: Directory to save plots
        x_values_func: Function to calculate x-axis values from data subset
        x_label: Label for x-axis
        y_values_func: Function to calculate y-axis values from data subset
        y_label: Label for y-axis
        color_values_func: Function to calculate color values from data subset
        color_label: Label for colorbar
        color_map: Matplotlib colormap name
        color_vmin: Minimum value for color scale (None for auto)
        color_vmax: Maximum value for color scale (None for auto)
        metric_name: Name for output directory
        plot_title_suffix: Additional text for plot title
    """
    
    # Ensure plot directory exists
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(exist_ok=True, parents=True)
    
    # Use all available methods from cell_barcodes
    methods = sorted(cell_barcodes.keys())
    
    for method in methods:
        # Get cells for this method
        cells = cell_barcodes[method]
        if len(cells) == 0:
            print(f"  Skipping {method} - 0 cells")
            continue
        
        # Subset adata to these cells
        cell_mask = adata.obs.index.isin(cells)
        
        # Get data based on stratification
        if stratify_by == 'sample':
            # Use all cells for this method
            group_data = [adata[cell_mask]]
            group_labels = ['All cells']
        elif stratify_by == 'biological_sample':
            # Group by biological sample
            groups = adata.obs.loc[cell_mask, 'biological_sample'].dropna().unique().tolist()
            group_data = []
            group_labels = []
            for group in groups:
                group_mask = cell_mask & (adata.obs['biological_sample'] == group)
                if np.sum(group_mask) > 0:
                    group_data.append(adata[group_mask])
                    group_labels.append(group)
        else:  # well
            # Skip well-level scatter plots (too many plots)
            continue
        
        if len(group_data) == 0:
            continue
        
        # Create figure
        n_groups = len(group_data)
        fig_width = max(8, 5 * min(n_groups, 3))
        fig_height = 5 * ((n_groups + 2) // 3)
        
        fig, axes = plt.subplots(
            nrows=(n_groups + 2) // 3,
            ncols=min(n_groups, 3),
            figsize=(fig_width, fig_height),
            squeeze=False
        )
        axes = axes.flatten()
        
        # Plot each group
        for i, (data, label) in enumerate(zip(group_data, group_labels)):
            ax = axes[i]
            
            # Get the data - always pass method for consistency
            x_values = x_values_func(data, method)
            y_values = y_values_func(data, method)
            color_values = color_values_func(data, method)
            
            # Sample if too many points
            max_points = 10000
            if len(x_values) > max_points:
                # Random sample while preserving distribution
                sample_idx = np.random.choice(len(x_values), size=max_points, replace=False)
                x_values_plot = x_values[sample_idx]
                y_values_plot = y_values[sample_idx]
                color_values_plot = color_values[sample_idx]
                plot_label = f"(showing {max_points:,} of {len(x_values):,} cells)"
            else:
                x_values_plot = x_values
                y_values_plot = y_values
                color_values_plot = color_values
                plot_label = f"({len(x_values):,} cells)"
            
            # Determine color scale limits
            vmin = color_vmin if color_vmin is not None else np.min(color_values)
            vmax = color_vmax if color_vmax is not None else np.max(color_values)
            
            # Create scatter plot
            scatter = ax.scatter(
                x_values_plot,
                y_values_plot,
                c=color_values_plot,
                cmap=color_map,
                s=1,
                alpha=0.5,
                vmin=vmin,
                vmax=vmax
            )
            
            # Set log scale
            ax.set_xscale('log')
            ax.set_yscale('log')
            
            # Labels and title
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.set_title(f'{label}\n{plot_label}')
            
            # Add colorbar
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label(color_label, rotation=270, labelpad=15)
            
            # Calculate and add correlation
            if len(x_values) > 1:
                # Calculate on log-transformed data
                log_x = np.log10(x_values + 1)
                log_y = np.log10(y_values + 1)
                corr, _ = pearsonr(log_x, log_y)
                ax.text(0.05, 0.95, f'r = {corr:.3f}', 
                       transform=ax.transAxes, fontsize=10,
                       verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            # Grid
            ax.grid(True, alpha=0.3)
        
        # Hide empty subplots
        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)
        
        # Overall title
        title = f'{x_label} vs {y_label} (colored by {color_label})\n{method} - {sample_id}'
        if plot_title_suffix:
            title += f' - {plot_title_suffix}'
        fig.suptitle(title, fontsize=14)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save with clean directory structure (no separate correlation subdirectory)
        metric_dir = plot_dir / metric_name / method / stratify_by / "linear"
        metric_dir.mkdir(parents=True, exist_ok=True)
        
        filename = "plot.png"
        filepath = metric_dir / filename
        plt.savefig(filepath, bbox_inches='tight')
        plt.close()
        
        print(f"    Saved {metric_name}/{method}/{stratify_by}/linear/{filename}")


def plot_umi_vs_genes_with_cell_quality(adata, cell_barcodes, sample_id, stratify_by, plot_dir, qc_cell_lists_path=None):
    """Create side-by-side UMI vs genes plots: one colored by mito%, one by cell quality posterior.
    
    Args:
        adata: AnnData object
        cell_barcodes: Dict of cell barcodes by method
        sample_id: Sample ID
        stratify_by: Stratification level
        plot_dir: Output directory for plots
        qc_cell_lists_path: Path to pre-calculated QC cell lists TSV with GMM posteriors
    """
    
    # Plot colored by mitochondrial %
    print("\nCreating UMI vs genes scatter plots colored by % mitochondrial...")
    _plot_umi_scatter_generic(
        adata=adata,
        cell_barcodes=cell_barcodes,
        sample_id=sample_id,
        stratify_by=stratify_by,
        plot_dir=plot_dir,
        x_values_func=lambda data, method: data.obs['total_counts'].values,
        x_label='Total UMI Counts',
        y_values_func=lambda data, method: data.obs['n_genes_by_counts'].values,
        y_label='Number of Genes Detected',
        color_values_func=lambda data, method: data.obs['pct_counts_mt'].values,
        color_label='% Mitochondrial',
        color_map='plasma',
        color_vmin=0,
        color_vmax=None,  # Will use actual max
        metric_name='umi_vs_genes_scatter_mito',
        plot_title_suffix=''
    )
    
    print("UMI vs genes scatter plots completed")


def plot_umi_vs_guides_scatter(adata, cell_barcodes, sample_id, stratify_by, plot_dir, source, processing, guide_cutoffs, pool=None):
    """Create scatter plots of UMI counts vs number of guides at different cutoffs, colored by pct mitochondrial."""
    print("\nCreating UMI vs guides scatter plots colored by % mitochondrial...")
    print(f"  Creating scatter plots for {len(cell_barcodes)} methods with {len(guide_cutoffs)} cutoffs")
    
    # Check if guide data exists
    if 'guide_counts' not in adata.obsm:
        print("  No guide data found in adata.obsm, skipping guide scatter plots")
        return
    
    # Create plots for each cutoff
    for cutoff_spec in guide_cutoffs:
        # Get file suffix and label
        cutoff_suffix = get_cutoff_suffix(cutoff_spec)
        cutoff_suffix_display, cutoff_label = get_cutoff_labels(cutoff_spec)
        print(f"\n  Processing cutoff {cutoff_spec} ({cutoff_label})...")
        
        def calculate_guides_at_cutoff(data, method_name):
            """Calculate number of guides detected at the specified cutoff."""
            guide_counts = data.obsm['guide_counts']
            binary_guides = apply_guide_threshold(guide_counts, cutoff_spec, data.obs, stratify_by, method_name)
            guides_per_cell = np.array(binary_guides.sum(axis=1)).flatten()
            return guides_per_cell
        
        _plot_umi_scatter_generic(
            adata=adata,
            cell_barcodes=cell_barcodes,
            sample_id=sample_id,
            stratify_by=stratify_by,
            plot_dir=plot_dir,
            x_values_func=lambda data, method: data.obs['total_counts'].values,
            x_label='Total UMI Counts',
            y_values_func=calculate_guides_at_cutoff,
            y_label=f'Number of Guides Detected ({cutoff_label})',
            color_values_func=lambda data, method: data.obs['pct_counts_mt'].values,
            color_label='% Mitochondrial',
            color_map='plasma',
            color_vmin=0,
            color_vmax=None,
            metric_name=f'umi_vs_guides_scatter_{cutoff_suffix}',
            plot_title_suffix=cutoff_label
        )
    
    print("\nUMI vs guides scatter plots completed")


def main():
    parser = argparse.ArgumentParser(description='Calculate QC metrics after cell calling by biological sample or well')
    parser.add_argument('--h5ad', required=True, help='Annotated h5ad file')
    parser.add_argument('--cell-calling-dir', required=True, help='Directory with cell calling results')
    parser.add_argument('--read-stats', required=True, help='GEX read statistics TSV file')
    parser.add_argument('--guide-read-stats', required=True, help='Guide read statistics TSV file')
    parser.add_argument('--sample-id', required=True, help='Full sample ID (format: pool:sample)')
    parser.add_argument('--output', required=True, help='Output TSV file (will be used as base path for category-specific files)')
    parser.add_argument('--config', required=True, help='Config YAML file')
    parser.add_argument('--stratify-by', choices=['sample', 'biological_sample', 'well'], default='biological_sample',
                        help='Stratify metrics by sample, biological sample, or well (default: biological_sample)')
    parser.add_argument('--plot-dir', help='Directory to save plots')
    parser.add_argument('--pool', help='Pool name for plot naming')
    parser.add_argument('--source', required=True, choices=['main', 'undetermined', 'all'],
                        help='Data source (main, undetermined, or all)')
    parser.add_argument('--processing', required=True, choices=['raw', 'trimmed', 'recovered', 'merged'],
                        help='Processing state (raw, trimmed, recovered, or merged)')
    parser.add_argument('--qc-cell-lists', help='Pre-calculated QC cell lists TSV file (for plotting GMM posteriors)')
    parser.add_argument('--gmm-thresholds', required=True, help='Pre-calculated GMM thresholds per cell TSV file')
    parser.add_argument('--threads', type=int, required=True, help='Number of threads/cores to use')
    parser.add_argument('--per-cell-plot-method', required=True, help='Cell calling method for per-cell distribution plots (use "all" for all methods)')
    
    args = parser.parse_args()
    
    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)
    
    print(f"Loading h5ad file: {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    print(f"Loaded {adata.n_obs} barcodes x {adata.n_vars} genes")
    
    # Ensure we have the necessary QC metrics
    required_columns = ['pct_counts_mt', 'total_counts', 'n_genes_by_counts']
    missing_columns = [col for col in required_columns if col not in adata.obs.columns]
    if missing_columns:
        raise ValueError(f"Required QC metrics missing from h5ad file: {missing_columns}. Run calculate_qc_metrics.py first.")
    
    # Load read statistics
    print(f"Loading read statistics: {args.read_stats}")
    read_stats_df = pd.read_csv(args.read_stats, sep='\t')
    
    if len(read_stats_df) == 0:
        raise ValueError(f"No data found in {args.read_stats}")
    
    # Load guide read statistics
    print(f"Loading guide read statistics: {args.guide_read_stats}")
    guide_stats_df = pd.read_csv(args.guide_read_stats, sep='\t')
    
    if len(guide_stats_df) == 0:
        raise ValueError(f"No data found in {args.guide_read_stats}")
    
    # Rename guide columns to avoid conflicts using new naming convention
    guide_stats_df = guide_stats_df.rename(columns={col: f'guide_{col}' for col in guide_stats_df.columns if col != 'Sample ID'})
    
    total_reads = int(read_stats_df['total_reads'].iloc[0])
    print(f"Total reads: {total_reads:,}")
    
    # Load guide reference
    if 'qc_analysis' not in config:
        raise ValueError("'qc_analysis' section missing from config")
    if 'guide_reference' not in config['input_paths']:
        raise ValueError("'guide_reference' not specified in input_paths config")
    guide_ref_path = config['input_paths']['guide_reference']
    
    if 'guide_counts' in adata.obsm:
        guide_to_genes = load_guide_reference(guide_ref_path)
    else:
        raise ValueError("No guide counts found in adata.obsm - this is required for guide analysis")
    
    # Get guide cutoffs from config
    if 'guide_cutoffs' not in config['qc_analysis']:
        raise ValueError("'guide_cutoffs' not specified in qc_analysis config")
    guide_cutoffs = config['qc_analysis']['guide_cutoffs']
    print(f"Using guide cutoffs: {guide_cutoffs}")
    
    # Extract posterior levels needed for GMM analysis
    posterior_levels = [50]  # Always include 50% as default
    for cutoff_spec in guide_cutoffs:
        if isinstance(cutoff_spec, str):
            if cutoff_spec.startswith('gmm_') or cutoff_spec.startswith('posterior_'):
                if cutoff_spec.startswith('gmm_'):
                    prob_str = cutoff_spec[4:]
                else:
                    prob_str = cutoff_spec[10:]
                try:
                    prob_level = int(prob_str)
                    if prob_level not in posterior_levels:
                        posterior_levels.append(prob_level)
                except ValueError:
                    raise ValueError(f"Invalid posterior probability specification: {cutoff_spec}")
    
    print(f"GMM posterior levels to calculate: {sorted(posterior_levels)}")
    
    # Get guide mixture model parameters from config
    mixture_config = config.get('guide_mixture_model', {})
    min_umi_threshold = mixture_config.get('min_umi_threshold', 2)
    print(f"Using mixture model parameters: min_umi_threshold={min_umi_threshold}")
    
    # Load cell barcodes for each method
    print(f"Loading cell barcodes from: {args.cell_calling_dir}")
    cell_barcodes = load_cell_barcodes(args.cell_calling_dir, args.sample_id)
    
    # Load QC cell lists for compromised cell percentage calculation
    qc_cell_lists_df = None
    default_method = config['cell_calling']['default_method']
    if args.qc_cell_lists:
        print(f"Loading QC cell lists from: {args.qc_cell_lists}")
        qc_cell_lists_df = pd.read_csv(args.qc_cell_lists, sep='\t')
        print(f"Loaded GMM posteriors for {len(qc_cell_lists_df):,} cells")
    
    # Get unique groups to stratify by
    if args.stratify_by == 'sample':
        # For sample level, just create one group with all data
        groups = ['Overall']
        stratify_column = None
        group_label = 'sample'
        print(f"Calculating metrics at sample level for {args.sample_id}")
    elif args.stratify_by == 'well':
        stratify_column = 'well'
        groups = adata.obs[stratify_column].dropna().unique().tolist()
        print(f"Found {len(groups)} wells")
        group_label = 'well'
    else:  # biological_sample
        stratify_column = 'biological_sample'
        groups = adata.obs[stratify_column].dropna().unique().tolist()
        print(f"Found {len(groups)} biological samples")
        group_label = 'biological_sample'
    
    # Setup plotting parameters
    source = args.source
    processing = args.processing
    plot_dir = Path(args.plot_dir) if args.plot_dir else Path(args.output).parent
    
    # Load pre-calculated GMM thresholds
    threshold_output = Path(args.gmm_thresholds)
    print(f"\n=== Loading pre-calculated GMM thresholds from {threshold_output} ===")
    threshold_table = pd.read_csv(threshold_output, sep='\t', index_col=0)
    
    # Add threshold columns to full adata.obs for use in plotting
    for col in threshold_table.columns:
        # Initialize column with NaN for all cells (means "not applicable")
        adata.obs[col] = np.nan
        # Add the threshold values for cells that have them
        common_indices = adata.obs.index.intersection(threshold_table.index)
        adata.obs.loc[common_indices, col] = threshold_table.loc[common_indices, col]
    print(f"Added {len(threshold_table.columns)} threshold columns to adata.obs")
    
    # Now generate plots using the pre-calculated thresholds
    print("\n=== Generating plots and fitting mixture models ===")
    # Create combined read stats dictionary
    combined_read_stats = read_stats_df.iloc[0].to_dict()
    combined_read_stats.update(guide_stats_df.iloc[0].to_dict())
    
    mixture_metrics = plot_per_cell_distributions(adata, cell_barcodes, groups, args.stratify_by, 
                               plot_dir, args.sample_id, source, processing, guide_cutoffs, 
                               pool=args.pool, read_stats=combined_read_stats,
                               plot_method_filter=args.per_cell_plot_method)
    
    # Initialize results list
    all_results = []
    
    # Calculate metrics for each group
    print("\n=== Calculating metrics ===")
    for group in groups:
        print(f"\nProcessing {group_label}: {group}")
        
        # Create mask for this group
        if group == 'Overall':
            group_mask = np.ones(adata.n_obs, dtype=bool)
        else:
            group_mask = adata.obs[stratify_column] == group
        
        print(f"  Barcodes in group: {np.sum(group_mask):,}")
        
        # Calculate metrics for this group
        results = {
            'sample_id': args.sample_id,
        }
        
        # Add stratification-specific columns
        if args.stratify_by == 'sample':
            # Sample level - include read statistics
            for col in read_stats_df.columns:
                if col != 'Sample ID':  # Skip the ID column
                    results[col] = read_stats_df[col].iloc[0]
            
            # Add guide read statistics columns
            for col in guide_stats_df.columns:
                if col != 'Sample ID':  # Skip the ID column
                    results[col] = guide_stats_df[col].iloc[0]
            
            # Calculate guide read percentage
            if 'total_reads' in results and 'guide_total_reads' in results:
                gex_total_reads = results['total_reads']
                guide_total_reads = results['guide_total_reads']
                total_reads_combined = gex_total_reads + guide_total_reads
                if total_reads_combined > 0:
                    results['guide_read_pct'] = (guide_total_reads / total_reads_combined) * 100
                else:
                    results['guide_read_pct'] = 0.0
        elif args.stratify_by == 'well':
            results['well'] = group
            # Add biological sample info for wells
            bio_samples = adata.obs.loc[group_mask, 'biological_sample'].dropna().unique()
            results['biological_sample'] = bio_samples[0] if len(bio_samples) == 1 else ','.join(bio_samples)
        else:  # biological_sample
            results['biological_sample'] = group
            # Add well info for biological samples
            wells = adata.obs.loc[group_mask, 'well'].dropna().unique()
            results['well'] = ','.join(sorted(wells))
        
        # Add all the per-method metrics
        # Get guide total reads if available (only at sample level)
        guide_total_reads = None
        if args.stratify_by == 'sample' and 'guide_total_reads' in guide_stats_df.columns:
            guide_total_reads = guide_stats_df['guide_total_reads'].iloc[0]
        
        method_results = calculate_metrics_for_group(adata, cell_barcodes, group_mask, total_reads, group, guide_to_genes, guide_cutoffs, guide_total_reads, mixture_metrics, qc_cell_lists_df, default_method, args.stratify_by)
        results.update(method_results)
        
        all_results.append(results)
    
    # Create DataFrame
    results_df = pd.DataFrame(all_results)
    
    # Save the full results file directly (no categories)
    results_df.to_csv(args.output, sep='\t', index=False)
    
    print(f"\nQC metrics by {args.stratify_by} saved to: {args.output}")
    
    if args.stratify_by == 'sample':
        print(f"Total rows: 1 (sample-level metrics)")
    else:
        print(f"Total rows: {len(results_df)} ({args.stratify_by}s)")
    
    # Print summary
    print(f"\nSummary of cell counts by {args.stratify_by}:")
    # Just use the first method for summary
    first_method = list(cell_barcodes.keys())[0]
    for idx, row in results_df.iterrows():
        if args.stratify_by == 'sample':
            group_name = args.sample_id
        elif args.stratify_by == 'well':
            group_name = f"Well {row['well']}"
        else:  # biological_sample
            group_name = row['biological_sample']
        n_cells = row[f'n_cells_{first_method}']
        print(f"  {group_name}: {int(n_cells):,} cells ({first_method})")
    
    # Generate additional UMI vs genes scatter plots (both mito% and cell quality posterior if available)
    plot_umi_vs_genes_with_cell_quality(adata, cell_barcodes, args.sample_id, args.stratify_by, plot_dir, 
                                qc_cell_lists_path=args.qc_cell_lists)
    
    # Generate UMI vs guides scatter plots for each cutoff
    if 'guide_counts' in adata.obsm:
        plot_umi_vs_guides_scatter(adata, cell_barcodes, args.sample_id, args.stratify_by,
                                 plot_dir, source, processing, guide_cutoffs, pool=args.pool)
        
        # Generate GEX UMIs vs Guide UMIs scatter plot
        print("\nCreating GEX UMIs vs Guide UMIs scatter plots colored by % mitochondrial...")
        print(f"  Creating scatter plots for {len(cell_barcodes)} methods: {', '.join(sorted(cell_barcodes.keys()))}")
        
        _plot_umi_scatter_generic(
            adata=adata,
            cell_barcodes=cell_barcodes,
            sample_id=args.sample_id,
            stratify_by=args.stratify_by,
            plot_dir=plot_dir,
            x_values_func=lambda data, method: data.obs['total_counts'].values,
            x_label='Total UMI Counts',
            y_values_func=lambda data, method: np.array(data.obsm['guide_counts'].sum(axis=1)).flatten(),
            y_label='Guide UMIs',
            color_values_func=lambda data, method: data.obs['pct_counts_mt'].values,
            color_label='% Mitochondrial',
            color_map='plasma',
            color_vmin=0,
            color_vmax=None,
            metric_name='gex_vs_guide_umis_scatter',
            plot_title_suffix=''
        )
        
        print("GEX vs Guide UMIs scatter plots completed")
        
        # Generate Guide UMIs vs Number of Guides scatter plots for each cutoff
        print("\nCreating Guide UMIs vs Number of Guides scatter plots colored by % mitochondrial...")
        print(f"  Creating scatter plots for {len(cell_barcodes)} methods with {len(guide_cutoffs)} cutoffs: {guide_cutoffs}")
        
        for cutoff_spec in guide_cutoffs:
            # Get file suffix and label
            cutoff_suffix = get_cutoff_suffix(cutoff_spec)
            cutoff_suffix_display, cutoff_label = get_cutoff_labels(cutoff_spec)
            if cutoff_suffix is None:
                continue
            print(f"  Processing cutoff {cutoff_spec} ({cutoff_label})...")
            
            def calculate_guides_at_cutoff(data, method_name):
                """Calculate number of guides detected at the specified cutoff."""
                guide_counts = data.obsm['guide_counts']
                binary_guides = apply_guide_threshold(guide_counts, cutoff_spec, data.obs, args.stratify_by, method_name)
                guides_per_cell = np.array(binary_guides.sum(axis=1)).flatten()
                return guides_per_cell
            
            # Use the now-generic scatter plot function with guide UMIs as x-axis
            _plot_umi_scatter_generic(
                adata=adata,
                cell_barcodes=cell_barcodes,
                sample_id=args.sample_id,
                stratify_by=args.stratify_by,
                plot_dir=plot_dir,
                x_values_func=lambda data, method: np.array(data.obsm['guide_counts'].sum(axis=1)).flatten(),
                x_label='Guide UMIs',
                y_values_func=calculate_guides_at_cutoff,
                y_label=f'Number of Guides Detected ({cutoff_label})',
                color_values_func=lambda data, method: data.obs['pct_counts_mt'].values,
                color_label='% Mitochondrial',
                color_map='plasma',
                color_vmin=0,
                color_vmax=None,
                metric_name=f'guide_umis_vs_num_guides_scatter_{cutoff_suffix}',
                plot_title_suffix=cutoff_label
            )
        
        print("Guide UMIs vs Number of Guides scatter plots completed")
    
    # Clean up memory
    del adata
    gc.collect()


if __name__ == "__main__":
    main()