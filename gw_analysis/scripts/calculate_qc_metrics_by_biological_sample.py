#!/usr/bin/env python3
"""
Calculate comprehensive QC metrics at three stratification levels.

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
from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import yaml
import scipy.sparse
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
# Import shared guide utility function
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.pipeline_utils import calculate_guides_per_cell, calculate_guide_metrics_for_cells


def plot_violin_with_mean(data, x, y, order, ax, log_scale=False):
    """Helper function to create violin plot with mean markers.
    
    Args:
        data: DataFrame with data to plot
        x: Column name for x-axis grouping
        y: Column name for y-axis values
        order: Order of groups on x-axis
        ax: Matplotlib axis to plot on
        log_scale: Whether to use log scale
        
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
    
    return nonzero_fractions





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


def load_guide_reference(guide_ref_path):
    """Load guide reference and create mapping from guide sequence to list of target genes."""
    print(f"Loading guide reference from: {guide_ref_path}")
    
    # Read the TSV file
    guide_df = pd.read_csv(guide_ref_path, sep='\t')
    
    # Create mapping from guide sequence to list of genes
    guide_to_genes = {}
    
    for _, row in guide_df.iterrows():
        guide_id = row['ID']  # Use ID instead of sequence
        gene_str = row['gene']
        
        # Skip if gene is NaN or empty
        if pd.isna(gene_str) or gene_str == '':
            continue
            
        # Split multi-gene targets by "_AND_"
        if '_AND_' in gene_str:
            genes = gene_str.split('_AND_')
        else:
            genes = [gene_str]
        
        guide_to_genes[guide_id] = genes
    
    print(f"Loaded {len(guide_to_genes)} guide IDs mapping to {len(set(sum(guide_to_genes.values(), [])))} unique genes")
    
    return guide_to_genes


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


def calculate_guide_metrics_for_barcodes(adata, barcode_set, guide_to_genes=None, cutoffs=None):
    """Calculate guide metrics for a given set of barcodes at multiple cutoffs."""
    if cutoffs is None:
        raise ValueError("Guide cutoffs must be provided")
    
    # Convert barcode_set to list for the unified function
    if barcode_set:
        barcode_list = list(barcode_set)
    else:
        # Empty set - no cells
        barcode_list = []
    
    if len(barcode_list) == 0:
        raise ValueError("No barcodes provided to calculate_guide_metrics_for_barcodes")
    
    # Extract guide data from adata
    # Create a temporary AnnData with just guide counts
    guide_adata = sc.AnnData(X=adata.obsm['guide_counts'])
    guide_adata.obs_names = adata.obs_names
    if 'guide_names' in adata.uns:
        guide_adata.uns['guide_names'] = adata.uns['guide_names']
    
    # Use the unified function
    unified_metrics = calculate_guide_metrics_for_cells(
        guide_adata=guide_adata,
        cell_barcodes=barcode_list,
        guide_cutoffs=cutoffs,
        guide_to_genes=guide_to_genes
    )
    
    # Map the unified metrics to the expected format
    # The unified function returns more detailed metrics, extract what we need
    metrics = {
        'guide_umis_per_cell': unified_metrics['guide_umis_per_cell']
    }
    
    for cutoff in cutoffs:
        metrics[f'guides_per_cell_cutoff{cutoff}'] = unified_metrics[f'guides_per_cell_cutoff{cutoff}']
        metrics[f'fraction_cells_with_guides_cutoff{cutoff}'] = unified_metrics[f'fraction_cells_with_guides_cutoff{cutoff}']
        metrics[f'umis_per_guide_per_cell_cutoff{cutoff}'] = unified_metrics[f'umis_per_guide_per_cell_cutoff{cutoff}']
        metrics[f'total_guides_detected_cutoff{cutoff}'] = unified_metrics[f'total_guides_detected_cutoff{cutoff}']
        metrics[f'total_guide_detections_cutoff{cutoff}'] = unified_metrics[f'total_guide_detections_cutoff{cutoff}']
        metrics[f'total_targeted_genes_cutoff{cutoff}'] = unified_metrics[f'total_targeted_genes_cutoff{cutoff}']
    
    return metrics


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


def calculate_metrics_for_group(adata, cell_barcodes, group_mask, total_reads, group_name, guide_to_genes=None, guide_cutoffs=None, guide_total_reads=None):
    """Calculate all metrics for a specific group (sample/biological sample/well)."""
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
        guide_metrics = calculate_guide_metrics_for_barcodes(adata[group_mask], group_barcodes, guide_to_genes, guide_cutoffs)
        
        # Calculate reads per UMI for cells
        # Note: For biological samples, we can't calculate a meaningful reads_per_umi
        # because we don't know what fraction of total reads came from each bio sample
        # Only meaningful for Overall
        if group_name == 'Overall':
            if cell_metrics['total_umis'] == 0:
                raise ValueError("No UMIs found in cells for reads per UMI calculation")
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
        
        for cutoff in guide_cutoffs:
            results[f'guides_per_cell_cutoff{cutoff}_{method}'] = guide_metrics[f'guides_per_cell_cutoff{cutoff}']
            results[f'fraction_cells_with_guides_cutoff{cutoff}_{method}'] = guide_metrics[f'fraction_cells_with_guides_cutoff{cutoff}']
            results[f'umis_per_guide_per_cell_cutoff{cutoff}_{method}'] = guide_metrics[f'umis_per_guide_per_cell_cutoff{cutoff}']
            
            # Guide diversity
            results[f'total_guides_detected_cutoff{cutoff}_{method}'] = guide_metrics[f'total_guides_detected_cutoff{cutoff}']
            results[f'total_guide_detections_cutoff{cutoff}_{method}'] = guide_metrics[f'total_guide_detections_cutoff{cutoff}']
            results[f'total_targeted_genes_cutoff{cutoff}_{method}'] = guide_metrics[f'total_targeted_genes_cutoff{cutoff}']
    
    return results


def plot_per_cell_distributions(adata, cell_barcodes, groups, stratify_by, output_dir, sample_id, source, processing, guide_cutoffs=None, pool=None, read_stats=None):
    """Generate violin plots for per-cell metrics as individual PNG files."""
    # Skip plotting for wells (too many groups)
    if stratify_by == 'well':
        print("Skipping per-cell distribution plots for well-level analysis (too many groups)")
        return
    
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
    if has_guides and guide_cutoffs is None:
        raise ValueError("Guide cutoffs must be provided when guide data exists")
    
    # Process each cell calling method
    for method, method_barcodes in cell_barcodes.items():
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
                plot_violin_with_mean(plot_df, 'group', 'value', sorted_groups, ax, log_scale=(scale_type == 'log'))
            
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
                plt.savefig(filepath, dpi=300, bbox_inches='tight')
                plt.close()
                
                
                print(f"    Saved {filename} ({scale_type} scale)")
        
        # Plot guide metrics if available
        if has_guides:
            for cutoff in guide_cutoffs:
                # Calculate guides per cell for this cutoff
                guides_per_cell_data = []
                guides_per_cell_labels = []
                guides_group_means = {}
                
                for i, (data, label) in enumerate(zip(group_data, group_labels)):
                    # Get guide counts for these cells
                    guide_counts = data.obsm['guide_counts']
                    # Convert to binary (1 if count >= cutoff, 0 otherwise)
                    binary_guides = (guide_counts >= cutoff).astype(int)
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
                    plot_violin_with_mean(plot_df, 'group', 'value', sorted_groups, ax, log_scale=(scale_type == 'log'))
                
                    # Customize plot
                    ax.set_xlabel(stratify_by.replace('_', ' ').title())
                    ax.set_ylabel(f'Guides per cell (cutoff={cutoff})' + (' (log scale)' if scale_type == 'log' else ''))
                    
                    # Add guide read count and cell count to title
                    n_cells = len(method_barcodes)
                    if read_stats and 'guide_total_reads' in read_stats:
                        guide_total_reads = read_stats['guide_total_reads']
                        guide_reads_per_cell = guide_total_reads // n_cells if n_cells > 0 else 0
                        title = f'Guides per Cell Distribution (cutoff={cutoff}) by {stratify_by.replace("_", " ").title()} ({method}) - {scale_type.capitalize()} Scale\nCells: {n_cells:,} | Guide reads/cell: {guide_reads_per_cell:,} | Percentiles: 1%-99%'
                    else:
                        title = f'Guides per Cell Distribution (cutoff={cutoff}) by {stratify_by.replace("_", " ").title()} ({method}) - {scale_type.capitalize()} Scale\nCells: {n_cells:,} | Percentiles: 1%-99%'
                    ax.set_title(title)
                    
                    # Rotate x labels if many groups
                    if len(sorted_groups) > 5:
                        plt.xticks(rotation=45, ha='right')
                    
                    # Adjust layout
                    plt.tight_layout()
                    
                    # Save plot with clean directory structure
                    # per_cell/{source}_{processing}/{pool}/{sample}/guides_per_cell_cutoff{cutoff}/{method}/{stratify_by}/plot_{scale_type}.png
                    metric_name = f"guides_per_cell_cutoff{cutoff}"
                    metric_dir = plot_dir / metric_name / method / stratify_by / scale_type
                    metric_dir.mkdir(parents=True, exist_ok=True)
                    
                    filename = "plot.png"
                    filepath = metric_dir / filename
                    plt.savefig(filepath, dpi=300, bbox_inches='tight')
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
                plot_violin_with_mean(plot_df, 'group', 'value', sorted_groups, ax, log_scale=(scale_type == 'log'))
            
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
                plt.savefig(filepath, dpi=300, bbox_inches='tight')
                plt.close()
                
                
                print(f"    Saved {filename} ({scale_type} scale)")
            
            # Plot UMIs per guide per cell for each cutoff
            for cutoff in guide_cutoffs:
                umis_per_guide_data = []
                umis_per_guide_labels = []
                umis_per_guide_means = {}
                
                for i, (data, label) in enumerate(zip(group_data, group_labels)):
                    # Get guide counts for these cells
                    guide_counts = data.obsm['guide_counts']
                    # Convert to binary (1 if count >= cutoff, 0 otherwise)
                    binary_guides = (guide_counts >= cutoff).astype(int)
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
                        plot_violin_with_mean(plot_df, 'group', 'value', sorted_groups_umis, ax, log_scale=(scale_type == 'log'))
                    
                        # Customize plot
                        ax.set_xlabel(stratify_by.replace('_', ' ').title())
                        ax.set_ylabel(f'UMIs per guide per cell (cutoff={cutoff})' + (' (log scale)' if scale_type == 'log' else ''))
                        
                        # Add guide read count and cell count to title
                        n_cells = len(method_barcodes)
                        if read_stats and 'guide_total_reads' in read_stats:
                            guide_total_reads = read_stats['guide_total_reads']
                            guide_reads_per_cell = guide_total_reads // n_cells if n_cells > 0 else 0
                            title = f'UMIs per Guide per Cell Distribution (cutoff={cutoff}) by {stratify_by.replace("_", " ").title()} ({method}) - {scale_type.capitalize()} Scale\nCells: {n_cells:,} | Guide reads/cell: {guide_reads_per_cell:,} | Percentiles: 1%-99%'
                        else:
                            title = f'UMIs per Guide per Cell Distribution (cutoff={cutoff}) by {stratify_by.replace("_", " ").title()} ({method}) - {scale_type.capitalize()} Scale\nCells: {n_cells:,} | Percentiles: 1%-99%'
                        ax.set_title(title)
                        
                        # Rotate x labels if many groups
                        if len(sorted_groups_umis) > 5:
                            plt.xticks(rotation=45, ha='right')
                        
                        # Adjust layout
                        plt.tight_layout()
                        
                        # Save plot with category subdirectory
                        # Save plot with clean directory structure
                        metric_name = f"umis_per_guide_per_cell_cutoff{cutoff}"
                        metric_dir = plot_dir / metric_name / method / stratify_by / scale_type
                        metric_dir.mkdir(parents=True, exist_ok=True)
                        
                        filename = "plot.png"
                        filepath = metric_dir / filename
                        plt.savefig(filepath, dpi=300, bbox_inches='tight')
                        plt.close()
                        
                        
                        print(f"    Saved {filename} ({scale_type} scale)")
    
    print("Per-cell distribution plots completed")


def plot_umi_vs_genes_scatter(adata, cell_barcodes, sample_id, stratify_by, plot_dir, source, processing, pool=None):
    """Create scatter plots of UMI counts vs number of genes, colored by pct mitochondrial."""
    from scipy.stats import pearsonr
    import matplotlib.cm as cm
    
    print("\nCreating UMI vs genes scatter plots colored by % mitochondrial...")
    
    # Ensure plot directory exists
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(exist_ok=True, parents=True)
    
    # Use all available methods from cell_barcodes
    methods = sorted(cell_barcodes.keys())
    print(f"  Creating scatter plots for {len(methods)} methods: {', '.join(methods)}")
    
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
            
            # Get the data
            total_counts = data.obs['total_counts'].values
            n_genes = data.obs['n_genes_by_counts'].values
            pct_mito = data.obs['pct_counts_mt'].values
            
            # Sample if too many points
            max_points = 10000
            if len(total_counts) > max_points:
                # Random sample while preserving distribution
                sample_idx = np.random.choice(len(total_counts), size=max_points, replace=False)
                total_counts_plot = total_counts[sample_idx]
                n_genes_plot = n_genes[sample_idx]
                pct_mito_plot = pct_mito[sample_idx]
                plot_label = f"(showing {max_points:,} of {len(total_counts):,} cells)"
            else:
                total_counts_plot = total_counts
                n_genes_plot = n_genes
                pct_mito_plot = pct_mito
                plot_label = f"({len(total_counts):,} cells)"
            
            # Create scatter plot
            scatter = ax.scatter(
                total_counts_plot,
                n_genes_plot,
                c=pct_mito_plot,
                cmap='plasma',
                s=1,
                alpha=0.5,
                vmin=0,
                vmax=np.percentile(pct_mito, 95)  # Cap colormap at 95th percentile
            )
            
            # Set log scale
            ax.set_xscale('log')
            ax.set_yscale('log')
            
            # Labels and title
            ax.set_xlabel('Total UMI Counts')
            ax.set_ylabel('Number of Genes Detected')
            ax.set_title(f'{label}\n{plot_label}')
            
            # Add colorbar
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('% Mitochondrial', rotation=270, labelpad=15)
            
            # Calculate and add correlation
            if len(total_counts) > 1:
                # Calculate on log-transformed data
                log_counts = np.log10(total_counts + 1)
                log_genes = np.log10(n_genes + 1)
                corr, _ = pearsonr(log_counts, log_genes)
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
        fig.suptitle(f'UMI Counts vs Number of Genes (colored by % Mitochondrial)\n{method} - {sample_id}', 
                     fontsize=14)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save with clean directory structure (no separate correlation subdirectory)
        metric_name = "umi_vs_genes_scatter"
        metric_dir = plot_dir / metric_name / method / stratify_by / "linear"
        metric_dir.mkdir(parents=True, exist_ok=True)
        
        filename = "plot.png"
        filepath = metric_dir / filename
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        
        
        print(f"    Saved {filename}")
    
    print("UMI vs genes scatter plots completed")


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
    if 'guide_reference' not in config['qc_analysis']:
        raise ValueError("'guide_reference' not specified in qc_analysis config")
    guide_ref_path = config['qc_analysis']['guide_reference']
    if 'guide_counts' in adata.obsm:
        guide_to_genes = load_guide_reference(guide_ref_path)
    else:
        raise ValueError("No guide counts found in adata.obsm - this is required for guide analysis")
    
    # Get guide cutoffs from config
    if 'guide_cutoffs' not in config['qc_analysis']:
        raise ValueError("'guide_cutoffs' not specified in qc_analysis config")
    guide_cutoffs = config['qc_analysis']['guide_cutoffs']
    print(f"Using guide cutoffs: {guide_cutoffs}")
    
    # Load cell barcodes for each method
    print(f"Loading cell barcodes from: {args.cell_calling_dir}")
    cell_barcodes = load_cell_barcodes(args.cell_calling_dir, args.sample_id)
    
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
    
    # Initialize results list
    all_results = []
    
    # Calculate metrics for each group
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
        
        method_results = calculate_metrics_for_group(adata, cell_barcodes, group_mask, total_reads, group, guide_to_genes, guide_cutoffs, guide_total_reads)
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
    
    # Use provided parameters
    source = args.source
    processing = args.processing
    
    # Generate per-cell distribution plots
    plot_dir = Path(args.plot_dir) if args.plot_dir else Path(args.output).parent
    
    # Create combined read stats dictionary
    combined_read_stats = read_stats_df.iloc[0].to_dict()
    combined_read_stats.update(guide_stats_df.iloc[0].to_dict())
    
    plot_per_cell_distributions(adata, cell_barcodes, groups, args.stratify_by, 
                               plot_dir, args.sample_id, source, processing, guide_cutoffs, 
                               pool=args.pool, read_stats=combined_read_stats)
    
    # Generate UMI vs genes scatter plots
    plot_umi_vs_genes_scatter(adata, cell_barcodes, args.sample_id, args.stratify_by,
                             plot_dir, source, processing, pool=args.pool)
    
    # Clean up memory
    del adata
    import gc
    gc.collect()


if __name__ == "__main__":
    main()