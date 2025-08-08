#!/usr/bin/env python3
"""
Visualize aggregated QC metrics from consolidated TSV files.

Creates horizontal bar plots for sample/biological sample metrics and
heatmaps for well-level metrics. Uses directory structure instead of JSON metadata.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re
import sys
from multiprocessing import Pool, cpu_count
from functools import partial

# No longer need plot_config since categories come from the input file structure


def create_horizontal_point_plot(df, metric_col, stratify_by, log_scale=False, show_context=True):
    """Create horizontal point plot for a metric."""
    # Get group identifier column
    if stratify_by == 'sample':
        group_col = 'sample_id'
    elif stratify_by == 'biological_sample':
        group_col = 'biological_sample'
    else:
        group_col = stratify_by
    
    # Sort by metric value
    sorted_df = df.sort_values(metric_col, ascending=True)
    
    # Create display labels
    if stratify_by == 'sample':
        # For sample level, just show sample_id as-is (already contains pool:sample)
        display_labels = sorted_df[group_col]
    elif stratify_by == 'biological_sample' and 'sample_id' in sorted_df.columns:
        # For biological sample level, show full sample_id (which includes pool:sample)
        display_labels = sorted_df.apply(lambda row: f"{row['sample_id']}: {row[group_col]}", axis=1)
    else:
        display_labels = sorted_df[group_col]
    
    # Extract method from metric_col using same logic as save_plot_with_structure
    method = None
    parts = metric_col.split('_')
    for i, part in enumerate(parts):
        if part and part[0].isupper() and i > 0:  # Methods start with capital
            method_parts = parts[i:]
            if len(method_parts) > 1 or any(p.isupper() for p in method_parts):
                method = '_'.join(method_parts)
                break
    
    # If this is a cell-based metric (has a method), add cell counts to labels
    if method:
        n_cells_col = f'n_cells_{method}'
        if n_cells_col in sorted_df.columns:
            display_labels = sorted_df.apply(
                lambda row: f"{display_labels.iloc[sorted_df.index.get_loc(row.name)]} (n={int(row[n_cells_col]):,})",
                axis=1
            )
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, max(6, len(sorted_df) * 0.3)))
    
    # Create point plot (scatter)
    positions = np.arange(len(sorted_df))
    ax.scatter(sorted_df[metric_col], positions, s=100, alpha=0.7, edgecolors='black', linewidth=1)
    
    # Note: Pool coloring removed to match JSON version which uses simple scatter plots
    
    # Set labels
    ax.set_yticks(positions)
    ax.set_yticklabels(display_labels)
    
    # Set simple title and labels to match JSON version
    ax.set_xlabel(metric_col)
    
    # Add context to title if available
    title = f'{metric_col} by {stratify_by.replace("_", " ").title()}'
    
    # For cell-based metrics, show which method's cells are used
    if show_context and metric_col.startswith('n_cells_'):
        # Extract method from column name (e.g., n_cells_EmptyDrops_FDR001)
        method = metric_col.replace('n_cells_', '')
        title += f'\nMethod: {method}'
    
    ax.set_title(title)
    
    # Add vertical grid lines for easier reading
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # Log scale if requested
    if log_scale and sorted_df[metric_col].min() > 0:
        ax.set_xscale('log')
    
    # Value labels removed for cleaner visualization
    
    # Set x-axis limits with some padding
    if not sorted_df[metric_col].isna().all():
        data_min = sorted_df[metric_col].min()
        data_max = sorted_df[metric_col].max()
        data_range = data_max - data_min
        
        if log_scale and data_min > 0:
            # For log scale, use multiplicative padding
            ax.set_xlim(data_min * 0.8, data_max * 1.2)
        else:
            # For linear scale, add 10% padding on each side
            # Make sure padding is enough to include text labels
            padding = max(data_range * 0.1, data_max * 0.05) if data_range > 0 else data_max * 0.05
            ax.set_xlim(data_min - padding, data_max + padding)
    
    # Adjust layout
    plt.tight_layout()
    return fig


def create_scatter_plot(df, x_col, y_col, stratify_by, log_scale=False):
    """Create scatter plot for two metrics."""
    # Get group identifier column
    if stratify_by == 'sample':
        group_col = 'sample_id'
    elif stratify_by == 'biological_sample':
        group_col = 'biological_sample'
    else:
        group_col = stratify_by
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create scatter plot
    ax.scatter(df[x_col], df[y_col], s=100, alpha=0.7, edgecolors='black', linewidth=1)
    
    # Add labels for each point
    for idx, row in df.iterrows():
        label = row[group_col]
        if stratify_by == 'sample':
            # For sample level, sample_id already contains pool:sample
            label = row[group_col]
        elif stratify_by == 'biological_sample' and 'sample_id' in df.columns:
            # For biological sample, show full sample_id to avoid confusion
            label = f"{row['sample_id']}: {label}"
        ax.annotate(label, (row[x_col], row[y_col]), xytext=(5, 5), 
                   textcoords='offset points', fontsize=8, alpha=0.7)
    
    # Set labels
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(f'{y_col} vs {x_col}')
    
    # Add grid
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # Log scale if requested
    if log_scale:
        if df[x_col].min() > 0:
            ax.set_xscale('log')
        if df[y_col].min() > 0:
            ax.set_yscale('log')
    
    plt.tight_layout()
    return fig


def create_well_heatmap(df, metric_col):
    """Create heatmap for well-level metrics."""
    # Create empty 96-well plate (8 rows x 12 columns)
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    cols = range(1, 13)
    
    # Initialize plate matrix
    plate_data = pd.DataFrame(np.nan, index=rows, columns=cols)
    
    # Fill in values from data
    for idx, row in df.iterrows():
        well = row['well']
        if pd.notna(well) and len(well) >= 2:
            # Parse well position (e.g., 'A1', 'B12')
            well_row = well[0]
            well_col = int(well[1:])
            
            if well_row in rows and 1 <= well_col <= 12:
                plate_data.loc[well_row, well_col] = row[metric_col]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create mask for NaN values
    mask = plate_data.isna()
    
    # Plot heatmap
    sns.heatmap(plate_data, mask=mask, cmap='YlOrRd', center=None,
                square=True, linewidths=0.5, cbar_kws={"shrink": 0.8},
                ax=ax, annot=True, fmt='.1f', annot_kws={'size': 8})
    
    # Set title and labels
    ax.set_title(f'{metric_col} - 96-Well Plate Heatmap')
    ax.set_xlabel('Column')
    ax.set_ylabel('Row')
    
    return fig


def save_plot_with_structure(fig, metric_col, stratify_by, source, processing, 
                           output_base_dir, scale='linear', category=None):
    """Save plot using clean directory structure."""
    # Extract method from metric_col if present
    # Pattern: metric_name_Method_Name
    method = None
    metric_name = metric_col
    
    # Look for method patterns in the metric name
    parts = metric_col.split('_')
    for i, part in enumerate(parts):
        if part and part[0].isupper() and i > 0:  # Methods start with capital
            # Found potential method start
            method_parts = parts[i:]
            # Check if this looks like a method (has underscores or all caps)
            if len(method_parts) > 1 or any(p.isupper() for p in method_parts):
                method = '_'.join(method_parts)
                metric_name = '_'.join(parts[:i])
                break
    
    # Determine which consolidated category based on whether metric has a method
    if method:
        # Cell-based metrics with methods
        plot_dir = output_base_dir / 'consolidated_cell_based' / f'{source}_{processing}' / stratify_by / metric_name / method / scale
    else:
        # General metrics without methods
        plot_dir = output_base_dir / 'consolidated_general' / f'{source}_{processing}' / stratify_by / metric_name / scale
    
    plot_dir.mkdir(parents=True, exist_ok=True)
    
    # Standard filename
    filename = 'plot.png'
    filepath = plot_dir / filename
    
    # Save plot
    fig.savefig(filepath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    
    return filepath


def process_metric(args_tuple):
    """Worker function to process a single metric."""
    metric_col, df, stratify_by, source, processing, output_base_dir, category = args_tuple
    
    results = []
    
    # Define metrics that should be plotted as scatter plots vs mean reads per cell
    scatter_plot_metrics = [
        'mean_umis_per_cell',
        'median_umis_per_cell',
        'mean_genes_per_cell',
        'median_genes_per_cell',
        # Background quality metrics (mirror cell metrics)
        'mean_umis_per_noncell',
        'median_umis_per_noncell',
        'mean_genes_per_noncell',
        'median_genes_per_noncell',
        # Saturation and guide metrics
        'reads_per_umi_in_cells',
        'total_umis_in_cells',
        'guide_umis_per_cell',
        'guides_per_cell_cutoff',
        'fraction_cells_with_guides_cutoff',
        'umis_per_guide_per_cell_cutoff',
        'total_guides_detected_cutoff',
        'total_guide_detections_cutoff',
        'total_targeted_genes_cutoff',
    ]
    
    # Define metrics that should be plotted as scatter plots vs total reads (group-level)
    scatter_vs_reads_metrics = [
        'total_umis',
        'reads_per_umi',
        'pcr_duplication_rate_pct',
        'guide_total_umis',
        'guide_pcr_duplication_rate_pct',
        'guide_reads_per_umi',
    ]
    
    # Define metrics to skip entirely (not plot at all)
    skip_metrics = [
        'reads_with_correctable_barcodes',
        'mapped_reads',
        'guide_expected_cells',
        'guide_mapped_reads',
        'guide_min_umi_threshold',
        'guide_reads_with_correctable_barcodes',
        'guide_sequencing_saturation_pct',
        'pct_mapped_reads_with_correctable_barcodes',
        'min_umi_threshold',
    ]
    
    # Check if this metric should be skipped
    if metric_col in skip_metrics:
        return results
    
    # Check if this metric should be a scatter plot
    is_scatter_metric = any(metric_col.startswith(scatter_metric) for scatter_metric in scatter_plot_metrics)
    is_scatter_vs_reads = metric_col in scatter_vs_reads_metrics
    
    if stratify_by == 'well':
        # Create heatmap
        fig = create_well_heatmap(df, metric_col)
        filepath = save_plot_with_structure(fig, metric_col, stratify_by, 
                                          source, processing, output_base_dir, 'linear', category)
        results.append((metric_col, 'heatmap', filepath))
    elif is_scatter_metric and 'total_reads' in df.columns:
        # Create scatter plot with mean reads per cell on x-axis
        # First need to find corresponding n_cells column
        # Extract method from metric name if present
        method_match = None
        if '_' in metric_col:
            parts = metric_col.split('_')
            for i, part in enumerate(parts):
                if part and part[0].isupper() and i > 0:
                    method_match = '_'.join(parts[i:])
                    break
        
        # Find corresponding n_cells column
        n_cells_col = None
        if method_match:
            n_cells_col = f'n_cells_{method_match}'
            if n_cells_col not in df.columns:
                # Try alternative naming
                n_cells_col = f'n_cells_called_{method_match}'
        
        # If we found the n_cells column, calculate mean reads per cell
        if n_cells_col and n_cells_col in df.columns:
            # Create a copy of df to add the calculated column
            plot_df = df.copy()
            plot_df['mean_reads_per_cell'] = plot_df['total_reads'] / plot_df[n_cells_col]
            plot_df['mean_reads_per_cell'] = plot_df['mean_reads_per_cell'].replace([np.inf, -np.inf], np.nan)
            
            # Create scatter plots (both linear and log scale)
            for scale, log_scale in [('linear', False), ('log', True)]:
                fig = create_scatter_plot(plot_df, 'mean_reads_per_cell', metric_col, stratify_by, log_scale)
                filepath = save_plot_with_structure(fig, metric_col, stratify_by,
                                                  source, processing, output_base_dir, scale, category)
                results.append((metric_col, scale, filepath))
        # No else - if we can't find n_cells, we simply skip this metric
    elif is_scatter_vs_reads:
        # Determine which Total reads column to use
        if metric_col.startswith('guide_') and 'guide_total_reads' in df.columns:
            x_col = 'guide_total_reads'
        elif 'total_reads' in df.columns:
            x_col = 'total_reads'
        else:
            # Skip if we can't find the appropriate Total reads column
            return results
            
        # Create scatter plot with appropriate Total reads on x-axis
        for scale, log_scale in [('linear', False), ('log', True)]:
            fig = create_scatter_plot(df, x_col, metric_col, stratify_by, log_scale)
            filepath = save_plot_with_structure(fig, metric_col, stratify_by,
                                              source, processing, output_base_dir, scale, category)
            results.append((metric_col, scale, filepath))
    elif not is_scatter_vs_reads:  # Only create point plots if NOT a scatter plot metric
        # Create horizontal point plots (both linear and log scale)
        for scale, log_scale in [('linear', False), ('log', True)]:
            # Skip log scale for certain metrics
            skip_patterns = [
                'fraction_', 'pct_', 'percent', 'expected_cells', 
                'ratio', 'coefficient', 'score', 'index',
                'n_cells_', 'total_cells', 'cell_count'
            ]
            if log_scale and any(skip in metric_col.lower() for skip in skip_patterns):
                continue
            
            # Skip if all values are zero, negative, or missing for log scale
            valid_values = df[metric_col].dropna()
            if log_scale and (len(valid_values) == 0 or (valid_values <= 0).all()):
                continue
            
            fig = create_horizontal_point_plot(df, metric_col, stratify_by, log_scale)
            filepath = save_plot_with_structure(fig, metric_col, stratify_by,
                                              source, processing, output_base_dir, scale, category)
            results.append((metric_col, scale, filepath))
    
    return results


def main():
    parser = argparse.ArgumentParser(description='Visualize aggregated QC metrics')
    parser.add_argument('--input', required=True, help='Input TSV file with aggregated metrics or base path for category files')
    parser.add_argument('--stratify_by', required=True, 
                       choices=['sample', 'biological_sample', 'well', 'pool'],
                       help='Level to stratify plots by')
    parser.add_argument('--output_dir', required=True, help='Base output directory for plots')
    parser.add_argument('--source', required=True, 
                       choices=['main', 'undetermined', 'all'],
                       help='Data source type')
    parser.add_argument('--processing', required=True,
                       choices=['raw', 'recovered', 'merged'],
                       help='Processing type')
    parser.add_argument('--metrics', nargs='+', 
                       help='Specific metrics to plot (default: all)')
    parser.add_argument('--threads', type=int, default=1,
                       help='Number of threads to use for parallel processing')
    
    args = parser.parse_args()
    
    # Load single consolidated file
    print(f"Loading data from {args.input}")
    df = pd.read_csv(args.input, sep='\t')
    print(f"Loaded {len(df)} rows with {len(df.columns)} columns")
    
    # Get metric columns (exclude identifier columns)
    id_cols = ['pool', 'sample_id', 'biological_sample', 'well', 'row', 'column']
    metric_cols = [col for col in df.columns if col not in id_cols]
    
    # Filter to specific metrics if requested
    if args.metrics:
        metric_cols = [col for col in metric_cols if any(m in col for m in args.metrics)]
    
    print(f"\\nGenerating plots for {len(metric_cols)} metrics")
    
    output_base_dir = Path(args.output_dir)
    
    # Prepare arguments for parallel processing
    plot_args = [
        (metric_col, df, args.stratify_by, args.source, args.processing, output_base_dir, None)
        for metric_col in metric_cols
    ]
    
    # Use multiprocessing to generate plots in parallel
    n_cores = min(args.threads, cpu_count())
    print(f"Using {n_cores} cores for parallel processing")
    
    with Pool(n_cores) as pool:
        # Process metrics in parallel with progress tracking
        results = []
        for i, result in enumerate(pool.imap_unordered(process_metric, plot_args), 1):
            results.extend(result)
            print(f"\\r[{i}/{len(metric_cols)}] Completed {i} metrics...", end='', flush=True)
    
    print(f"\\n\\nGenerated {len(results)} plots")
    print(f"Plots saved to:")
    print(f"  - {output_base_dir / 'consolidated_general' / f'{args.source}_{args.processing}'}")
    print(f"  - {output_base_dir / 'consolidated_cell_based' / f'{args.source}_{args.processing}'}")


if __name__ == "__main__":
    main()