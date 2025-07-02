#!/usr/bin/env python3
"""
Cell Calling Plots Generator (No JSON version)

Generates comprehensive visualization of cell calling results using directory structure.
"""

import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import scanpy as sc

# Set plot style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

def load_cell_calling_results(cell_calling_dir, sample_id):
    """Load cell calling results from saved files."""
    results_dir = Path(cell_calling_dir)
    
    # Load summary results
    summary_file = results_dir / "results.tsv"
    summary_df = pd.read_csv(summary_file, sep='\\t')
    
    # Load EmptyDrops detailed results
    emptydrops_file = results_dir / "emptydrops_results.tsv"
    emptydrops_df = pd.read_csv(emptydrops_file, sep='\\t')
    
    # Reconstruct methods_results dictionary
    methods_results = {}
    for _, row in summary_df.iterrows():
        method = row['method']
        threshold = row['threshold_used']
        
        # Load cell barcodes for this method
        barcode_file = results_dir / f"{sample_id}_{method}_cell_barcodes.txt"
        with open(barcode_file, 'r') as f:
            cell_barcodes = set(line.strip() for line in f)
        
        methods_results[method] = {
            'cell_barcodes': cell_barcodes,
            'threshold': threshold,
            'n_cells': row['n_cells_called'],
            'umis_in_cells_pct': row['umis_in_cells_pct']
        }
    
    return methods_results, emptydrops_df


def load_read_statistics(read_stats_file):
    """Load read statistics from TSV file."""
    df = pd.read_csv(read_stats_file, sep='\t')
    if df.empty:
        raise ValueError(f"Read statistics file is empty: {read_stats_file}")
    # Return the first row as a dictionary
    return df.iloc[0].to_dict()


def create_umi_distribution_plot_by_method(adata, methods_results, sample_id,
                                          output_dir, plot_dir, pool, source, processing, read_stats=None):
    """Create separate UMI distribution plots for each cell calling method."""
    # Calculate total UMIs per barcode
    total_counts = np.array(adata.X.sum(axis=1)).flatten()
    
    # Get all methods
    methods = sorted(methods_results.keys())
    n_methods = len(methods)
    
    # Create figure with subplots (3x3 grid to accommodate 7 methods)
    n_cols = 3
    n_rows = (n_methods + n_cols - 1) // n_cols  # Ceiling division
    
    fig = plt.figure(figsize=(18, 6 * n_rows))
    
    # Create consistent bins for all plots
    bins = np.logspace(0, np.log10(max(total_counts)), 50)
    
    # Create a subplot for each method
    for idx, method in enumerate(methods):
        ax = plt.subplot(n_rows, n_cols, idx + 1)
        
        # Get cell barcodes for this method
        cell_barcodes = methods_results[method]['cell_barcodes']
        n_cells = methods_results[method]['n_cells']
        
        # Create masks for cells vs non-cells using vectorized operations
        barcode_series = pd.Series(range(len(adata.obs_names)), index=adata.obs_names)
        cell_indices = barcode_series[barcode_series.index.isin(cell_barcodes)].values
        cell_mask = np.zeros(len(adata.obs_names), dtype=bool)
        cell_mask[cell_indices] = True
        cell_umis = total_counts[cell_mask]
        noncell_umis = total_counts[~cell_mask]
        
        # Plot histograms
        ax.hist(noncell_umis, bins=bins, alpha=0.5, 
                label=f'Non-cells (n={len(noncell_umis):,})', 
                color='gray', density=True)
        ax.hist(cell_umis, bins=bins, alpha=0.7, 
                label=f'Cells (n={len(cell_umis):,})', 
                color='blue', density=True)
        
        # Add vertical line at threshold if applicable
        threshold = methods_results[method]['threshold']
        if threshold > 0 and method not in ['Expected_Cells', 'EmptyDrops_FDR0001', 'EmptyDrops_FDR001', 'EmptyDrops_FDR005']:
            ax.axvline(x=threshold, color='red', linestyle='--', alpha=0.7, 
                      label=f'Threshold: {threshold:.0f}')
        
        # Customize subplot
        ax.set_xscale('log')
        ax.set_xlabel('UMI Count', fontsize=10)
        ax.set_ylabel('Density', fontsize=10)
        ax.set_title(f'{method}\n{n_cells:,} cells', fontsize=12, fontweight='bold')
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(True, alpha=0.3)
        
        # Set consistent axis limits
        ax.set_xlim(1, max(total_counts))
    
    # Remove empty subplots if any
    for idx in range(n_methods, n_rows * n_cols):
        fig.delaxes(plt.subplot(n_rows, n_cols, idx + 1))
    
    # Add total reads to title
    total_reads = read_stats['total_reads']
    title = f'UMI Distribution by Method - {sample_id}\nTotal Reads: {total_reads:,}'
    
    plt.suptitle(title, fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    # Save the plot
    plot_path = Path(plot_dir) / 'cell_calling' / f'{source}_{processing}' / sample_id / 'umi_distribution_by_method' / 'linear'
    plot_path.mkdir(parents=True, exist_ok=True)
    
    plot_filename = 'plot.png'
    full_plot_path = plot_path / plot_filename
    plt.savefig(full_plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nUMI distribution by method plot saved to: {full_plot_path}")
    return str(full_plot_path)


def create_barcode_rank_plot_by_method(adata, methods_results, sample_id, 
                                       output_dir, plot_dir, pool, source, processing, read_stats=None):
    """Create separate barcode rank plots for each cell calling method."""
    # Calculate total UMIs per barcode
    total_counts = np.array(adata.X.sum(axis=1)).flatten()
    total_counts_sorted = np.sort(total_counts)[::-1]
    ranks = np.arange(1, len(total_counts_sorted) + 1)
    
    # Get all methods
    methods = sorted(methods_results.keys())
    n_methods = len(methods)
    
    # Create figure with subplots (3x3 grid to accommodate 7 methods)
    n_cols = 3
    n_rows = (n_methods + n_cols - 1) // n_cols  # Ceiling division
    
    fig = plt.figure(figsize=(18, 6 * n_rows))
    
    # Create a subplot for each method
    for idx, method in enumerate(methods):
        ax = plt.subplot(n_rows, n_cols, idx + 1)
        
        # Plot all barcodes as background
        ax.loglog(ranks, total_counts_sorted, 'k-', alpha=0.2, linewidth=0.5, label='All barcodes')
        
        # Get cell barcodes for this method
        cell_barcodes = methods_results[method]['cell_barcodes']
        n_cells = methods_results[method]['n_cells']
        threshold = methods_results[method]['threshold']
        umis_in_cells_pct = methods_results[method]['umis_in_cells_pct']
        
        if len(cell_barcodes) > 0:
            # Get indices of cells
            barcode_to_idx = {bc: i for i, bc in enumerate(adata.obs_names)}
            cell_indices = [barcode_to_idx[bc] for bc in cell_barcodes if bc in barcode_to_idx]
            
            if cell_indices:
                cell_counts = total_counts[cell_indices]
                # Calculate ranks for the cells
                cell_ranks = []
                for count in cell_counts:
                    # Find rank of this count value
                    rank = np.searchsorted(-total_counts_sorted, -count) + 1
                    cell_ranks.append(rank)
                cell_ranks = np.array(cell_ranks)
                
                # Plot cells in color
                ax.scatter(cell_ranks, cell_counts, alpha=0.8, s=4, 
                          label=f'Cells (n={n_cells:,})', color='red', zorder=5)
                
                # Add threshold line if applicable
                if threshold > 0 and method not in ['Expected_Cells']:
                    ax.axhline(y=threshold, color='blue', linestyle='--', alpha=0.7, 
                              label=f'Threshold: {threshold:.0f}')
                    
                    # For BarcodeRanks methods, also add vertical line at the rank position
                    if method in ['BarcodeRanks_Knee', 'BarcodeRanks_Inflection']:
                        # Find the rank where this threshold is reached
                        # Use searchsorted to find where threshold would be inserted in sorted array
                        threshold_rank = np.searchsorted(-total_counts_sorted, -threshold) + 1
                        ax.axvline(x=threshold_rank, color='blue', linestyle='--', alpha=0.7,
                                  label=f'Rank: {threshold_rank:,}')
        
        # Customize subplot
        ax.set_xlabel('Barcode Rank', fontsize=10)
        ax.set_ylabel('UMI Counts', fontsize=10)
        ax.set_title(f'{method}\n{n_cells:,} cells ({umis_in_cells_pct:.1f}% UMIs)', 
                    fontsize=12, fontweight='bold')
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(True, alpha=0.3)
        
        # Set consistent axis limits
        ax.set_xlim(1, len(ranks))
        ax.set_ylim(1, max(total_counts_sorted))
    
    # Remove empty subplots if any
    for idx in range(n_methods, n_rows * n_cols):
        fig.delaxes(plt.subplot(n_rows, n_cols, idx + 1))
    
    # Add total reads to title
    total_reads = read_stats['total_reads']
    plt.suptitle(f'Barcode Rank Plots by Method - {sample_id}\nTotal Reads: {total_reads:,}', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    # Save the plot
    plot_path = Path(plot_dir) / 'cell_calling' / f'{source}_{processing}' / sample_id / 'barcode_rank_by_method' / 'linear'
    plot_path.mkdir(parents=True, exist_ok=True)
    
    plot_filename = 'plot.png'
    full_plot_path = plot_path / plot_filename
    plt.savefig(full_plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nBarcode rank by method plot saved to: {full_plot_path}")
    return str(full_plot_path)


def create_cell_calling_summary_plot(adata, methods_results, emptydrops_df, sample_id, 
                                    output_dir, plot_dir, pool, source, processing, read_stats=None):
    """Create cell calling summary plots (comparisons and statistics)."""
    # Calculate total UMIs per barcode
    total_counts = np.array(adata.X.sum(axis=1)).flatten()
    total_counts_sorted = np.sort(total_counts)[::-1]
    ranks = np.arange(1, len(total_counts_sorted) + 1)
    
    # Create figure with 4 subplots (2x2 grid)
    fig = plt.figure(figsize=(15, 10))
    
    # ===== Plot 1: Cell Count Comparison =====
    ax1 = plt.subplot(2, 2, 1)
    
    methods = sorted(methods_results.keys())
    cell_counts = [methods_results[m]['n_cells'] for m in methods]
    colors_bar = plt.cm.Set3(np.linspace(0, 1, len(methods)))
    
    bars = ax1.bar(range(len(methods)), cell_counts, color=colors_bar)
    ax1.set_xlabel('Method')
    ax1.set_ylabel('Number of Cells Called')
    ax1.set_title('Cell Calling Method Comparison', fontsize=14, fontweight='bold')
    ax1.set_xticks(range(len(methods)))
    ax1.set_xticklabels([m.replace('_', '\n') for m in methods], rotation=45, ha='right')
    
    # Add value labels with thresholds
    for i, (bar, count) in enumerate(zip(bars, cell_counts)):
        method = methods[i]
        threshold = methods_results[method]['threshold']
        
        # Format threshold based on method type
        if method.startswith('EmptyDrops'):
            threshold_str = f"FDR {threshold}"
        elif threshold > 0:
            threshold_str = f"Thr: {threshold:.0f}"
        else:
            threshold_str = ""
        
        # Add cell count label
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(cell_counts)*0.01,
                f'{count:,}', ha='center', va='bottom', fontsize=9)
        
        # Add threshold label below if exists
        if threshold_str:
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height()/2,
                    threshold_str, ha='center', va='center', fontsize=8,
                    color='darkred', weight='bold')
    
    # ===== Plot 2: Percentage of UMIs in Cells =====
    ax2 = plt.subplot(2, 2, 2)
    
    umis_percentages = [methods_results[m]['umis_in_cells_pct'] for m in methods]
    bars = ax2.bar(range(len(methods)), umis_percentages, color=colors_bar)
    ax2.set_xlabel('Method')
    ax2.set_ylabel('Percentage of UMIs in Cells')
    ax2.set_title('UMI Recovery by Method', fontsize=14, fontweight='bold')
    ax2.set_xticks(range(len(methods)))
    ax2.set_xticklabels([m.replace('_', '\n') for m in methods], rotation=45, ha='right')
    ax2.set_ylim(0, 100)
    
    # Add value labels
    for bar, pct in zip(bars, umis_percentages):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{pct:.1f}%', ha='center', va='bottom', fontsize=9)
    
    # ===== Plot 3: Summary Statistics =====
    ax3 = plt.subplot(2, 2, 3)
    ax3.axis('off')
    
    # Create summary text
    summary_text = f"Sample: {sample_id}\n"
    summary_text += f"Pool: {pool}\n"
    summary_text += f"Source: {source}, Processing: {processing}\n\n"
    summary_text += f"Total barcodes: {len(total_counts):,}\n"
    summary_text += f"Total UMIs: {np.sum(total_counts):,}\n\n"
    
    # Add method details
    summary_text += "Method Summary:\n"
    for method in methods[:7]:  # Show all methods
        n_cells = methods_results[method]['n_cells']
        umis_pct = methods_results[method]['umis_in_cells_pct']
        summary_text += f"{method}: {n_cells:,} cells ({umis_pct:.1f}% UMIs)\n"
    
    ax3.text(0.1, 0.9, summary_text, transform=ax3.transAxes,
            fontsize=11, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # ===== Plot 4: Read Recovery Explanation =====
    ax4 = plt.subplot(2, 2, 4)
    ax4.axis('off')
    
    # Create explanation text
    explanation_text = "UMI Recovery Calculation\n\n"
    explanation_text += "UMI Recovery % = \n"
    explanation_text += "Total UMIs in called cells\n"
    explanation_text += "─────────────────────────  × 100\n"
    explanation_text += "Total UMIs in all barcodes\n\n"
    explanation_text += "Measures the percentage of UMIs\n"
    explanation_text += "from true cells vs empty droplets."
    
    ax4.text(0.5, 0.5, explanation_text, transform=ax4.transAxes,
            fontsize=10, verticalalignment='center', horizontalalignment='center',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
    
    # Add title with total reads
    total_reads = read_stats['total_reads']
    plt.suptitle(f'Cell Calling Summary - {sample_id}\nTotal Reads: {total_reads:,}', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    # Save using directory structure matching Snakemake expectations
    plot_path = Path(plot_dir) / 'cell_calling' / f'{source}_{processing}' / sample_id / 'cell_calling_summary' / 'linear'
    plot_path.mkdir(parents=True, exist_ok=True)
    
    plot_filename = 'plot.png'
    full_plot_path = plot_path / plot_filename
    plt.savefig(full_plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nCell calling summary plot saved to: {full_plot_path}")
    return str(full_plot_path)


def main():
    parser = argparse.ArgumentParser(description='Generate cell calling plots from analysis results')
    parser.add_argument('--h5ad_file', required=True, help='Path to h5ad file')
    parser.add_argument('--cell_calling_dir', required=True, help='Cell calling results directory')
    parser.add_argument('--read_stats', required=True, help='Read statistics TSV file')
    parser.add_argument('--sample-id', required=True, help='Full sample ID (format: pool:sample)')
    parser.add_argument('--plot_dir', required=True, help='Plot output directory')
    parser.add_argument('--source', required=True, help='Source type (main/undetermined/all)')
    parser.add_argument('--processing', required=True, help='Processing type (raw/recovered/merged)')
    
    args = parser.parse_args()
    
    # Extract pool and sample from sample_id
    pool, sample = args.sample_id.split(':')
    
    # Load h5ad file
    print(f"Loading h5ad file: {args.h5ad_file}")
    adata = sc.read_h5ad(args.h5ad_file)
    print(f"Loaded matrix: {adata.n_obs} barcodes x {adata.n_vars} genes")
    
    # Load cell calling results
    print(f"Loading cell calling results from: {args.cell_calling_dir}")
    methods_results, emptydrops_df = load_cell_calling_results(args.cell_calling_dir, args.sample_id)
    print(f"Loaded results for {len(methods_results)} methods")
    
    # Load read statistics
    print(f"Loading read statistics from: {args.read_stats}")
    read_stats = load_read_statistics(args.read_stats)
    total_reads = read_stats['total_reads']
    print(f"Total reads: {total_reads:,}")
    
    # Generate plots
    print("Generating cell calling plots...")
    
    # Generate the barcode rank plots by method
    print("Generating barcode rank plots by method...")
    by_method_path = create_barcode_rank_plot_by_method(
        adata, methods_results, args.sample_id,
        args.cell_calling_dir, args.plot_dir,
        pool, args.source, args.processing, read_stats
    )
    
    # Generate the UMI distribution plots by method
    print("Generating UMI distribution plots by method...")
    umi_dist_path = create_umi_distribution_plot_by_method(
        adata, methods_results, args.sample_id,
        args.cell_calling_dir, args.plot_dir,
        pool, args.source, args.processing, read_stats
    )
    
    # Generate the summary plots
    print("Generating cell calling summary plots...")
    summary_path = create_cell_calling_summary_plot(
        adata, methods_results, emptydrops_df, args.sample_id,
        args.cell_calling_dir, args.plot_dir, 
        pool, args.source, args.processing, read_stats
    )
    
    print(f"✅ Cell calling plots completed successfully!")
    
    # Clean up memory
    del adata
    import gc
    gc.collect()


if __name__ == "__main__":
    main()