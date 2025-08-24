#!/usr/bin/env python3
"""
UMI saturation analysis using FASTQ subsampling and kb count quantification.
Provides accurate saturation curves by subsampling at the sequencing level.
"""

import argparse
import subprocess
import sys
import os
import tempfile
import time
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import gc
from pathlib import Path
from scipy import sparse
import scanpy as sc
from datetime import datetime

# Set global matplotlib parameters
plt.rcParams['savefig.dpi'] = 100

# Timing utilities
def timer(func):
    """Decorator to time function execution"""
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        elapsed = time.time() - start
        print(f"    ⏱️  {func.__name__}: {elapsed:.2f}s")
        return result
    return wrapper

class Timer:
    """Context manager for timing code blocks"""
    def __init__(self, name):
        self.name = name
        self.start = None
    
    def __enter__(self):
        self.start = time.time()
        return self
    
    def __exit__(self, *args):
        elapsed = time.time() - self.start
        print(f"    ⏱️  {self.name}: {elapsed:.2f}s")
# Import shared guide utility functions
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.pipeline_utils import get_paired_sample
from scripts.guide_analysis import (
    perform_guide_mixture_analysis, 
    plot_mixture_analysis,
    calculate_guide_metrics_for_cells,
    map_guide_counts_to_gex_barcodes
)


def get_sample_pool_from_id(config, sample_id):
    """Get the pool name for a sample ID."""
    sample_info_file = config['sample_info_file']
    if not os.path.exists(sample_info_file):
        raise FileNotFoundError(f"Sample info file not found: {sample_info_file}")
    
    sample_df = pd.read_excel(sample_info_file)
    sample_row = sample_df[sample_df['sample_id'] == sample_id]
    
    if sample_row.empty:
        raise ValueError(f"Sample {sample_id} not found in sample_info.xlsx")
    
    return sample_row.iloc[0]['pool']


def load_read_statistics(read_stats_file):
    """Load read statistics from TSV file."""
    df = pd.read_csv(read_stats_file, sep='\t')
    if df.empty:
        raise ValueError(f"Read statistics file is empty: {read_stats_file}")
    return df.iloc[0]['total_reads']

def load_adata_from_filtered(matrix_dir, is_guide=False):
    """Load AnnData from filtered h5ad files."""
    try:
        h5ad_file = Path(matrix_dir) / "adata.h5ad"
        
        if not h5ad_file.exists():
            raise FileNotFoundError(f"H5AD file not found: {h5ad_file}")
            
        # Load h5ad file (maintains correct cells x genes format)
        adata = sc.read_h5ad(h5ad_file)
        
        # Sum layers to get total counts for GEX samples
        if not is_guide and 'mature' in adata.layers and 'nascent' in adata.layers and 'ambiguous' in adata.layers:
            adata.X = adata.layers['mature'] + adata.layers['nascent'] + adata.layers['ambiguous']
        
        return adata
        
    except Exception as e:
        print(f"Error loading h5ad file: {e}", file=sys.stderr)
        raise

def load_adata_all_barcodes(matrix_dir, is_guide=False):
    """Load all barcodes from unfiltered AnnData (no filtering)."""
    try:
        h5ad_file = Path(matrix_dir) / "counts_unfiltered" / "adata.h5ad"
        
        if not h5ad_file.exists():
            raise FileNotFoundError(f"Unfiltered H5AD file not found: {h5ad_file}")
            
        print(f"Loading unfiltered data from {h5ad_file}")
        
        # Load all data into memory
        adata = sc.read_h5ad(h5ad_file)
        print(f"  Loaded {adata.n_obs} barcodes x {adata.n_vars} features")
        
        # Sum layers to get total counts for GEX samples
        if not is_guide and 'mature' in adata.layers and 'nascent' in adata.layers and 'ambiguous' in adata.layers:
            adata.X = adata.layers['mature'] + adata.layers['nascent'] + adata.layers['ambiguous']
            print(f"  Summed mature + nascent + ambiguous layers")
        
        return adata
        
    except Exception as e:
        print(f"Error loading unfiltered h5ad file: {e}", file=sys.stderr)
        raise

def calculate_umi_stats(adata, cell_barcodes_dict=None, is_guide=False, guide_cutoffs=None):
    """Calculate UMI statistics for an AnnData object, optionally for specific cell subsets."""
    stats = {}
    
    # Total UMIs per barcode
    umis_per_barcode = np.array(adata.X.sum(axis=1)).flatten()
    
    # Total features detected per barcode (genes or guides)
    features_per_barcode = np.array((adata.X > 0).sum(axis=1)).flatten()
    
    # Total unique UMIs
    total_umis = np.sum(umis_per_barcode)
    
    # Number of barcodes with UMIs
    barcodes_with_umis = np.sum(umis_per_barcode > 0)
    
    feature_type = 'guides' if is_guide else 'genes'
    
    stats['total_umis'] = total_umis
    stats['total_barcodes'] = adata.n_obs
    stats['barcodes_with_umis'] = barcodes_with_umis
    stats['mean_umis_per_barcode'] = np.mean(umis_per_barcode)
    stats['median_umis_per_barcode'] = np.median(umis_per_barcode)
    stats[f'mean_{feature_type}_per_barcode'] = np.mean(features_per_barcode)
    stats[f'median_{feature_type}_per_barcode'] = np.median(features_per_barcode)
    
    # Calculate stats for each cell calling method if provided
    if cell_barcodes_dict:
        for method, cell_barcode_list in cell_barcodes_dict.items():
            # Vectorized barcode intersection
            cell_mask = adata.obs_names.isin(cell_barcode_list)
            valid_barcodes = adata.obs_names[cell_mask].tolist()
            
            if len(valid_barcodes) == 0:
                continue
            
            # For guide samples, calculate metrics for all cutoffs
            if is_guide and guide_cutoffs:
                # Create a subset adata with only the guide counts for these cells
                guide_adata_subset = sc.AnnData(
                    X=adata.obsm['guide_counts'][cell_mask],
                    obs=adata.obs.loc[valid_barcodes].copy()
                )
                
                with Timer(f"    Guide metrics for {method} ({len(valid_barcodes)} cells)"):
                    guide_metrics = calculate_guide_metrics_for_cells(
                        guide_adata=guide_adata_subset,
                        cell_barcodes=valid_barcodes,
                        guide_cutoffs=guide_cutoffs,
                        obs_data=adata.obs.loc[valid_barcodes],  # Pass subset of obs matching the cells
                        method_name=method,
                        stratify_by='sample'
                    )
                
                stats[f'{method}_n_cells'] = guide_metrics['n_cells_total']
                stats[f'{method}_mean_umis_per_cell'] = guide_metrics['guide_umis_per_cell']
                stats[f'{method}_median_umis_per_cell'] = guide_metrics['guide_umis_per_cell']
                
                cell_umis = umis_per_barcode[cell_mask]
                stats[f'{method}_total_umis_in_cells'] = np.sum(cell_umis)
                stats[f'{method}_fraction_umis_in_cells'] = np.sum(cell_umis) / total_umis
                
                for key, value in guide_metrics.items():
                    if key.startswith('mean_') or key.startswith('median_') or key.startswith('fraction_') or key.startswith('umis_'):
                        stats[f'{method}_{key}'] = value
                    
            else:
                # For GEX samples, use barcode-based indexing
                cell_umis = umis_per_barcode[cell_mask]
                cell_features = features_per_barcode[cell_mask]
                
                stats[f'{method}_n_cells'] = len(valid_barcodes)
                stats[f'{method}_mean_umis_per_cell'] = np.mean(cell_umis)
                stats[f'{method}_median_umis_per_cell'] = np.median(cell_umis)
                stats[f'{method}_mean_{feature_type}_per_cell'] = np.mean(cell_features)
                stats[f'{method}_median_{feature_type}_per_cell'] = np.median(cell_features)
                stats[f'{method}_total_umis_in_cells'] = np.sum(cell_umis)
                stats[f'{method}_fraction_umis_in_cells'] = np.sum(cell_umis) / total_umis
    
    # Clean up intermediate arrays
    del umis_per_barcode, features_per_barcode
    
    return stats


def run_saturation_analysis(config, sample_id, saturation_points, cell_barcodes_dict, source, processing, scratch_dir, is_guide=False, calculate_gmm=False, guide_cutoffs=None, gmm_plot_dir=None, n_threads=4):
    """Run UMI saturation analysis by loading pre-generated subsampled matrices."""
    analysis_start = time.time()
    print(f"\n{'='*60}")
    print(f"Starting saturation analysis for {sample_id}")
    print(f"{'='*60}")
    
    saturation_data = []
    gmm_thresholds_by_depth = {}  # Store GMM thresholds at each depth
    
    # Process each saturation fraction
    for fraction in saturation_points:
        print(f"\nProcessing fraction {fraction}...")
        fraction_start = time.time()
        
        if is_guide:
            sample_info_file = config['sample_info_file']
            paired_gex_sample = get_paired_sample(sample_id, 'guide', sample_info_file)
            
            if fraction == 1.0:
                gex_dir = Path(f"{scratch_dir}/{paired_gex_sample}/kb_all_{source}_{processing}")
                with Timer("Loading GEX data"):
                    gex_adata = load_adata_all_barcodes(gex_dir, False)
                
                guide_dir = Path(f"{scratch_dir}/{sample_id}/kb_guide_{source}_{processing}")
                with Timer("Loading guide data"):
                    guide_adata = load_adata_all_barcodes(guide_dir, True)
                
                with Timer("Mapping guide to GEX barcodes"):
                    mapped_counts = map_guide_counts_to_gex_barcodes(guide_adata, gex_adata.obs_names)
                
                adata = gex_adata.copy()
                adata.obsm['guide_counts'] = mapped_counts
                print(f"  Mapped guide data to {len(gex_adata)} GEX barcodes")
            else:
                gex_subsample_dir = f"{scratch_dir}/tmp/umi_sat_{paired_gex_sample}-{fraction}/kb_all"
                with Timer("Loading GEX data"):
                    gex_adata = load_adata_all_barcodes(Path(gex_subsample_dir), False)
                
                guide_subsample_dir = f"{scratch_dir}/tmp/umi_sat_{sample_id}-{fraction}/kb_guide"
                with Timer("Loading guide data"):
                    guide_adata = load_adata_all_barcodes(Path(guide_subsample_dir), True)
                
                with Timer("Mapping guide to GEX barcodes"):
                    mapped_counts = map_guide_counts_to_gex_barcodes(guide_adata, gex_adata.obs_names)
                
                adata = gex_adata.copy()
                adata.obsm['guide_counts'] = mapped_counts
                print(f"  Mapped guide data to {len(gex_adata)} GEX barcodes")
        else:
            if fraction == 1.0:
                kb_result_dir = Path(f"{scratch_dir}/{sample_id}/kb_all_{source}_{processing}")
                adata = load_adata_all_barcodes(kb_result_dir, False)
                print(f"  Using existing full matrix from {kb_result_dir}")
            else:
                subsample_dir = f"{scratch_dir}/tmp/umi_sat_{sample_id}-{fraction}/kb_all"
                adata = load_adata_all_barcodes(Path(subsample_dir), False)
                print(f"  Loaded subsampled matrix from {subsample_dir}")
        
        if calculate_gmm and is_guide and guide_cutoffs:
            with Timer("GMM threshold calculation"):
                # Extract posterior levels from guide cutoffs
                posterior_levels = []
                for cutoff in guide_cutoffs:
                    if isinstance(cutoff, str) and (cutoff.startswith('gmm_') or cutoff.startswith('posterior_')):
                        if cutoff.startswith('gmm_'):
                            level = int(cutoff[4:])
                        else:
                            level = int(cutoff[10:])
                        if level not in posterior_levels:
                            posterior_levels.append(level)
                
                # If no GMM cutoffs specified, skip GMM calculation
                if not posterior_levels:
                    posterior_levels = [50]  # Default fallback
                
                min_umi_threshold = config.get('guide_mixture_model', {}).get('min_umi_threshold', 2)
                
                # Calculate GMM thresholds for each method at this depth
                # Store plot data for default method to create combined plot later
                default_method = config['cell_calling']['default_method']
                
                for method_name, method_cells in cell_barcodes_dict.items():
                    cell_mask = adata.obs.index.isin(method_cells)
                    adata_method = adata[cell_mask]
                    
                    if len(adata_method) == 0:
                        continue
                    
                    # Calculate GMM thresholds for this method at this depth
                    with Timer(f"  GMM analysis for {method_name}"):
                        mixture_result = perform_guide_mixture_analysis(
                    guide_counts_data=adata_method.obsm['guide_counts'],
                    min_umi_threshold=min_umi_threshold,
                    subsample_size=5000,
                    n_threads=n_threads,
                    posterior_levels=posterior_levels
                        )
                    
                    # Extract and store thresholds
                    if mixture_result is not None:
                        method_thresholds = mixture_result['posterior_thresholds']
                        
                        # Store in depth-specific structure
                        if fraction not in gmm_thresholds_by_depth:
                            gmm_thresholds_by_depth[fraction] = {}
                        gmm_thresholds_by_depth[fraction][method_name] = method_thresholds
                        
                        # Only store plot data for default method (will plot all depths together later)
                        if method_name == default_method and gmm_plot_dir is not None:
                            if 'gmm_plot_data' not in locals():
                                gmm_plot_data = {}
                            gmm_plot_data[fraction] = {
                                'data': mixture_result['pooled_data'],
                                'params': mixture_result['mixture_params']
                            }
                        
                        # Add method-specific threshold columns to adata.obs
                        for level in posterior_levels:
                            if level in method_thresholds:
                                col_name = f'gmm_{level}_sample_{method_name}'
                                threshold_value = method_thresholds[level]
                                
                                # Initialize column if it doesn't exist
                                if col_name not in adata.obs.columns:
                                    adata.obs[col_name] = np.nan
                                
                                # Set threshold for this method's cells
                                adata.obs.loc[cell_mask, col_name] = threshold_value
        
        # Calculate statistics including cell-specific metrics
        with Timer("Calculating UMI statistics"):
            stats = calculate_umi_stats(adata, cell_barcodes_dict, is_guide, guide_cutoffs)
        
        # Add GMM thresholds to stats if they were calculated
        if calculate_gmm and is_guide and guide_cutoffs:
            for method_name in cell_barcodes_dict.keys():
                if fraction in gmm_thresholds_by_depth and method_name in gmm_thresholds_by_depth[fraction]:
                    method_thresholds = gmm_thresholds_by_depth[fraction][method_name]
                    for level in method_thresholds:
                        if method_name == config['cell_calling']['default_method']:
                            stats[f'gmm_{level}_threshold'] = method_thresholds[level]
        
        # Add fraction info
        stats['sample_id'] = sample_id
        stats['subsample_fraction'] = fraction
        stats['n_features'] = adata.n_vars
        
        saturation_data.append(stats)
        
        print(f"  Fraction {fraction}: {stats['total_umis']} UMIs, "
              f"{stats['barcodes_with_umis']} active barcodes")
        
        # Print cell-specific stats if available
        if cell_barcodes_dict:
            for method in cell_barcodes_dict.keys():
                n_cells = stats.get(f'{method}_n_cells', 0)
                mean_umis = stats.get(f'{method}_mean_umis_per_cell', 0)
                print(f"    {method}: {n_cells} cells, {mean_umis:.1f} avg UMIs/cell")
        
        # Clean up memory after processing each fraction
        del adata
        gc.collect()
        
        fraction_elapsed = time.time() - fraction_start
        print(f"  ✅ Fraction {fraction} completed in {fraction_elapsed:.2f}s")
    
    # Create combined GMM plot for default method if we have the data
    if calculate_gmm and gmm_plot_dir is not None and 'gmm_plot_data' in locals() and gmm_plot_data:
        print(f"Generating combined GMM plot for {default_method}...")
        
        # Create figure with subplots for each depth
        n_depths = len(gmm_plot_data)
        fig, axes = plt.subplots(n_depths, 2, figsize=(14, 6 * n_depths))
        
        # Handle single depth case
        if n_depths == 1:
            axes = axes.reshape(1, -1)
        
        # Sort by depth for consistent ordering
        sorted_depths = sorted(gmm_plot_data.keys())
        
        for idx, fraction in enumerate(sorted_depths):
            plot_data = gmm_plot_data[fraction]
            ax_pair = (axes[idx, 0], axes[idx, 1])
            
            # Use plot_mixture_analysis with provided axes
            plot_mixture_analysis(
                data=plot_data['data'],
                params=plot_data['params'],
                output_path=None,  # Not saving individual plots
                title_prefix=f"{sample_id} - {default_method} - {int(fraction*100)}% depth",
                min_umi_threshold=min_umi_threshold,
                posterior_levels=posterior_levels,
                axes=ax_pair
            )
            
            # Add a clear fraction label on the left side of each row
            axes[idx, 0].text(-0.18, 0.5, f'{int(fraction*100)}%\nDepth', 
                            transform=axes[idx, 0].transAxes,
                            fontsize=12, fontweight='bold', 
                            verticalalignment='center',
                            horizontalalignment='right',
                            color='darkblue')
        
        # Overall title with more descriptive text
        fig.suptitle(f'GMM Mixture Model Fits at Different Sequencing Depths - {default_method}', 
                    fontsize=16, y=1.01)
        
        # Save combined plot
        plot_subdir = gmm_plot_dir / f"gmm_{default_method}_combined" / "linear"
        plot_subdir.mkdir(parents=True, exist_ok=True)
        plot_path = plot_subdir / "plot.png"
        
        plt.tight_layout()
        plt.savefig(plot_path, bbox_inches='tight')
        plt.close()
        
        print(f"  Saved combined GMM plot: {plot_path}")
    
    # Print timing summary
    total_elapsed = time.time() - analysis_start
    print(f"\n{'='*60}")
    print(f"✅ Saturation analysis completed for {sample_id}")
    print(f"   Total time: {total_elapsed:.2f}s ({total_elapsed/60:.1f} minutes)")
    print(f"   Fractions processed: {len(saturation_points)}")
    print(f"   Average time per fraction: {total_elapsed/len(saturation_points):.2f}s")
    print(f"{'='*60}\n")
    
    if calculate_gmm:
        return saturation_data, gmm_thresholds_by_depth
    return saturation_data

def plot_saturation_curves(saturation_data, output_plot_path, sample_id, config, cell_methods=None, is_guide=False, pool=None, total_reads=None, gmm_thresholds_by_depth=None, guide_cutoffs=None):
    """Generate UMI saturation plots.
    
    Returns:
        dict: Plot metadata
    """
    df = pd.DataFrame(saturation_data)
    
    if len(df) == 0:
        print("No saturation data to plot")
        return
    
    # Determine feature type for labels
    feature_type = 'guides' if is_guide else 'genes'
    feature_type_cap = 'Guide' if is_guide else 'Gene'
    
    # Create plots with cell-specific metrics if available
    if cell_methods:
        if is_guide:
            # For guide samples: 2 basic + 1 UMI + 3 cutoff-dependent × 2 + 1 GMM threshold = 8 rows
            fig, axes = plt.subplots(8, 2, figsize=(15, 48))
        else:
            fig, axes = plt.subplots(3, 2, figsize=(15, 18))
    else:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Calculate actual read counts from subsample fractions
    df['num_reads'] = df['subsample_fraction'] * total_reads
    x_data = df['num_reads'] / 1e6  # Convert to millions
    x_label = 'Number of Reads (millions)'
    
    # UMI saturation curve
    axes[0, 0].plot(x_data, df['total_umis'], 'bo-')
    axes[0, 0].set_xlabel(x_label)
    axes[0, 0].set_ylabel('Total UMIs')
    title = f'{feature_type_cap} UMI Saturation Curve'
    axes[0, 0].set_title(title)
    axes[0, 0].grid(True, alpha=0.3)
    
    # Active barcodes saturation
    axes[0, 1].plot(x_data, df['barcodes_with_umis'], 'ro-')
    axes[0, 1].set_xlabel(x_label)
    axes[0, 1].set_ylabel('Barcodes with UMIs')
    axes[0, 1].set_title('Active Barcodes Saturation')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Mean UMIs per barcode
    axes[1, 0].plot(x_data, df['mean_umis_per_barcode'], 'go-')
    axes[1, 0].set_xlabel(x_label)
    axes[1, 0].set_ylabel('Mean UMIs per Barcode')
    axes[1, 0].set_title('Mean UMIs per Barcode')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Mean features per barcode
    feature_col = f'mean_{feature_type}_per_barcode'
    axes[1, 1].plot(x_data, df[feature_col], 'mo-')
    axes[1, 1].set_xlabel(x_label)
    axes[1, 1].set_ylabel(f'Mean {feature_type_cap}s per Barcode')
    axes[1, 1].set_title(f'Mean {feature_type_cap}s per Barcode')
    axes[1, 1].grid(True, alpha=0.3)
    
    # Cell-specific plots if available
    if cell_methods:
        # Use tab10 colormap for consistent colors across all methods (up to 10 methods)
        colors = plt.cm.tab10(np.linspace(0, 0.9, len(cell_methods)))
        
        # Average UMIs per cell for each method
        for i, method in enumerate(cell_methods):
            mean_umis_col = f'{method}_mean_umis_per_cell'
            if mean_umis_col in df.columns:
                axes[2, 0].plot(x_data, df[mean_umis_col], 
                               'o-', color=colors[i], label=method)
        
        axes[2, 0].set_xlabel(x_label)
        axes[2, 0].set_ylabel('Mean UMIs per Cell')
        axes[2, 0].set_title(f'Mean {feature_type_cap} UMIs per Cell by Method')
        axes[2, 0].grid(True, alpha=0.3)
        axes[2, 0].legend()
        
        # For non-guide samples, show average genes per cell
        if not is_guide:
            for i, method in enumerate(cell_methods):
                mean_genes_col = f'{method}_mean_{feature_type}_per_cell'
                if mean_genes_col in df.columns:
                    axes[2, 1].plot(x_data, df[mean_genes_col], 
                                   'o-', color=colors[i], label=method)
            
            axes[2, 1].set_xlabel(x_label)
            axes[2, 1].set_ylabel(f'Mean {feature_type_cap}s per Cell')
            axes[2, 1].set_title(f'Mean {feature_type_cap}s per Cell by Method')
            axes[2, 1].grid(True, alpha=0.3)
            axes[2, 1].legend()
        else:
            # For guide samples, plot cutoff-independent metric by method (50% GMM)
            default_cutoff = "gmm50"  # Use 50% GMM as default
            for i, method in enumerate(cell_methods):
                mean_guides_col = f'{method}_mean_guides_per_cell_{default_cutoff}'
                if mean_guides_col not in df.columns:
                    # Fallback to first available cutoff
                    for col in df.columns:
                        if col.startswith(f'{method}_mean_guides_per_cell_'):
                            mean_guides_col = col
                            break
                
                if mean_guides_col in df.columns:
                    axes[2, 1].plot(x_data, df[mean_guides_col], 
                                   'o-', color=colors[i], label=method)
            
            axes[2, 1].set_xlabel(x_label)
            axes[2, 1].set_ylabel('Mean Guides per Cell')
            axes[2, 1].set_title('Mean Guides per Cell by Method (50% GMM)')
            axes[2, 1].grid(True, alpha=0.3)
            axes[2, 1].legend()
        
        # Add guide-specific plots if this is a guide sample
        if is_guide and cell_methods:
            # Row 3: Cutoff-independent metric (Mean UMIs per cell by method)
            for i, method in enumerate(cell_methods):
                umis_col = f'{method}_mean_umis_per_cell'
                if umis_col in df.columns:
                    axes[2, 0].plot(x_data, df[umis_col], 
                                   'o-', color=colors[i], label=method)
            
            axes[2, 0].set_xlabel(x_label)
            axes[2, 0].set_ylabel('Mean Guide UMIs per Cell')
            axes[2, 0].set_title('Guide UMIs per Cell by Method')
            axes[2, 0].grid(True, alpha=0.3)
            axes[2, 0].legend()
            
            # Row 4: Mean guides per cell - Method comparison (50% GMM)
            default_cutoff = "gmm50"
            for i, method in enumerate(cell_methods):
                col = f'{method}_mean_guides_per_cell_{default_cutoff}'
                if col not in df.columns:
                    for c in df.columns:
                        if c.startswith(f'{method}_mean_guides_per_cell_'):
                            col = c
                            break
                if col in df.columns:
                    axes[3, 0].plot(x_data, df[col], 'o-', color=colors[i], label=method)
            
            axes[3, 0].set_xlabel(x_label)
            axes[3, 0].set_ylabel('Mean Guides per Cell')
            axes[3, 0].set_title('Mean Guides per Cell - Method Comparison (50% GMM)')
            axes[3, 0].grid(True, alpha=0.3)
            axes[3, 0].legend()
            
            # Row 4: Mean guides per cell - Cutoff comparison (default method)
            default_method = config['cell_calling']['default_method']
            cutoff_colors = ['purple', 'orange', 'brown', 'pink', 'red']
            plot_idx = 0
            for cutoff_spec in guide_cutoffs if guide_cutoffs else []:
                if isinstance(cutoff_spec, int):
                    suffix = f"cutoff{cutoff_spec}"
                    label = f"≥{cutoff_spec} UMIs"
                elif isinstance(cutoff_spec, str) and (cutoff_spec.startswith('gmm_') or cutoff_spec.startswith('posterior_')):
                    level = cutoff_spec[4:] if cutoff_spec.startswith('gmm_') else cutoff_spec[10:]
                    suffix = f"gmm{level}"
                    label = f"{level}% GMM"
                else:
                    continue
                
                col = f'{default_method}_mean_guides_per_cell_{suffix}'
                if col in df.columns and plot_idx < len(cutoff_colors):
                    axes[3, 1].plot(x_data, df[col], 'o-', 
                                   color=cutoff_colors[plot_idx], label=label)
                    plot_idx += 1
            
            axes[3, 1].set_xlabel(x_label)
            axes[3, 1].set_ylabel('Mean Guides per Cell')
            axes[3, 1].set_title(f'Mean Guides per Cell - Cutoff Comparison ({default_method})')
            axes[3, 1].grid(True, alpha=0.3)
            axes[3, 1].legend()
            
            # Row 5: % Cells with guides - Method comparison (50% GMM)
            for i, method in enumerate(cell_methods):
                col = f'{method}_fraction_cells_with_guides_{default_cutoff}'
                if col not in df.columns:
                    for c in df.columns:
                        if c.startswith(f'{method}_fraction_cells_with_guides_'):
                            col = c
                            break
                if col in df.columns:
                    axes[4, 0].plot(x_data, df[col] * 100, 'o-', color=colors[i], label=method)
            
            axes[4, 0].set_xlabel(x_label)
            axes[4, 0].set_ylabel('% Cells with Guides')
            axes[4, 0].set_title('% Cells with Guides - Method Comparison (50% GMM)')
            axes[4, 0].grid(True, alpha=0.3)
            axes[4, 0].legend()
            
            # Row 5: % Cells with guides - Cutoff comparison (default method)
            plot_idx = 0
            for cutoff_spec in guide_cutoffs if guide_cutoffs else []:
                if isinstance(cutoff_spec, int):
                    suffix = f"cutoff{cutoff_spec}"
                    label = f"≥{cutoff_spec} UMIs"
                elif isinstance(cutoff_spec, str) and (cutoff_spec.startswith('gmm_') or cutoff_spec.startswith('posterior_')):
                    level = cutoff_spec[4:] if cutoff_spec.startswith('gmm_') else cutoff_spec[10:]
                    suffix = f"gmm{level}"
                    label = f"{level}% GMM"
                else:
                    continue
                
                col = f'{default_method}_fraction_cells_with_guides_{suffix}'
                if col in df.columns and plot_idx < len(cutoff_colors):
                    axes[4, 1].plot(x_data, df[col] * 100, 'o-', 
                                   color=cutoff_colors[plot_idx], label=label)
                    plot_idx += 1
            
            axes[4, 1].set_xlabel(x_label)
            axes[4, 1].set_ylabel('% Cells with Guides')
            axes[4, 1].set_title(f'% Cells with Guides - Cutoff Comparison ({default_method})')
            axes[4, 1].grid(True, alpha=0.3)
            axes[4, 1].legend()
            
            # Row 6: UMIs per guide - Method comparison (50% GMM)
            for i, method in enumerate(cell_methods):
                col = f'{method}_umis_per_guide_per_cell_{default_cutoff}'
                if col not in df.columns:
                    for c in df.columns:
                        if c.startswith(f'{method}_umis_per_guide_per_cell_'):
                            col = c
                            break
                if col in df.columns:
                    axes[5, 0].plot(x_data, df[col], 'o-', color=colors[i], label=method)
            
            axes[5, 0].set_xlabel(x_label)
            axes[5, 0].set_ylabel('Mean UMIs per Guide')
            axes[5, 0].set_title('UMIs per Guide - Method Comparison (50% GMM)')
            axes[5, 0].grid(True, alpha=0.3)
            axes[5, 0].legend()
            
            # Row 6: UMIs per guide - Cutoff comparison (default method)
            plot_idx = 0
            for cutoff_spec in guide_cutoffs if guide_cutoffs else []:
                if isinstance(cutoff_spec, int):
                    suffix = f"cutoff{cutoff_spec}"
                    label = f"≥{cutoff_spec} UMIs"
                elif isinstance(cutoff_spec, str) and (cutoff_spec.startswith('gmm_') or cutoff_spec.startswith('posterior_')):
                    level = cutoff_spec[4:] if cutoff_spec.startswith('gmm_') else cutoff_spec[10:]
                    suffix = f"gmm{level}"
                    label = f"{level}% GMM"
                else:
                    continue
                
                col = f'{default_method}_umis_per_guide_per_cell_{suffix}'
                if col in df.columns and plot_idx < len(cutoff_colors):
                    axes[5, 1].plot(x_data, df[col], 'o-', 
                                   color=cutoff_colors[plot_idx], label=label)
                    plot_idx += 1
            
            axes[5, 1].set_xlabel(x_label)
            axes[5, 1].set_ylabel('Mean UMIs per Guide')
            axes[5, 1].set_title(f'UMIs per Guide - Cutoff Comparison ({default_method})')
            axes[5, 1].grid(True, alpha=0.3)
            axes[5, 1].legend()
            
    
    # Add GMM threshold plot in the main figure if available
    if gmm_thresholds_by_depth and is_guide and cell_methods:
        # Row 7: GMM Thresholds vs Reads (single panel)
        fractions = sorted(gmm_thresholds_by_depth.keys())
        reads_per_fraction = [f * total_reads / 1e6 for f in fractions]
        
        # Get default method for threshold plotting
        default_method = config['cell_calling']['default_method']
        
        # Plot thresholds for each posterior level
        threshold_colors = ['blue', 'green', 'red', 'orange', 'purple']
        plot_idx = 0
        
        for level in sorted(set().union(*[list(thresholds.get(default_method, {}).keys()) 
                                        for thresholds in gmm_thresholds_by_depth.values()])):
            thresholds_at_level = []
            for fraction in fractions:
                if (fraction in gmm_thresholds_by_depth and 
                    default_method in gmm_thresholds_by_depth[fraction] and
                    level in gmm_thresholds_by_depth[fraction][default_method]):
                    thresholds_at_level.append(gmm_thresholds_by_depth[fraction][default_method][level])
                else:
                    thresholds_at_level.append(None)
            
            # Filter out None values
            valid_reads = [r for r, t in zip(reads_per_fraction, thresholds_at_level) if t is not None]
            valid_thresholds = [t for t in thresholds_at_level if t is not None]
            
            if valid_thresholds and plot_idx < len(threshold_colors):
                axes[6, 0].plot(valid_reads, valid_thresholds, 
                               marker='o', color=threshold_colors[plot_idx], 
                               label=f'{level}% posterior', linewidth=2)
                plot_idx += 1
        
        axes[6, 0].set_xlabel(x_label)
        axes[6, 0].set_ylabel('GMM Threshold (UMIs)')
        axes[6, 0].set_title(f'GMM Thresholds vs Sequencing Depth ({default_method})')
        axes[6, 0].grid(True, alpha=0.3)
        axes[6, 0].legend()
        
        # Hide the second panel in row 7
        axes[6, 1].set_visible(False)
        
        # Hide row 8 panels
        axes[7, 0].set_visible(False)
        axes[7, 1].set_visible(False)
    
    plt.figure(fig.number)  # Switch back to main figure
    plt.tight_layout()
    
    # Create clean directory structure from the output path
    # Expected: .../saturation/main_raw/gex/pool3/gex_batch_1/pool3_gex_batch_1_umi_saturation.png
    # Convert to: .../saturation/main_raw/gex/pool3/gex_batch_1/umi_saturation/plot.png
    output_dir = output_plot_path.parent
    
    # For guide samples, add guide_umi_saturation subdirectory
    if is_guide:
        clean_dir = output_dir / 'guide_umi_saturation' / 'linear'
    else:
        clean_dir = output_dir / 'umi_saturation' / 'linear'
    
    clean_dir.mkdir(parents=True, exist_ok=True)
    clean_path = clean_dir / 'plot.png'
    
    plt.savefig(clean_path, bbox_inches='tight')
    plt.close()
    
    # Update output_plot_path for metadata
    output_plot_path = clean_path
    
    # Create plot metadata
    plot_metadata = {
        'filename': output_plot_path.name,
        'path': str(output_plot_path),  # Will be adjusted in main
        'plot_type': 'saturation_multi_panel',
        'category': 'saturation',
        'sample_type': 'guide' if is_guide else 'gex',
        'display_name': f'{"Guide" if is_guide else "Gene Expression"} UMI Saturation Analysis',
        'description': f'UMI saturation curves showing sequencing depth effects on {"guide" if is_guide else "gene"} detection',
        'pool': pool,
        'sample': sample_id,
        'panels': []
    }
    
    # Add panel descriptions
    plot_metadata['panels'].extend([
        'UMI Saturation Curve',
        'Unique Features Detected'
    ])
    
    if cell_methods:
        plot_metadata['panels'].extend([
            'Mean UMIs per Cell by Method'
        ])
        
        if is_guide:
            plot_metadata['panels'].extend([
                'Mean Guides per Cell by Method (Cutoff = 1 UMI)',
                'Mean Guides per Cell by Method (Cutoff = 2 UMIs)'
            ])
        else:
            plot_metadata['panels'].append('Mean Genes per Cell by Method')
            
        plot_metadata['methods_included'] = cell_methods
        
        if is_guide and guide_cutoffs:
            plot_metadata['panels'].extend([
                'Guide Detection Saturation by Cutoff',
                'Cell Coverage Saturation by Cutoff'
            ])
            plot_metadata['guide_cutoffs'] = guide_cutoffs
    
    return plot_metadata

# DEPRECATED: Diversity analysis moved to cell QC pipeline
# # def analyze_umi_diversity_with_cells(config, sample_id, output_dir, cell_barcodes_dict, is_guide=False):
#     """Analyze UMI diversity and duplication rates with cell-specific metrics."""
#     # Load full count matrix from existing results (unfiltered)
#     try:
#         if is_guide:
#             # Guide samples use kb_guide directory 
#             kb_result_dir = Path(f"../analysis_results/{get_sample_pool_from_id(config, sample_id)}/{sample_id}/kb_guide_main_raw")
#         else:
#             # GEX samples use kb_all directory
#             kb_result_dir = Path(f"../analysis_results/{get_sample_pool_from_id(config, sample_id)}/{sample_id}/kb_all_main_raw")
#         
#         adata = load_adata_with_filtering(kb_result_dir, cell_barcodes_dict, is_guide)
#     except Exception as e:
#         raise RuntimeError(f"Critical error loading count matrix: {e}. Cannot proceed with diversity analysis.") from e
#     
#     # Calculate per-barcode statistics
#     umis_per_barcode = np.array(adata.X.sum(axis=1)).flatten()
#     features_per_barcode = np.array((adata.X > 0).sum(axis=1)).flatten()
#     
#     # UMI diversity (unique features / total UMIs) approximation
#     # This is simplified - true UMI diversity would require UMI sequences
#     umi_diversity = features_per_barcode / (umis_per_barcode + 1)  # +1 to avoid division by zero
#     
#     feature_type = 'guides' if is_guide else 'genes'
#     feature_type_cap = 'Guide' if is_guide else 'Gene'
#     
#     # Create summary statistics
#     diversity_stats = {
#         'sample_id': sample_id,
#         'mean_umi_diversity': np.mean(umi_diversity[umis_per_barcode > 0]),
#         'median_umi_diversity': np.median(umi_diversity[umis_per_barcode > 0]),
#         'mean_umis_per_barcode': np.mean(umis_per_barcode),
#         'median_umis_per_barcode': np.median(umis_per_barcode),
#         f'mean_{feature_type}_per_barcode': np.mean(features_per_barcode),
#         f'median_{feature_type}_per_barcode': np.median(features_per_barcode),
#         'total_barcodes': adata.n_obs,
#         'active_barcodes': np.sum(umis_per_barcode > 0)
#     }
#     
#     # Add cell-specific diversity stats
#     # Create pandas Series for vectorized barcode lookup
#     import pandas as pd
#     barcode_series = pd.Series(range(adata.n_obs), index=adata.obs_names)
#     
#     for method, cell_barcode_list in cell_barcodes_dict.items():
#         # Find indices of called cells using vectorized operations
#         cell_indices = barcode_series.reindex(cell_barcode_list).dropna().astype(int).values
#         
#         if len(cell_indices) > 0:
#             cell_umis = umis_per_barcode[cell_indices]
#             cell_features = features_per_barcode[cell_indices]
#             cell_diversity = umi_diversity[cell_indices]
#             
#             diversity_stats[f'{method}_mean_umi_diversity'] = np.mean(cell_diversity)
#             diversity_stats[f'{method}_median_umi_diversity'] = np.median(cell_diversity)
#             diversity_stats[f'{method}_mean_umis_per_cell'] = np.mean(cell_umis)
#             diversity_stats[f'{method}_median_umis_per_cell'] = np.median(cell_umis)
#             diversity_stats[f'{method}_mean_{feature_type}_per_cell'] = np.mean(cell_features)
#             diversity_stats[f'{method}_median_{feature_type}_per_cell'] = np.median(cell_features)
#         else:
#             diversity_stats[f'{method}_mean_umi_diversity'] = 0
#             diversity_stats[f'{method}_median_umi_diversity'] = 0
#             diversity_stats[f'{method}_mean_umis_per_cell'] = 0
#             diversity_stats[f'{method}_median_umis_per_cell'] = 0
#             diversity_stats[f'{method}_mean_{feature_type}_per_cell'] = 0
#             diversity_stats[f'{method}_median_{feature_type}_per_cell'] = 0
#     
#     # Plot UMI diversity with cell-specific overlays
#     fig, axes = plt.subplots(2, 3, figsize=(18, 12))
#     
#     # UMI vs features scatter plot with sampling for performance
#     active_mask = umis_per_barcode > 0
#     active_indices = np.where(active_mask)[0]
#     
#     # Sample points for plotting to avoid slowdown with 150K+ points
#     max_plot_points = 10000
#     if len(active_indices) > max_plot_points:
#         sampled_indices = np.random.choice(active_indices, size=max_plot_points, replace=False)
#         plot_label = f'Sampled barcodes (n={max_plot_points:,})'
#     else:
#         sampled_indices = active_indices
#         plot_label = f'All barcodes (n={len(sampled_indices):,})'
#     
#     axes[0, 0].scatter(umis_per_barcode[sampled_indices], features_per_barcode[sampled_indices], 
#                       alpha=0.3, s=1, color='gray', label=plot_label)
#     
#     # Overlay cells from different methods (sample cells for performance, show all methods)
#     colors = plt.cm.tab10(np.linspace(0, 1, len(cell_barcodes_dict)))  # Generate enough colors
#     max_cell_points_per_method = max_plot_points // max(len(cell_barcodes_dict), 1)  # Distribute points across all methods
#     
#     for i, (method, cell_barcode_list) in enumerate(cell_barcodes_dict.items()):
#         cell_indices = barcode_series.reindex(cell_barcode_list).dropna().astype(int).values
#         if len(cell_indices) > 0:
#             # Sample cells for plotting if too many
#             if len(cell_indices) > max_cell_points_per_method:
#                 cell_indices = np.random.choice(cell_indices, size=max_cell_points_per_method, replace=False)
#                 cell_label = f'{method} (n={max_cell_points_per_method:,})'
#             else:
#                 cell_label = f'{method} (n={len(cell_indices):,})'
#             
#             axes[0, 0].scatter(umis_per_barcode[cell_indices], features_per_barcode[cell_indices],
#                               alpha=0.6, s=2, color=colors[i], label=cell_label)
#     
#     axes[0, 0].set_xlabel('UMIs per Barcode')
#     axes[0, 0].set_ylabel(f'{feature_type_cap}s per Barcode')
#     axes[0, 0].set_title(f'UMIs vs {feature_type_cap}s per Barcode')
#     axes[0, 0].set_xscale('log')
#     axes[0, 0].set_yscale('log')
#     axes[0, 0].legend()
#     
#     # UMI diversity histogram
#     axes[0, 1].hist(umi_diversity[umis_per_barcode > 0], bins=50, alpha=0.5, 
#                    color='gray', label='All barcodes')
#     for i, (method, cell_barcode_list) in enumerate(cell_barcodes_dict.items()):
#         cell_indices = barcode_series.reindex(cell_barcode_list).dropna().astype(int).values
#         if len(cell_indices) > 0:
#             axes[0, 1].hist(umi_diversity[cell_indices], bins=50, alpha=0.7,
#                            color=colors[i], label=f'{method} cells')
#     
#     axes[0, 1].set_xlabel(f'UMI Diversity ({feature_type_cap}s/UMIs)')
#     axes[0, 1].set_ylabel('Number of Barcodes')
#     axes[0, 1].set_title(f'{feature_type_cap} UMI Diversity Distribution')
#     axes[0, 1].legend()
#     
#     # UMI distribution
#     axes[0, 2].hist(umis_per_barcode[umis_per_barcode > 0], bins=50, alpha=0.5,
#                    color='gray', label='All barcodes')
#     for i, (method, cell_barcode_list) in enumerate(cell_barcodes_dict.items()):
#         cell_indices = barcode_series.reindex(cell_barcode_list).dropna().astype(int).values
#         if len(cell_indices) > 0:
#             axes[0, 2].hist(umis_per_barcode[cell_indices], bins=50, alpha=0.7,
#                            color=colors[i], label=f'{method} cells')
#     
#     axes[0, 2].set_xlabel('UMIs per Barcode')
#     axes[0, 2].set_ylabel('Number of Barcodes')
#     axes[0, 2].set_title(f'{feature_type_cap} UMI Count Distribution')
#     axes[0, 2].set_xscale('log')
#     axes[0, 2].legend()
#     
#     # Features distribution
#     axes[1, 0].hist(features_per_barcode[features_per_barcode > 0], bins=50, alpha=0.5,
#                    color='gray', label='All barcodes')
#     for i, (method, cell_barcode_list) in enumerate(cell_barcodes_dict.items()):
#         cell_indices = barcode_series.reindex(cell_barcode_list).dropna().astype(int).values
#         if len(cell_indices) > 0:
#             axes[1, 0].hist(features_per_barcode[cell_indices], bins=50, alpha=0.7,
#                            color=colors[i], label=f'{method} cells')
#     
#     axes[1, 0].set_xlabel(f'{feature_type_cap}s per Barcode')
#     axes[1, 0].set_ylabel('Number of Barcodes')
#     axes[1, 0].set_title(f'{feature_type_cap} Count Distribution')
#     axes[1, 0].legend()
#     
#     # Average UMIs per cell by method (bar plot)
#     methods = list(cell_barcodes_dict.keys())
#     avg_umis = [diversity_stats.get(f'{method}_mean_umis_per_cell', 0) for method in methods]
#     
#     axes[1, 1].bar(methods, avg_umis, color=colors[:len(methods)])
#     axes[1, 1].set_ylabel('Average UMIs per Cell')
#     axes[1, 1].set_title(f'Average {feature_type_cap} UMIs per Cell by Method')
#     axes[1, 1].tick_params(axis='x', rotation=45)
#     
#     # Number of cells by method (bar plot)
#     n_cells = [len(cell_barcodes_dict[method]) for method in methods]
#     
#     axes[1, 2].bar(methods, n_cells, color=colors[:len(methods)])
#     axes[1, 2].set_ylabel('Number of Cells')
#     cell_source = ' (from paired GEX)' if is_guide else ''
#     axes[1, 2].set_title(f'Number of Cells by Method{cell_source}')
#     axes[1, 2].tick_params(axis='x', rotation=45)
#     
#     plt.tight_layout()
#     output_suffix = '_guide' if is_guide else ''
#     plt.savefig(Path(output_dir) / f'{sample_id}{output_suffix}_umi_diversity_analysis.png', 
#                 dpi=100, bbox_inches='tight')
#     plt.close()
#     
#     return diversity_stats

def load_cell_calling_results(cell_calling_dir, sample_id):
    """Load cell calling results from previous analysis."""
    results_file = Path(cell_calling_dir) / 'results.tsv'
    
    if not results_file.exists():
        raise FileNotFoundError(f"Cell calling results not found: {results_file}")
    
    results_df = pd.read_csv(results_file, sep='\t')
    
    # Load cell barcode lists for each method
    cell_barcodes = {}
    for method in results_df['method'].unique():
        barcode_file = Path(cell_calling_dir) / f'{sample_id}_{method}_cell_barcodes.txt'
        if barcode_file.exists():
            with open(barcode_file, 'r') as f:
                cell_barcodes[method] = [line.strip() for line in f]
        else:
            cell_barcodes[method] = []
    
    return results_df, cell_barcodes

# Removed duplicate function - now using pipeline_utils.get_paired_sample

def determine_sample_type(config, sample_id):
    """Determine if a sample is GEX or guide library."""
    sample_info_file = config['sample_info_file']
    if not os.path.exists(sample_info_file):
        raise FileNotFoundError(f"Sample info file not found: {sample_info_file}")
    
    sample_df = pd.read_excel(sample_info_file)
    sample_row = sample_df[sample_df['sample_id'] == sample_id]
    
    if sample_row.empty:
        raise ValueError(f"Sample {sample_id} not found in sample_info.xlsx")
    
    return sample_row.iloc[0]['sample_type']


def main():
    parser = argparse.ArgumentParser(description="UMI saturation analysis for GEX and guide libraries")
    parser.add_argument("--sample-id", required=True, help="Full sample ID (format: pool:sample)")
    parser.add_argument("--config", required=True, help="Config YAML file")
    parser.add_argument("--read_stats", required=True, help="Read statistics TSV file")
    parser.add_argument("--output_tsv", required=True, help="Output TSV file path")
    parser.add_argument("--output_plot", required=True, help="Output plot file path")
    parser.add_argument("--cell_calling_dir", required=True, help="Cell calling results directory")
    parser.add_argument("--sample_type", help="Sample type (gex or guide). If not provided, will be auto-detected.")
    parser.add_argument("--source", required=True, choices=["main", "undetermined", "all"],
                        help="Data source (main, undetermined, or all)")
    parser.add_argument("--processing", required=True, choices=["raw", "trimmed", "recovered", "merged"],
                        help="Processing state (raw, trimmed, recovered, or merged)")
    parser.add_argument("--scratch-dir", required=True, help="Scratch directory base path")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing")
    
    args = parser.parse_args()
    
    # Create output directories
    output_tsv_path = Path(args.output_tsv)
    output_plot_path = Path(args.output_plot)
    output_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    output_plot_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)
    
    # Determine sample type
    if args.sample_type:
        sample_type = args.sample_type
    else:
        sample_type = determine_sample_type(config, args.sample_id)
    
    is_guide = sample_type == 'guide'
    
    print(f"Running UMI saturation analysis for {args.sample_id} (type: {sample_type})")
    
    # Load read statistics
    print(f"Loading read statistics from: {args.read_stats}")
    total_reads = load_read_statistics(args.read_stats)
    print(f"Total reads at 100% depth: {total_reads:,}")
    
    # Handle guide vs GEX samples differently
    if is_guide:
        # For guide samples, get cell barcodes from paired GEX sample
        paired_gex_sample = get_paired_sample(args.sample_id, 'guide', config['sample_info_file'])
        print(f"Found paired GEX sample: {paired_gex_sample}")
        
        # Load cell calling results from paired GEX sample
        print("Loading cell calling results from paired GEX sample...")
        cell_calling_results, cell_barcodes_dict = load_cell_calling_results(
            args.cell_calling_dir, paired_gex_sample
        )
        
        print(f"Using cell barcodes from paired GEX sample: {paired_gex_sample}")
    else:
        # For GEX samples, use their own cell calling results
        print("Loading cell calling results...")
        cell_calling_results, cell_barcodes_dict = load_cell_calling_results(args.cell_calling_dir, args.sample_id)
    
    # Get saturation points from config
    saturation_points = config['cell_calling']['saturation_points']
    print(f"Saturation points: {saturation_points}")
    
    print(f"Found cell calling methods: {list(cell_barcodes_dict.keys())}")
    for method, barcodes in cell_barcodes_dict.items():
        print(f"  {method}: {len(barcodes)} cells")
    
    # Get guide cutoffs from config if this is a guide sample
    guide_cutoffs = None
    if is_guide:
        guide_cutoffs = config['qc_analysis']['guide_cutoffs']
        print(f"Using guide cutoffs: {guide_cutoffs}")
    
    # Run saturation analysis with cell-specific metrics
    # Always calculate GMM for guide samples
    calculate_gmm = is_guide
    
    # Set up GMM plot directory if calculating GMM
    gmm_plot_dir = None
    if calculate_gmm:
        # Place GMM plots in the same parent directory as saturation plots
        # They'll be at the same level as guide_umi_saturation
        gmm_plot_dir = output_plot_path.parent
        # No need to create directory here - will be created per plot
    
    result = run_saturation_analysis(
        config, args.sample_id, saturation_points, cell_barcodes_dict, args.source, args.processing, 
        args.scratch_dir, is_guide, calculate_gmm, guide_cutoffs, gmm_plot_dir, args.threads
    )
    
    # Handle return value based on whether GMM was calculated
    if calculate_gmm:
        saturation_data, gmm_thresholds_by_depth = result
    else:
        saturation_data = result
        gmm_thresholds_by_depth = None
    
    if saturation_data:
        # Save saturation data
        saturation_df = pd.DataFrame(saturation_data)
        saturation_df.to_csv(output_tsv_path, sep='\t', index=False)
        
        # Generate saturation plots
        print("Generating saturation plots...")
        pool = get_sample_pool_from_id(config, args.sample_id)
        plot_metadata = plot_saturation_curves(saturation_data, output_plot_path, args.sample_id, config,
                                               list(cell_barcodes_dict.keys()), is_guide, pool=pool, total_reads=total_reads,
                                               gmm_thresholds_by_depth=gmm_thresholds_by_depth, guide_cutoffs=guide_cutoffs)
        
        # Save plot metadata
        if plot_metadata:
            # Adjust path to be relative to qc_report directory
            plot_metadata['path'] = str(Path(plot_metadata['path']).relative_to(output_plot_path.parent.parent))
    
    # DEPRECATED: Diversity analysis moved to cell QC pipeline
    # Diversity analysis is now done after cell calling as part of QC
    
    feature_type = 'guide' if is_guide else 'gene expression'
    print(f"{feature_type.capitalize()} UMI saturation analysis completed.")
    print(f"Results saved to: {output_tsv_path}")
    
    # Construct the actual plot path based on the logic in plot_saturation_curves
    actual_plot_dir = output_plot_path.parent
    if is_guide:
        actual_plot_path = actual_plot_dir / 'guide_umi_saturation' / 'linear' / 'plot.png'
    else:
        actual_plot_path = actual_plot_dir / 'umi_saturation' / 'linear' / 'plot.png'
    print(f"Plot saved to: {actual_plot_path}")
    
    # Mention GMM plots if they were generated
    if gmm_plot_dir and calculate_gmm:
        print(f"GMM mixture model plots saved to: {gmm_plot_dir}/gmm_*/")

if __name__ == "__main__":
    main()