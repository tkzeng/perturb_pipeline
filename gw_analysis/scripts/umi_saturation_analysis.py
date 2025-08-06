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
# Import shared guide utility functions
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.pipeline_utils import calculate_guides_per_cell, get_paired_sample, calculate_guide_metrics_for_cells


def get_sample_pool_from_id(config, sample_id):
    """Get the pool name for a sample ID."""
    sample_info_file = config.get('sample_info_file', 'sample_info.xlsx')
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
        # Create pandas Series for vectorized barcode lookup
        import pandas as pd
        barcode_series = pd.Series(range(adata.n_obs), index=adata.obs_names)
        
        for method, cell_barcode_list in cell_barcodes_dict.items():
            # Find indices of called cells using vectorized operations
            cell_indices = barcode_series.reindex(cell_barcode_list).dropna().astype(int).values
            
            # For guide samples, use unified function for ALL metrics
            if is_guide and guide_cutoffs is not None:
                # Use unified function to calculate metrics for ALL cells
                guide_metrics = calculate_guide_metrics_for_cells(
                    guide_adata=adata,
                    cell_barcodes=cell_barcode_list,
                    guide_cutoffs=guide_cutoffs,
                    guide_to_genes=None  # Not needed for basic metrics
                )
                
                # Basic cell stats
                stats[f'{method}_n_cells'] = guide_metrics['n_cells_total']
                stats[f'{method}_mean_umis_per_cell'] = guide_metrics['guide_umis_per_cell']
                stats[f'{method}_median_umis_per_cell'] = guide_metrics['guide_umis_per_cell']  # Using mean as approximation
                
                # For mean guides per cell, use cutoff=1 as proxy for "features per cell"
                stats[f'{method}_mean_{feature_type}_per_cell'] = guide_metrics['guides_per_cell_cutoff1']
                stats[f'{method}_median_{feature_type}_per_cell'] = guide_metrics['guides_per_cell_cutoff1']
                
                # Total UMIs in cells (need to calculate from found cells only)
                if len(cell_indices) > 0:
                    cell_umis = umis_per_barcode[cell_indices]
                    stats[f'{method}_total_umis_in_cells'] = np.sum(cell_umis)
                    stats[f'{method}_fraction_umis_in_cells'] = np.sum(cell_umis) / total_umis if total_umis > 0 else 0
                else:
                    stats[f'{method}_total_umis_in_cells'] = 0
                    stats[f'{method}_fraction_umis_in_cells'] = 0
                
                # Guide-specific metrics for each cutoff
                for cutoff in guide_cutoffs:
                    stats[f'{method}_mean_guides_per_cell_cutoff{cutoff}'] = guide_metrics[f'guides_per_cell_cutoff{cutoff}']
                    stats[f'{method}_median_guides_per_cell_cutoff{cutoff}'] = guide_metrics[f'guides_per_cell_cutoff{cutoff}']
                    stats[f'{method}_fraction_cells_with_guides_cutoff{cutoff}'] = guide_metrics[f'fraction_cells_with_guides_cutoff{cutoff}']
                    
            elif len(cell_indices) > 0:
                # For GEX samples, use original approach
                cell_umis = umis_per_barcode[cell_indices]
                cell_features = features_per_barcode[cell_indices]
                
                stats[f'{method}_n_cells'] = len(cell_barcode_list)
                stats[f'{method}_mean_umis_per_cell'] = np.mean(cell_umis)
                stats[f'{method}_median_umis_per_cell'] = np.median(cell_umis)
                stats[f'{method}_mean_{feature_type}_per_cell'] = np.mean(cell_features)
                stats[f'{method}_median_{feature_type}_per_cell'] = np.median(cell_features)
                stats[f'{method}_total_umis_in_cells'] = np.sum(cell_umis)
                stats[f'{method}_fraction_umis_in_cells'] = np.sum(cell_umis) / total_umis if total_umis > 0 else 0
            else:
                stats[f'{method}_n_cells'] = 0
                stats[f'{method}_mean_umis_per_cell'] = 0
                stats[f'{method}_median_umis_per_cell'] = 0
                stats[f'{method}_mean_{feature_type}_per_cell'] = 0
                stats[f'{method}_median_{feature_type}_per_cell'] = 0
                stats[f'{method}_total_umis_in_cells'] = 0
                stats[f'{method}_fraction_umis_in_cells'] = 0
                
                # Add zero values for guide metrics
                if is_guide:
                    if guide_cutoffs is None:
                        raise ValueError("Guide cutoffs must be provided for guide samples")
                    for cutoff in guide_cutoffs:
                        stats[f'{method}_mean_guides_per_cell_cutoff{cutoff}'] = 0
                        stats[f'{method}_median_guides_per_cell_cutoff{cutoff}'] = 0
                        stats[f'{method}_fraction_cells_with_guides_cutoff{cutoff}'] = 0
    
    # Clean up intermediate arrays
    del umis_per_barcode, features_per_barcode
    if cell_barcodes_dict and 'barcode_series' in locals():
        del barcode_series
    
    return stats


def run_saturation_analysis(config, sample_id, saturation_points, cell_barcodes_dict, source, processing, scratch_dir, is_guide=False, guide_cutoffs=None):
    """Run UMI saturation analysis by loading pre-generated subsampled matrices."""
    saturation_data = []
    
    # Process each saturation fraction
    for fraction in saturation_points:
        print(f"Processing fraction {fraction}...")
        
        if fraction == 1.0:
            # Use existing full matrix (unfiltered)
            try:
                if is_guide:
                    kb_result_dir = Path(f"{scratch_dir}/{sample_id}/kb_guide_{source}_{processing}")
                else:
                    kb_result_dir = Path(f"{scratch_dir}/{sample_id}/kb_all_{source}_{processing}")
                
                adata = load_adata_all_barcodes(kb_result_dir, is_guide)
                print(f"  Using existing full matrix from {kb_result_dir}")
            except Exception as e:
                raise RuntimeError(f"Failed to load existing matrix for fraction 1.0: {e}. Cannot proceed with saturation analysis.") from e
        else:
            # Load pre-generated subsampled matrix (unfiltered)
            try:
                # Temp is always under scratch/tmp
                if is_guide:
                    subsample_dir = f"{scratch_dir}/tmp/umi_sat_{sample_id}-{fraction}/kb_guide"
                else:
                    subsample_dir = f"{scratch_dir}/tmp/umi_sat_{sample_id}-{fraction}/kb_all"
                
                adata = load_adata_all_barcodes(Path(subsample_dir), is_guide)
                print(f"  Loaded subsampled matrix from {subsample_dir}")
            except Exception as e:
                raise RuntimeError(f"Failed to load subsampled matrix for fraction {fraction}: {e}. Make sure Snakemake has generated the subsampled data.") from e
        
        # Calculate statistics including cell-specific metrics
        stats = calculate_umi_stats(adata, cell_barcodes_dict, is_guide, guide_cutoffs)
        
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
    
    return saturation_data

def plot_saturation_curves(saturation_data, output_plot_path, sample_id, cell_methods=None, is_guide=False, guide_cutoffs=None, pool=None, total_reads=None):
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
            # For guide samples, add extra row for guide-per-cell metrics
            # Now we need 5 rows: 2 basic, 1 for UMI/genes, 2 for cutoff 1&2
            fig, axes = plt.subplots(5, 2, figsize=(15, 30))
        else:
            fig, axes = plt.subplots(3, 2, figsize=(15, 18))
    else:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # UMI saturation curve
    axes[0, 0].plot(df['subsample_fraction'], df['total_umis'], 'bo-')
    axes[0, 0].set_xlabel('Subsample Fraction')
    axes[0, 0].set_ylabel('Total UMIs')
    title = f'{feature_type_cap} UMI Saturation Curve'
    if total_reads:
        title += f'\nTotal Reads at 100%: {total_reads:,}'
    axes[0, 0].set_title(title)
    axes[0, 0].grid(True, alpha=0.3)
    
    # Active barcodes saturation
    axes[0, 1].plot(df['subsample_fraction'], df['barcodes_with_umis'], 'ro-')
    axes[0, 1].set_xlabel('Subsample Fraction')
    axes[0, 1].set_ylabel('Barcodes with UMIs')
    axes[0, 1].set_title('Active Barcodes Saturation')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Mean UMIs per barcode
    axes[1, 0].plot(df['subsample_fraction'], df['mean_umis_per_barcode'], 'go-')
    axes[1, 0].set_xlabel('Subsample Fraction')
    axes[1, 0].set_ylabel('Mean UMIs per Barcode')
    axes[1, 0].set_title('Mean UMIs per Barcode')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Mean features per barcode
    feature_col = f'mean_{feature_type}_per_barcode'
    axes[1, 1].plot(df['subsample_fraction'], df[feature_col], 'mo-')
    axes[1, 1].set_xlabel('Subsample Fraction')
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
                axes[2, 0].plot(df['subsample_fraction'], df[mean_umis_col], 
                               'o-', color=colors[i], label=method)
        
        axes[2, 0].set_xlabel('Subsample Fraction')
        axes[2, 0].set_ylabel('Mean UMIs per Cell')
        axes[2, 0].set_title(f'Mean {feature_type_cap} UMIs per Cell by Method')
        axes[2, 0].grid(True, alpha=0.3)
        axes[2, 0].legend()
        
        # For non-guide samples, show average genes per cell
        if not is_guide:
            for i, method in enumerate(cell_methods):
                mean_genes_col = f'{method}_mean_{feature_type}_per_cell'
                if mean_genes_col in df.columns:
                    axes[2, 1].plot(df['subsample_fraction'], df[mean_genes_col], 
                                   'o-', color=colors[i], label=method)
            
            axes[2, 1].set_xlabel('Subsample Fraction')
            axes[2, 1].set_ylabel(f'Mean {feature_type_cap}s per Cell')
            axes[2, 1].set_title(f'Mean {feature_type_cap}s per Cell by Method')
            axes[2, 1].grid(True, alpha=0.3)
            axes[2, 1].legend()
        else:
            # For guide samples, we'll use this space for cutoff 1
            for i, method in enumerate(cell_methods):
                mean_guides_col = f'{method}_mean_guides_per_cell_cutoff1'
                if mean_guides_col in df.columns:
                    axes[2, 1].plot(df['subsample_fraction'], df[mean_guides_col], 
                                   'o-', color=colors[i], label=method)
            
            axes[2, 1].set_xlabel('Subsample Fraction')
            axes[2, 1].set_ylabel('Mean Guides per Cell')
            axes[2, 1].set_title('Mean Guides per Cell by Method (Cutoff = 1 UMI)')
            axes[2, 1].grid(True, alpha=0.3)
            axes[2, 1].legend()
        
        # Add guide-specific plots if this is a guide sample
        if is_guide and cell_methods:
            # Add Mean Guides per Cell for cutoff 2
            for i, method in enumerate(cell_methods):
                mean_guides_col = f'{method}_mean_guides_per_cell_cutoff2'
                if mean_guides_col in df.columns:
                    axes[3, 0].plot(df['subsample_fraction'], df[mean_guides_col], 
                                   'o-', color=colors[i], label=method)
            
            axes[3, 0].set_xlabel('Subsample Fraction')
            axes[3, 0].set_ylabel('Mean Guides per Cell')
            axes[3, 0].set_title('Mean Guides per Cell by Method (Cutoff = 2 UMIs)')
            axes[3, 0].grid(True, alpha=0.3)
            axes[3, 0].legend()
            
            # Plot guides per cell for different cutoffs
            cutoff_colors = ['purple', 'orange', 'brown', 'pink', 'gray']
            
            # Mean guides per cell at different cutoffs
            if guide_cutoffs is None:
                raise ValueError("Guide cutoffs must be provided for guide plots")
            
            # Use EmptyDrops_FDR001 method specifically for single-method plots
            single_method = 'EmptyDrops_FDR001'
            if single_method not in cell_methods:
                raise ValueError(f"Required method {single_method} not found in cell calling results. Available methods: {cell_methods}")
            
            for j, cutoff in enumerate(guide_cutoffs):
                col_name = f'{single_method}_mean_guides_per_cell_cutoff{cutoff}'
                if col_name in df.columns:
                    axes[4, 0].plot(df['subsample_fraction'], df[col_name], 
                                   'o-', color=cutoff_colors[j % len(cutoff_colors)], label=f'Cutoff {cutoff}')
            
            axes[4, 0].set_xlabel('Subsample Fraction')
            axes[4, 0].set_ylabel('Mean Guides per Cell')
            axes[4, 0].set_title(f'Guide Detection Saturation by Cutoff ({single_method})')
            axes[4, 0].grid(True, alpha=0.3)
            axes[4, 0].legend()
            
            # Fraction of cells with guides at different cutoffs
            method = 'EmptyDrops_FDR001'
            if method not in cell_methods:
                raise ValueError(f"Required method {method} not found in cell calling results. Available methods: {cell_methods}")
            
            for j, cutoff in enumerate(guide_cutoffs):
                col_name = f'{method}_fraction_cells_with_guides_cutoff{cutoff}'
                if col_name in df.columns:
                    axes[4, 1].plot(df['subsample_fraction'], df[col_name], 
                                   'o-', color=cutoff_colors[j % len(cutoff_colors)], label=f'Cutoff {cutoff}')
            
            axes[4, 1].set_xlabel('Subsample Fraction')
            axes[4, 1].set_ylabel('Fraction of Cells with Guides')
            axes[4, 1].set_title(f'Cell Coverage Saturation by Cutoff ({method})')
            axes[4, 1].grid(True, alpha=0.3)
            axes[4, 1].legend()
            
            # Hide the unused subplot for guide samples
            axes[3, 1].set_visible(False)
    
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
    
    plt.savefig(clean_path, dpi=300, bbox_inches='tight')
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
#                 dpi=300, bbox_inches='tight')
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
    sample_info_file = config.get('sample_info_file', 'sample_info.xlsx')
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
        paired_gex_sample = get_paired_sample(args.sample_id, 'guide', config.get('sample_info_file', 'sample_info.xlsx'))
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
    saturation_points = config.get('cell_calling', {}).get('defaults', {}).get('saturation_points', [0.1, 0.25, 0.5, 0.75, 1.0])
    print(f"Saturation points: {saturation_points}")
    
    print(f"Found cell calling methods: {list(cell_barcodes_dict.keys())}")
    for method, barcodes in cell_barcodes_dict.items():
        print(f"  {method}: {len(barcodes)} cells")
    
    # Get guide cutoffs from config if this is a guide sample
    guide_cutoffs = None
    if is_guide:
        guide_cutoffs = config.get('qc_analysis', {}).get('guide_cutoffs', [1, 2, 3, 4, 5])
        print(f"Using guide cutoffs: {guide_cutoffs}")
    
    # Run saturation analysis with cell-specific metrics
    saturation_data = run_saturation_analysis(
        config, args.sample_id, saturation_points, cell_barcodes_dict, args.source, args.processing, 
        args.scratch_dir, is_guide, guide_cutoffs
    )
    
    if saturation_data:
        # Save saturation data
        saturation_df = pd.DataFrame(saturation_data)
        saturation_df.to_csv(output_tsv_path, sep='\t', index=False)
        
        # Generate saturation plots
        print("Generating saturation plots...")
        pool = get_sample_pool_from_id(config, args.sample_id)
        plot_metadata = plot_saturation_curves(saturation_data, output_plot_path, args.sample_id, 
                                               list(cell_barcodes_dict.keys()), is_guide, guide_cutoffs, pool=pool, total_reads=total_reads)
        
        # Save plot metadata
        if plot_metadata:
            # Adjust path to be relative to qc_report directory
            plot_metadata['path'] = str(Path(plot_metadata['path']).relative_to(output_plot_path.parent.parent))
    
    # DEPRECATED: Diversity analysis moved to cell QC pipeline
    # Diversity analysis is now done after cell calling as part of QC
    
    feature_type = 'guide' if is_guide else 'gene expression'
    print(f"{feature_type.capitalize()} UMI saturation analysis completed.")
    print(f"Results saved to: {output_tsv_path}")
    print(f"Plot saved to: {output_plot_path}")

if __name__ == "__main__":
    main()