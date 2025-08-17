#!/usr/bin/env python3
"""
Generate UMAP plots from preprocessed standard_analyses output.
Creates individual PNG files for dashboard integration.
Supports both full dataset UMAPs and subset-specific UMAPs.
"""

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from pathlib import Path

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.pipeline_utils import log_print


def plot_umap_set(adata, output_dir, umap_colors, umap_key='X_umap', set_name='full'):
    """Helper function to plot a set of UMAP visualizations.
    
    Args:
        adata: AnnData object with UMAP coordinates
        output_dir: Directory to save plots
        umap_colors: List of columns to color UMAPs by
        umap_key: Key in obsm containing UMAP coordinates
        set_name: Name of this UMAP set (for logging)
    """
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    log_print(f"  📊 Plotting {set_name} with {len(umap_colors)} metrics...")
    
    # Temporarily swap UMAP coordinates if using a subset
    original_umap = None
    if umap_key != 'X_umap':
        original_umap = adata.obsm['X_umap'].copy() if 'X_umap' in adata.obsm else None
        adata.obsm['X_umap'] = adata.obsm[umap_key]
    
    # Filter to only existing columns
    umap_colors = [col for col in umap_colors if col in adata.obs.columns]
    
    missing_cols = [col for col in umap_colors if col not in adata.obs.columns]
    if missing_cols:
        log_print(f"⚠️  WARNING: Missing columns for UMAP plotting: {missing_cols}")
    
    log_print(f"🗺️ Plotting UMAP with available columns: {umap_colors}")
    
    # Set scanpy settings for rasterized output
    sc.settings.set_figure_params(dpi_save=150, facecolor='white')
    
    # Generate individual plots
    for i, color in enumerate(umap_colors):
        log_print(f"    📊 Generating UMAP {i+1}/{len(umap_colors)}: {color}")
        
        # Create metric directory
        metric_dir = output_path / color
        metric_dir.mkdir(parents=True, exist_ok=True)
        
        # Create figure
        fig, ax = plt.subplots(figsize=(6, 5))
        
        # Check if this is a target gene category column
        is_target_gene_column = (color in adata.obs.columns and 
                               adata.obs[color].dtype.name == 'category' and
                               'target_genes' in adata.obs[color].cat.categories)
        
        # Check if it's a categorical column
        is_categorical = color in adata.obs.columns and (
            adata.obs[color].dtype.name == 'category' or
            adata.obs[color].dtype == 'object'
        )
        
        # Determine legend setting
        if is_target_gene_column:
            legend_setting = True
        elif is_categorical:
            # For other categorical columns, show legend if < 20 categories
            n_categories = len(adata.obs[color].unique())
            legend_setting = n_categories < 20
        else:
            legend_setting = False
        
        # Plot UMAP with appropriate settings
        plot_kwargs = {
            'ax': ax,
            'show': False,
            'legend_loc': 'on data' if is_target_gene_column else 'right margin',
            'legend_fontsize': 'small',
            'frameon': False,
            'size': 40,
            'alpha': 0.8,
            'title': color.replace('_', ' ').title()
        }
        
        # Add vmin/vmax for continuous variables (not categorical)
        if not is_categorical:
            plot_kwargs['vmin'] = 'p5'
            plot_kwargs['vmax'] = 'p95'
        
        sc.pl.umap(adata, color=color, **plot_kwargs)
        
        # Ensure points are rasterized (not vector) for smaller file size
        for child in ax.get_children():
            if hasattr(child, 'set_rasterized'):
                child.set_rasterized(True)
        
        # Save figure
        output_file = metric_dir / "plot.png"
        plt.tight_layout()
        fig.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        
        log_print(f"      ✅ Saved: {output_file}")
    
    # Restore original UMAP if we swapped it
    if original_umap is not None:
        adata.obsm['X_umap'] = original_umap
    
    log_print(f"    ✅ Generated {len(umap_colors)} plots for {set_name}")


def save_umap_plots(adata, output_dir):
    """Save all UMAP plots including full and subset visualizations.
    
    Args:
        adata: AnnData object with UMAP coordinates
        output_dir: Output directory for plots  
    """
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Set scanpy settings for rasterized output
    sc.settings.set_figure_params(dpi_save=150, facecolor='white')
    
    # Build list of general metrics (used for all_cells_scores directories)
    general_metrics = [
        # Basic metrics
        "pct_counts_mt", "pct_counts_ribo", "rep", "library", "lib_for_batch",
        "guide_umi_counts", "total_counts", "guides_per_cell", "n_guides_total",
        
        # Gene expression scores
        "G2M_score", "S_score"
    ]
    
    # Add categories from unified scoring
    if "unified_scoring_categories" in adata.uns:
        general_metrics.extend(adata.uns["unified_scoring_categories"])
    
    # Add general scores (those without subset suffixes)
    if "unified_scoring_scores" in adata.uns:
        all_scores = adata.uns["unified_scoring_scores"]
        de_subsets = adata.uns.get("de_subsets", [])
        
        # Only include scores that don't have subset suffixes
        for score in all_scores:
            is_subset_specific = any(score.endswith(f"_{subset}") for subset in de_subsets)
            if not is_subset_specific:
                general_metrics.append(score)
    
    # Skip leiden clustering results (redundant with clusters_N which are already included)
    # leiden_cols = [col for col in adata.obs.columns if col.startswith('leiden_')]
    # general_metrics.extend(leiden_cols)
    
    # Filter to existing columns
    general_metrics = [col for col in general_metrics if col in adata.obs.columns]
    
    log_print(f"📊 General metrics to plot: {len(general_metrics)} columns")
    
    # Step 1: Collect all UMAP datasets to plot
    umap_datasets = []
    
    # Add full dataset UMAP (renamed to all_cells_umap)
    umap_datasets.append(("all_cells_umap", adata, 'X_umap'))
    log_print(f"📊 Added full dataset UMAP as 'all_cells_umap'")
    
    # Add UMAP subsets from stored information
    umap_subsets = adata.uns.get('umap_subsets', [])
    if len(umap_subsets) > 0:
        log_print(f"📊 Found {len(umap_subsets)} UMAP subsets: {umap_subsets}")
        
        for subset_name in umap_subsets:
            subset_key = f'X_umap_{subset_name}'
            if subset_key not in adata.obsm:
                log_print(f"  ⚠️  Skipping subset '{subset_name}': no UMAP coordinates found")
                continue
                
            # Check which cells have coordinates (not NaN)
            mask = ~np.isnan(adata.obsm[subset_key][:, 0])
            n_cells = mask.sum()
            
            if n_cells == 0:
                log_print(f"  ⚠️  Skipping subset '{subset_name}': no valid UMAP coordinates")
                continue
            
            # Create subset adata with only cells that have coordinates
            subset_adata = adata[mask].copy()
            subset_adata.obsm['X_umap'] = adata.obsm[subset_key][mask]
            
            umap_datasets.append((subset_name, subset_adata, 'X_umap'))
            log_print(f"  ✅ Added UMAP subset '{subset_name}' with {n_cells:,} cells")
    
    # Get DE subsets from stored information
    de_subsets = adata.uns.get("de_subsets", [])
    de_subset_scores = adata.uns.get("de_subset_scores", {})
    
    log_print("=" * 60)
    log_print(f"📊 Processing {len(umap_datasets)} UMAP datasets")
    if len(de_subsets) > 0:
        log_print(f"📊 With DE subsets: {de_subsets}")
    else:
        log_print("📊 No DE subsets found")
    log_print("=" * 60)
    
    # Step 2: For each UMAP dataset, plot all DE subsets
    for umap_name, umap_adata, umap_key in umap_datasets:
        log_print(f"\n🗺️ Processing UMAP: {umap_name}")
        umap_dir = output_path / umap_name
        
        # Always plot general scores in all_cells_scores directory
        general_dir = umap_dir / "all_cells_scores"
        log_print(f"  📊 Plotting general scores in {umap_name}/all_cells_scores/")
        plot_umap_set(umap_adata, general_dir, general_metrics, umap_key, f"{umap_name}_general")
        
        # Plot each DE subset's specific scores
        for de_subset in de_subsets:
            de_dir = umap_dir / de_subset
            subset_scores = de_subset_scores.get(de_subset, [])
            
            if len(subset_scores) == 0:
                log_print(f"  ⚠️  No scores found for DE subset '{de_subset}'")
                continue
            
            # Clean score names (remove subset suffix if present)
            cleaned_scores = []
            for score in subset_scores:
                if score.endswith(f"_{de_subset}"):
                    cleaned_score = score[:-len(f"_{de_subset}")]
                else:
                    cleaned_score = score
                
                # Only include if column exists in this UMAP subset's data
                if cleaned_score in umap_adata.obs.columns:
                    cleaned_scores.append(cleaned_score)
            
            if len(cleaned_scores) > 0:
                log_print(f"  📊 Plotting DE subset '{de_subset}' in {umap_name}/{de_subset}/ with {len(cleaned_scores)} scores")
                plot_umap_set(umap_adata, de_dir, cleaned_scores, umap_key, f"{umap_name}_{de_subset}")
            else:
                log_print(f"  ⚠️  No valid scores for DE subset '{de_subset}' in UMAP '{umap_name}'")
    
    # Create completion sentinel
    log_print("=" * 60)
    complete_file = output_path.parent / f"{output_path.name}.complete"
    complete_file.touch()
    log_print(f"✅ UMAP plotting complete: {complete_file}")


def main():
    parser = argparse.ArgumentParser(description='Generate UMAP plots from preprocessed data')
    parser.add_argument('--input-file', required=True,
                        help='Input preprocessed h5ad file from standard_analyses')
    parser.add_argument('--output-dir', required=True,
                        help='Output directory for plots')
    args = parser.parse_args()
    
    # Load preprocessed data
    log_print(f"Loading preprocessed data from: {args.input_file}")
    adata = sc.read_h5ad(args.input_file)
    
    # Check for UMAP coordinates
    if 'X_umap' not in adata.obsm:
        raise ValueError("No UMAP coordinates found in data. Run standard_analyses first.")
    
    # Generate plots
    save_umap_plots(adata, args.output_dir)
    
    log_print("✅ UMAP plotting complete!")


if __name__ == "__main__":
    main()