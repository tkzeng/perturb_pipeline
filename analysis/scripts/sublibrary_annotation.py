"""
Sublibrary Annotation Script for GW_PERTURB Pipeline
===================================================

This script performs per-sublibrary annotation and quality control, including:
1. Gene annotation (from GENCODE database)
2. Basic QC metrics (mitochondrial %, counts, genes)
3. Sample mapping (Parse Bio barcode to sample/well)
4. Guide processing (filter and assign guides)
5. Cell quality annotation (no filtering - uses cells from filtered kallisto output)

Steps 6-8 (normalization, dimensionality reduction, differential expression)
are performed after merging sublibraries in a separate script.
"""

import os
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from scipy import io
import gc
import sys
import yaml
import re
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pipeline_utils import (
    log_print, 
    add_comprehensive_gene_annotations,
    add_mitochondrial_metrics,
    map_cells_to_samples_with_plate,
    filter_guides_by_reference,
    align_datasets_by_cells,
    add_guide_data,
    get_guide_gex_pairings
)
# Removed unused imports from config.py - all configuration comes from config.yaml via Snakemake




# Removed get_expected_cells_from_sample_info - now extracted inline with other sample info


def filter_and_prepare_gex(adata, sample_id, expected_cells, no_filter=False, min_umi_filter=2):
    """Filter GEX data to top 3x expected cells + random 5x background.
    
    Args:
        adata: Unfiltered AnnData object
        sample_id: Sample identifier
        expected_cells: Expected number of cells
        no_filter: If True, skip filtering and return all barcodes
        min_umi_filter: Minimum UMI count to keep a barcode (applied even if no_filter=True)
        
    Returns:
        Filtered (or unfiltered) AnnData object
    """
    # Sum layers to get total counts
    adata.X = adata.layers['mature'] + adata.layers['nascent'] + adata.layers['ambiguous']
    log_print("   Summed mature + nascent + ambiguous layers")
    log_print(f"   adata.X sparsity: {scipy.sparse.issparse(adata.X)}, type: {type(adata.X)}")
    
    # Apply minimal UMI filter (always applied, even with no_filter=True)
    total_counts = np.array(adata.X.sum(axis=1)).flatten()
    keep_mask = total_counts >= min_umi_filter
    n_removed = (~keep_mask).sum()
    
    if n_removed > 0:
        log_print(f"   Removing {n_removed:,} barcodes with < {min_umi_filter} UMIs ({n_removed/len(total_counts)*100:.1f}%)")
        adata = adata[keep_mask, :].copy()
        # Recalculate total_counts for filtered data
        total_counts = total_counts[keep_mask]
        log_print(f"   Kept {adata.n_obs:,} barcodes with >= {min_umi_filter} UMIs")
    
    if no_filter:
        log_print(f"📊 Preparing GEX sample {sample_id} (NO FILTERING)")
        log_print(f"   Keeping all {adata.n_obs:,} barcodes")
        return adata
    
    # Original filtering code
    n_top = 3 * expected_cells  # High-quality cells for analysis
    n_random_low = max(100000, 5 * expected_cells)  # Background droplets for EmptyDrops
    
    log_print(f"📊 Filtering GEX sample {sample_id}")
    log_print(f"   Expected cells: {expected_cells}")
    log_print(f"   Keeping top {n_top} cells + random {n_random_low} low-UMI barcodes")
    
    # Sort by UMI counts (total_counts already calculated above)
    sorted_indices = np.argsort(total_counts)[::-1]  # High to low
    
    # Get top cells by UMI count
    n_top_actual = min(n_top, adata.n_obs)
    top_indices = sorted_indices[:n_top_actual]
    
    # Get random sample from all remaining barcodes
    remaining_indices = sorted_indices[n_top_actual:]
    
    n_random_actual = min(n_random_low, len(remaining_indices))
    if n_random_actual > 0:
        np.random.seed(42)  # For reproducibility
        random_low_indices = np.random.choice(remaining_indices, size=n_random_actual, replace=False)
    else:
        random_low_indices = np.array([], dtype=int)
    
    # Combine indices
    combined_indices = np.concatenate([top_indices, random_low_indices])
    
    log_print(f"   Selected {len(top_indices)} top cells + {len(random_low_indices)} random remaining cells")
    log_print(f"   Total: {len(combined_indices)} cells ({len(combined_indices)/adata.n_obs:.1%} of original)")
    
    # Filter adata
    adata_filtered = adata[combined_indices, :].copy()
    
    # Clean up original adata
    del adata
    gc.collect()
    
    return adata_filtered


def filter_and_prepare_guide(adata, sample_id, gex_barcodes, no_filter=False):
    """Filter guide data to match GEX cell barcodes.
    
    Args:
        adata: Unfiltered guide AnnData object
        sample_id: Sample identifier
        gex_barcodes: List of cell barcodes from paired GEX sample
        no_filter: If True, return all barcodes (but still need to align with GEX)
        
    Returns:
        Filtered AnnData object
    """
    if no_filter:
        log_print(f"📊 Preparing guide sample {sample_id} (NO FILTERING)")
        log_print(f"   Total guide barcodes: {adata.n_obs:,}")
    else:
        log_print(f"📊 Filtering guide sample {sample_id}")
    
    log_print(f"   Using {len(gex_barcodes)} barcodes from paired GEX")
    
    # Filter to cells present in GEX sample
    common_barcodes = [bc for bc in gex_barcodes if bc in adata.obs.index]
    log_print(f"   Found {len(common_barcodes)} common barcodes between GEX and guide")
    
    # Filter adata
    adata_filtered = adata[common_barcodes, :].copy()
    log_print(f"   Filtered: {adata_filtered.shape[0]} cells x {adata_filtered.shape[1]} guides")
    
    # Clean up original adata
    del adata
    gc.collect()
    
    return adata_filtered


def save_mtx_files(adata, output_dir, sample_type="gex"):
    """Save MTX format files needed for cell calling analysis.
    
    Args:
        adata: AnnData object to save
        output_dir: Directory to save files
        sample_type: "gex" or "guide" to determine file naming
    """
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    if sample_type == "gex":
        # Save total matrix (transpose to genes x cells for R compatibility)
        matrix_file = output_path / "cells_x_genes.total.mtx"
        io.mmwrite(matrix_file, adata.X.T)
        
        # Save individual layer matrices
        for layer_name in ['mature', 'nascent', 'ambiguous']:
            if layer_name in adata.layers:
                layer_file = output_path / f"cells_x_genes.{layer_name}.mtx"
                io.mmwrite(layer_file, adata.layers[layer_name].T)
        
        # Save barcodes and features
        barcodes_file = output_path / "cells_x_genes.barcodes.txt"
        features_file = output_path / "cells_x_genes.genes.txt"
        features_names_file = output_path / "cells_x_genes.genes.names.txt"
        
    else:  # guide
        # Save matrix
        matrix_file = output_path / "cells_x_features.mtx"
        io.mmwrite(matrix_file, adata.X.T)
        
        # Save barcodes and features
        barcodes_file = output_path / "cells_x_features.barcodes.txt"
        features_file = output_path / "cells_x_features.genes.txt"
    
    # Save barcode and feature files
    pd.DataFrame(adata.obs.index).to_csv(barcodes_file, sep='\t', header=False, index=False)
    pd.DataFrame(adata.var.index).to_csv(features_file, sep='\t', header=False, index=False)
    
    if sample_type == "gex":
        pd.DataFrame(adata.var.index).to_csv(features_names_file, sep='\t', header=False, index=False)
    
    log_print(f"   Saved MTX files to {output_path}")


def generate_qc_plots(adata, output_pdf, prefix=""):
    """Generate QC plots for sublibrary."""
    with PdfPages(output_pdf) as pdf:
        # Plot 1: n_genes vs n_counts colored by mito%
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        scatter = ax.scatter(adata.obs['n_genes'], adata.obs['total_counts'], 
                           c=adata.obs['pct_counts_mt'], cmap='viridis', 
                           alpha=0.5, s=1)
        ax.set_xlabel('Number of genes')
        ax.set_ylabel('Total counts')
        ax.set_title(f'{prefix} - Genes vs Counts (colored by mito%)')
        plt.colorbar(scatter, ax=ax, label='% mitochondrial')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # Plot 2: Distribution plots
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        
        # Total counts distribution
        axes[0,0].hist(adata.obs['total_counts'], bins=50, edgecolor='black')
        axes[0,0].set_xlabel('Total counts')
        axes[0,0].set_ylabel('Number of cells')
        axes[0,0].set_title('Total counts per cell')
        
        # Number of genes distribution
        axes[0,1].hist(adata.obs['n_genes'], bins=50, edgecolor='black')
        axes[0,1].set_xlabel('Number of genes')
        axes[0,1].set_ylabel('Number of cells')
        axes[0,1].set_title('Genes per cell')
        
        # Mitochondrial % distribution
        axes[1,0].hist(adata.obs['pct_counts_mt'], bins=50, edgecolor='black')
        axes[1,0].set_xlabel('% mitochondrial')
        axes[1,0].set_ylabel('Number of cells')
        axes[1,0].set_title('Mitochondrial content')
        
        # Guides per cell distribution
        if 'guides_per_cell' in adata.obs:
            axes[1,1].hist(adata.obs['guides_per_cell'], bins=range(0, int(adata.obs['guides_per_cell'].max()) + 2), 
                         edgecolor='black')
            axes[1,1].set_xlabel('Guides per cell')
            axes[1,1].set_ylabel('Number of cells')
            axes[1,1].set_title('Guide distribution')
        else:
            axes[1,1].text(0.5, 0.5, 'No guide data', ha='center', va='center')
            axes[1,1].set_xlim(0, 1)
            axes[1,1].set_ylim(0, 1)
        
        plt.suptitle(f'{prefix} - QC Distributions')
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()


def add_guide_statistics(adata, cutoff=5):
    """Add guide statistics to the combined data.
    
    Args:
        adata: Combined AnnData with guide_counts in obsm
        cutoff: UMI cutoff for guide assignment
    """
    if "guide_counts" not in adata.obsm:
        log_print("⚠️  No guide_counts found in obsm, skipping guide statistics")
        return
    
    guide_matrix = adata.obsm["guide_counts"]
    
    # Calculate guide statistics using sparse operations
    log_print("📊 Calculating guide statistics...")
    
    # For sparse matrices, use efficient operations
    if scipy.sparse.issparse(guide_matrix):
        log_print(f"   Processing sparse guide matrix...")
        # Count guides per cell (guides >= cutoff)
        guide_binary = (guide_matrix >= cutoff).astype(np.int8)
        log_print(f"   guide_binary sparsity: {scipy.sparse.issparse(guide_binary)}, type: {type(guide_binary)}")
        adata.obs["guides_per_cell"] = np.array(guide_binary.sum(axis=1)).flatten()
        
        # Sum UMI counts per cell
        adata.obs["guide_umi_counts"] = np.array(guide_matrix.sum(axis=1)).flatten()
        
        # Store binary guide assignment matrix as sparse
        adata.obsm["guide_assignment"] = guide_binary
    else:
        # Handle dense matrices (though this should rarely happen)
        adata.obs["guides_per_cell"] = (guide_matrix >= cutoff).sum(axis=1)
        adata.obs["guide_umi_counts"] = guide_matrix.sum(axis=1)
        adata.obsm["guide_assignment"] = scipy.sparse.csr_matrix(
            (guide_matrix >= cutoff).astype(np.int8)
        )
    
    log_print(f"✅ Guide statistics calculated: {adata.shape[0]} cells")
    log_print(f"   Mean guides per cell: {adata.obs['guides_per_cell'].mean():.2f}")
    log_print(f"   % cells with guides: {(adata.obs['guides_per_cell'] > 0).mean() * 100:.1f}%")


def main():
    """Main pipeline execution for per-sublibrary annotation."""
    # Parse arguments
    parser = argparse.ArgumentParser(description="Filter and annotate individual sublibraries")
    parser.add_argument("--gex-kb-dir", required=True, help="Path to GEX kallisto output directory")
    parser.add_argument("--guide-kb-dir", required=True, help="Path to guide kallisto output directory")
    parser.add_argument("--output-dir", required=True, help="Output directory for filtered and annotated files")
    # parser.add_argument("--qc-report", required=True, help="Path for QC report PDF")  # TODO: Implement in future
    parser.add_argument("--cutoff", type=int, default=5,
                        help="Cutoff for guide assignment")
    parser.add_argument("--config", required=True, 
                        help="Path to config.yaml file")
    parser.add_argument("--sample-id", required=True,
                        help="Sample ID (GEX sample ID in format pool:sample)")
    parser.add_argument("--guide-sample-id", required=True,
                        help="Guide sample ID")
    parser.add_argument("--skip-filtering", action="store_true",
                        help="Skip barcode filtering, keep all barcodes (recommended for cell calling)")
    parser.add_argument("--source", required=True, choices=["main", "undetermined", "all"],
                        help="Data source (main, undetermined, or all)")
    parser.add_argument("--processing", required=True, choices=["raw", "trimmed", "recovered", "merged"],
                        help="Processing state (raw, trimmed, recovered, or merged)")
    parser.add_argument("--species", required=True, choices=["human", "mouse"], help="Species for genome/MT/cell-cycle annotation (human or mouse)")

    args = parser.parse_args()

    # Start pipeline
    log_print("\n" + "=" * 80)
    log_print("🚀 STARTING SUBLIBRARY FILTERING AND ANNOTATION PIPELINE")
    
    # Extract pool from sample_id
    if ':' not in args.sample_id:
        raise ValueError(f"Invalid sample_id format: {args.sample_id}. Expected 'pool:sample'")
    pool = args.sample_id.split(':')[0]
    
    log_print(f"📊 Processing: Pool {pool}, GEX {args.sample_id}, Guide {args.guide_sample_id}")
    log_print("=" * 80)

    # Load config
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    # Extract all config values at the boundary
    log_print("📋 Extracting configuration parameters...")
    
    # Sample info file
    sample_info_file = config['sample_info_file']
    
    # Cell filtering parameters
    min_umi_filter = config['sublibrary_filtering']['min_umi_filter']
    
    # Gene annotation reference files
    ribosomal_genes_path = config['input_paths']['ribosomal_genes']
    gene_database_path = config['input_paths']['gene_database']
    cell_cycle_genes_path = config['input_paths']['cell_cycle_genes']
    
    # Barcode file for validation
    parsebio_barcodes_path = config['input_paths']['parsebio_barcodes']
    
    # Plate mapping file
    plate_maps_file = config.get('plate_maps_file')
    
    # Load barcode file once
    parsebio_barcodes = pd.read_csv(parsebio_barcodes_path)
    
    # Load sample info once and extract all needed information
    log_print("📋 Loading sample information...")
    if not os.path.exists(sample_info_file):
        raise FileNotFoundError(f"Sample info file not found: {sample_info_file}")
    
    sample_df = pd.read_excel(sample_info_file)
    sample_row = sample_df[sample_df['sample_id'] == args.sample_id]
    
    if sample_row.empty:
        raise ValueError(f"Sample {args.sample_id} not found in {sample_info_file}")
    
    # Extract expected cells
    if 'expected_cells' not in sample_df.columns:
        raise ValueError(f"'expected_cells' column not found in {sample_info_file}")
    expected_cells = sample_row.iloc[0]['expected_cells']
    if pd.isna(expected_cells):
        raise ValueError(f"Expected cells value is missing for sample {args.sample_id}")
    expected_cells = int(expected_cells)
    
    # Extract plate mapping info (will be used later)
    plate_name = sample_row.iloc[0].get('sample_to_well_mapping')
    
    log_print(f"   Sample: {args.sample_id}")
    log_print(f"   Expected cells: {expected_cells}")
    log_print(f"   Plate mapping: {plate_name if not pd.isna(plate_name) else 'None'}")
    
    # Run the pipeline
    try:
        # 1. LOAD AND FILTER GEX DATA
        log_print("\n📥 1. LOADING AND FILTERING GEX DATA...")
        gex_h5ad_path = Path(args.gex_kb_dir) / "counts_unfiltered_modified" / "adata.h5ad"
        log_print(f"   Loading from: {gex_h5ad_path}")
        
        adata_gex = sc.read_h5ad(gex_h5ad_path)
        log_print(f"   Loaded {adata_gex.shape[0]} cells x {adata_gex.shape[1]} genes")
        log_print(f"   GEX layers sparsity - mature: {scipy.sparse.issparse(adata_gex.layers['mature'])}, nascent: {scipy.sparse.issparse(adata_gex.layers['nascent'])}, ambiguous: {scipy.sparse.issparse(adata_gex.layers['ambiguous'])}")
        
        # Filter using the expected cells we extracted above
        adata_gex = filter_and_prepare_gex(adata_gex, args.sample_id, expected_cells, args.skip_filtering, min_umi_filter)
        gc.collect()  # Clean up after filtering
        
        # 2. LOAD AND FILTER GUIDE DATA
        log_print("\n📥 2. LOADING AND FILTERING GUIDE DATA...")
        guide_h5ad_path = Path(args.guide_kb_dir) / "counts_unfiltered_modified" / "adata.h5ad"
        log_print(f"   Loading from: {guide_h5ad_path}")
        
        adata_guide = sc.read_h5ad(guide_h5ad_path)
        log_print(f"   Loaded {adata_guide.shape[0]} cells x {adata_guide.shape[1]} guides")
        log_print(f"   Guide X sparsity: {scipy.sparse.issparse(adata_guide.X)}, type: {type(adata_guide.X)}")
        
        # Filter guide to match GEX barcodes
        adata_guide = filter_and_prepare_guide(adata_guide, args.guide_sample_id, adata_gex.obs.index.tolist(), args.skip_filtering)
        gc.collect()  # Clean up after filtering
        
        # 3. COMBINE GEX AND GUIDE DATA
        log_print("\n🔗 3. COMBINING GEX AND GUIDE DATA...")
        adata = add_guide_data(adata_gex, adata_guide)
        
        # Add guide statistics
        add_guide_statistics(adata, cutoff=args.cutoff)
        
        # Clean up separate objects after combining
        del adata_gex, adata_guide
        gc.collect()

        # 4. GENE ANNOTATION
        log_print("\n🏷️  4. ADDING GENE ANNOTATIONS...")
        add_comprehensive_gene_annotations(adata, ribosomal_genes_path, gene_database_path, cell_cycle_genes_path, species = args.species)

        # 5. QC METRICS
        log_print("\n📊 5. CALCULATING QC METRICS...")
        add_mitochondrial_metrics(adata)
        
        # 6. SAMPLE MAPPING AND ANNOTATION
        log_print("\n🏷️  6. ADDING SAMPLE ANNOTATIONS...")
        # Add pool and sample metadata
        adata.obs['pool'] = pool
        adata.obs['sample_id'] = args.sample_id
        
        # Use the plate_name we already extracted earlier
        if pd.isna(plate_name):
            log_print("⚠️  No plate mapping specified, skipping cell-to-sample mapping")
        else:
            # Map cells to samples based on Parse Bio barcodes
            log_print(f"   Using plate mapping: {plate_name}")
            if not plate_maps_file:
                raise ValueError("plate_maps_file must be specified in config for plate mapping")
            map_cells_to_samples_with_plate(adata, plate_name, plate_maps_file, parsebio_barcodes)
        
        # Process guide assignments (currently a stub function)
        # filter_guides_by_reference(adata)  # Skip - not needed for this pipeline

        # 7. SAVE OUTPUTS
        log_print("\n💾 7. SAVING OUTPUTS...")
        
        # Create output directory
        output_path = Path(args.output_dir)
        output_path.mkdir(exist_ok=True, parents=True)
        
        # Save filtered and annotated h5ad
        h5ad_path = output_path / "adata.h5ad"
        adata.write_h5ad(h5ad_path)
        log_print(f"   Saved annotated h5ad: {h5ad_path}")
        
        # Save MTX files for cell calling
        save_mtx_files(adata, output_path, sample_type="gex")
        
        # Generate QC plots - COMMENTED OUT: Should happen after cell calling
        # generate_qc_plots(adata, args.qc_report, prefix=f"pool{args.pool}_sample{args.sample_id}")
        # log_print(f"   Saved QC report: {args.qc_report}")

        log_print("\n✅ SUBLIBRARY FILTERING AND ANNOTATION COMPLETED SUCCESSFULLY!")

    except Exception as e:
        log_print(f"\n❌ ERROR: Pipeline failed with error: {str(e)}")
        raise


if __name__ == "__main__":
    main()