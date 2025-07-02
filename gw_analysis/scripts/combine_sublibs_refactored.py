import anndata as ad
import argparse
import numpy as np
import scanpy as sc
import pandas as pd
import gc
import time
from pathlib import Path
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config import *
from scripts.pipeline_utils import *

def get_pools_gex_libraries(pool_names, sample_info_file):
    """Get all GEX libraries for a list of pools from sample_info.xlsx
    
    Args:
        pool_names: List of pool names to process (empty list means all pools)
        sample_info_file: Path to sample_info.xlsx file
        
    Returns:
        Dictionary of {sample_id: sample_id} for all GEX samples in the specified pools
    """
    # Load sample information
    sample_df = pd.read_excel(sample_info_file)
    
    # If empty list, use all pools
    if not pool_names:
        pool_names = sample_df['pool'].unique().tolist()
        log_print(f"No pools specified, using all pools: {pool_names}")
    
    # Filter for GEX samples in specified pools
    pool_gex = sample_df[
        (sample_df['pool'].isin(pool_names)) & 
        (sample_df['sample_type'] == 'gex')
    ]
    
    if pool_gex.empty:
        raise ValueError(f"No GEX samples found for pools {pool_names}")
    
    # Create library dictionary with standardized IDs for GEX/guide pairing
    # Extract standardized ID from sample_id (remove gex/guide prefix)
    gex_libs = {}
    for _, row in pool_gex.iterrows():
        sample_id = row['sample_id']
        # Extract standardized ID: gex_1 -> 1, gex_batch_2 -> batch_2
        if sample_id.startswith('gex_'):
            standardized_id = sample_id[4:]  # Remove 'gex_' prefix
        elif sample_id.startswith('gex'):
            standardized_id = sample_id[3:]  # Remove 'gex' prefix
        else:
            standardized_id = sample_id  # Fallback
        gex_libs[sample_id] = standardized_id
    
    log_print(f"Found {len(gex_libs)} GEX libraries across pools {pool_names}: {list(gex_libs.keys())}")
    return gex_libs

def create_unified_gex_file(pool_names, sample_info_file, output_file="combined_gex.h5ad"):
    """Create single combined_gex.h5ad with total counts in X and nascent data in layers (if available)
    
    Args:
        pool_names: List of pool names to process (empty list means all pools)
        sample_info_file: Path to sample_info.xlsx file
        output_file: Output file path
    """
    
    # Note: Snakemake handles dependency tracking, so we always regenerate when called
    if Path(output_file).exists():
        log_print(f"🔄 {output_file} exists - regenerating as requested by Snakemake...")
    
    log_print(f"📦 Creating unified GEX file for pools {pool_names}...")
    
    # Get GEX libraries for specified pools
    gex_libs = get_pools_gex_libraries(pool_names, sample_info_file)
    
    # Process libraries for total counts (standard processing)
    log_print("📊 === Processing libraries for total counts ===")
    adata_total_list = []
    gex_lib_names = list(gex_libs.keys())
    
    for lib in gex_lib_names[::-1]:  # Process in reverse order
        adata = process_single_library(lib, nascent=False)
        adata_total_list.append(adata)
    
    if not adata_total_list:
        raise RuntimeError("No libraries processed successfully for total counts")
    
    # Concatenate total count libraries
    log_print("🔗 Concatenating total count libraries...")
    adata_total = ad.concat(adata_total_list)
    log_print(f"📊 Total counts shape: {adata_total.shape}")
    
    # Log biological sample distribution if available
    if 'biological_sample' in adata_total.obs.columns:
        sample_counts = adata_total.obs['biological_sample'].value_counts()
        log_print(f"📊 Combined biological samples: {len(sample_counts)} unique samples")
        for sample, count in sample_counts.head(10).items():
            log_print(f"  {sample}: {count} cells")
        if len(sample_counts) > 10:
            log_print(f"  ... and {len(sample_counts) - 10} more samples")
    else:
        log_print("📋 No biological sample annotations found")
    
    # Free memory from intermediate list
    del adata_total_list
    gc.collect()
    
    # Get cells from total dataset for filtering nascent data
    cells_total = set(adata_total.obs.index)
    log_print(f"📊 Total dataset has {len(cells_total)} cells for reference")
    
    # Try to process libraries for nascent counts (thresh1 processing) - OPTIONAL
    log_print("\n🧬 === Attempting to process libraries for nascent counts ===")
    adata_nascent_list = []
    total_nascent_cells_before = 0
    total_nascent_cells_after = 0
    nascent_processing_successful = True
    
    try:
        for i, lib in enumerate(gex_lib_names[::-1]):  # Process in reverse order
            log_print(f"[{i+1}/{len(gex_lib_names)}] 🔄 Processing {lib} for nascent...")
            try:
                adata = process_single_library(lib, nascent=True)
                total_nascent_cells_before += len(adata)
                
                # Subset to cells present in total dataset BEFORE adding to list
                log_print(f"  🔍 Finding common cells with total dataset...")
                cells_nascent = set(adata.obs.index)
                cells_common = cells_total & cells_nascent
                
                if len(cells_common) > 0:
                    log_print(f"  Subsetting from {len(adata)} to {len(cells_common)} cells...", flush=True)
                    adata_subset = adata[list(cells_common)]
                    total_nascent_cells_after += len(adata_subset)
                    adata_nascent_list.append(adata_subset)
                    log_print(f"  ✅ Added to list. Running total: {total_nascent_cells_after} cells", flush=True)
                    
                    # Free memory from full dataset
                    log_print(f"  Cleaning up memory for {lib}...", flush=True)
                    del adata
                    gc.collect()
                else:
                    log_print(f"  ⚠️ Warning: No common cells found between total and nascent data for {lib}")
                    nascent_processing_successful = False
                    break
                    
            except (FileNotFoundError, KeyError) as e:
                log_print(f"  ⚠️ Warning: Could not process nascent data for {lib}: {e}")
                nascent_processing_successful = False
                break
                
    except (FileNotFoundError, KeyError, OSError) as e:
        log_print(f"⚠️ Warning: Nascent processing failed: {e}")
        log_print("  📋 Nascent processing is optional - continuing with total counts only")
        nascent_processing_successful = False
    
    # Handle nascent data if available
    if nascent_processing_successful and adata_nascent_list:
        log_print(f"\n📊 Nascent processing summary:", flush=True)
        log_print(f"  Total cells before filtering: {total_nascent_cells_before:,}", flush=True)
        log_print(f"  Total cells after filtering: {total_nascent_cells_after:,}", flush=True)
        log_print(f"  Reduction: {(1 - total_nascent_cells_after/total_nascent_cells_before)*100:.1f}%", flush=True)
        
        # Concatenate nascent count libraries
        log_print(f"\n🔄 Concatenating {len(adata_nascent_list)} nascent libraries...", flush=True)
        concat_start_time = time.time()
        adata_nascent = ad.concat(adata_nascent_list)
        concat_time = time.time() - concat_start_time
        log_print(f"✅ Concatenation completed in {concat_time:.2f} seconds", flush=True)
        log_print(f"📊 Nascent counts shape: {adata_nascent.shape}", flush=True)
        
        # Free memory from intermediate list
        log_print("🧹 Cleaning up concatenation memory...", flush=True)
        del adata_nascent_list
        gc.collect()
        
        # Add nascent layers to total dataset
        add_nascent_layers_to_total(adata_total, adata_nascent)
        
        # Free memory from nascent dataset after layer copying
        log_print("🧹 Cleaning up nascent dataset memory...", flush=True)
        del adata_nascent
        gc.collect()
        
    else:
        log_print("📊 No nascent data available - continuing with total counts only")
        # Add placeholder nascent metadata
        total_counts = np.array(adata_total.X.sum(axis=1)).flatten()
        adata_total.obs["total_nascent_counts"] = 0
        adata_total.obs["nascent_fraction"] = 0.0
    
    log_print(f"\n📊 Final unified dataset:", flush=True)
    log_print(f"  Shape: {adata_total.shape}", flush=True)
    log_print(f"  Layers: {list(adata_total.layers.keys()) if adata_total.layers else 'None'}", flush=True)
    
    log_print(f"\n💾 Writing to {output_file}...", flush=True)
    write_start = time.time()
    adata_total.write_h5ad(output_file)
    write_time = time.time() - write_start
    log_print(f"✅ File written in {write_time:.2f} seconds", flush=True)
    
    log_print(f"Unified GEX file for pools {pool_names} created successfully!")
    return True

def add_nascent_layers_to_total(adata_total, adata_nascent):
    """Add nascent layers to total dataset (helper function)"""
    
    # ⚠️  CRITICAL: TOTAL RNA IS AUTHORITATIVE FOR CELL CALLING ⚠️
    # ⚠️  KEEP ALL TOTAL CELLS, ADD NASCENT DATA WHERE AVAILABLE ⚠️
    log_print("\n🔄 Aligning nascent data to total dataset (total is authoritative)...", flush=True)
    cells_total_list = list(adata_total.obs.index)
    cells_nascent_set = set(adata_nascent.obs.index)
    common_cells_ordered = [cell for cell in cells_total_list if cell in cells_nascent_set]
    
    log_print(f"📊 Total cells (authoritative): {len(cells_total_list):,}", flush=True)
    log_print(f"📊 Nascent cells available: {adata_nascent.shape[0]:,}", flush=True)
    log_print(f"📊 Cells with nascent data: {len(common_cells_ordered):,}", flush=True)
    log_print(f"📊 Cells missing nascent data: {len(cells_total_list) - len(common_cells_ordered):,}", flush=True)
    
    if len(common_cells_ordered) == 0:
        raise RuntimeError("CRITICAL: No common cells between total and nascent data - this indicates a fundamental data integrity issue")
    
    # KEEP ALL TOTAL CELLS (no subsetting of total dataset)
    # Only subset nascent to matching cells
    log_print("🔄 Subsetting nascent data to match total cells...", flush=True)
    adata_nascent_common = adata_nascent[common_cells_ordered]
    log_print(f"✅ Total dataset unchanged: {adata_total.shape[0]} cells (authoritative)")
    log_print(f"✅ Nascent subset: {adata_nascent_common.shape[0]} cells")
    
    # Ensure genes are aligned between datasets
    genes_total = set(adata_total.var.index)
    genes_nascent = set(adata_nascent_common.var.index)
    # Use total dataset gene order as authoritative (it went through cell filtering)
    common_genes = [gene for gene in adata_total.var.index if gene in genes_nascent]
    
    log_print(f"Total genes: {len(genes_total)}")
    log_print(f"Nascent genes: {len(genes_nascent)}")
    log_print(f"Common genes: {len(common_genes)}")
    
    if len(common_genes) == 0:
        raise RuntimeError("CRITICAL: No common genes between total and nascent data - this indicates a fundamental data integrity issue")
    
    # Filter datasets to common genes
    adata_total_subset = adata_total[:, common_genes]
    adata_nascent_common = adata_nascent_common[:, common_genes]
    
    # Update adata_total in place
    adata_total._inplace_subset_var(adata_total.var.index.isin(common_genes))
    
    # Create unified dataset with nascent data in layers (in-place modification)
    log_print("\n🔄 Creating unified dataset with nascent layers...", flush=True)
    
    # Rename existing total layers to be explicit
    layers_added = 0
    if "mature" in adata_total.layers:
        adata_total.layers["mature_total"] = adata_total.layers.pop("mature")
        log_print("✅ Renamed mature → mature_total", flush=True)
        layers_added += 1
    
    if "ambiguous" in adata_total.layers:
        adata_total.layers["ambiguous_total"] = adata_total.layers.pop("ambiguous")
        log_print("✅ Renamed ambiguous → ambiguous_total", flush=True)
        layers_added += 1
    
    if "nascent" in adata_total.layers:
        adata_total.layers["nascent_total"] = adata_total.layers.pop("nascent")
        log_print("✅ Renamed nascent → nascent_total", flush=True)
        layers_added += 1
    
    # Create nascent layers for ALL total cells (with zeros for missing cells)
    import scipy.sparse as sp
    n_total_cells = adata_total.shape[0]
    n_genes = len(common_genes)
    
    # Get indices of cells that have nascent data
    cell_to_idx = {cell: i for i, cell in enumerate(adata_total.obs.index)}
    nascent_indices = [cell_to_idx[cell] for cell in adata_nascent_common.obs.index]
    
    log_print(f"📊 Creating nascent layers for {n_total_cells} total cells ({len(nascent_indices)} with data)")
    
    # Add nascent layers from nascent dataset (with zeros for missing cells)
    nascent_layer_names = ["mature", "ambiguous", "nascent"]
    for layer_name in nascent_layer_names:
        if layer_name in adata_nascent_common.layers:
            # Get nascent data and convert to coo format for efficient construction
            nascent_data = adata_nascent_common.layers[layer_name]
            if not sp.issparse(nascent_data):
                nascent_data = sp.csr_matrix(nascent_data)
            nascent_coo = nascent_data.tocoo()
            
            # Map cell indices from nascent to total dataset coordinates
            row_mapping = np.array(nascent_indices)
            new_rows = row_mapping[nascent_coo.row]
            
            # Create sparse matrix directly with correct coordinates
            full_layer = sp.csr_matrix(
                (nascent_coo.data, (new_rows, nascent_coo.col)), 
                shape=(n_total_cells, n_genes)
            )
            
            # Store in total dataset
            adata_total.layers[f"{layer_name}_nascent"] = full_layer
            log_print(f"✅ Added {layer_name}_nascent layer (zeros for {n_total_cells - len(nascent_indices)} cells)")
            layers_added += 1
    
    # Add metadata about nascent vs total counts
    log_print("🔄 Calculating nascent vs total count metrics...", flush=True)
    total_counts = np.array(adata_total.X.sum(axis=1)).flatten()
    
    # Calculate total nascent counts from nascent layers if they exist
    if "mature_nascent" in adata_total.layers:
        nascent_counts = (adata_total.layers["mature_nascent"] + 
                         adata_total.layers["ambiguous_nascent"] + 
                         adata_total.layers["nascent_nascent"]).sum(axis=1)
        nascent_counts = np.array(nascent_counts).flatten()
        
        adata_total.obs["total_nascent_counts"] = nascent_counts
        adata_total.obs["nascent_fraction"] = nascent_counts / (total_counts + 1)  # +1 to avoid division by zero
        log_print("✅ Added nascent count metadata", flush=True)
        log_print(f"  Nascent fraction range: {adata_total.obs['nascent_fraction'].min():.3f} - {adata_total.obs['nascent_fraction'].max():.3f}", flush=True)
    else:
        # No nascent layers available
        adata_total.obs["total_nascent_counts"] = 0
        adata_total.obs["nascent_fraction"] = 0.0
        log_print("✅ Added placeholder nascent metadata (no nascent data)", flush=True)

def main():
    parser = argparse.ArgumentParser(description='Create unified GEX data for specified pools')
    parser.add_argument('--pools', nargs='*', default=[], help='Pool names to process (empty means all pools)')
    parser.add_argument('--sample-info', default='../references/sample_info.xlsx', help='Path to sample_info.xlsx file')
    parser.add_argument('--output', required=True, help='Output file path')
    args = parser.parse_args()
    
    success = create_unified_gex_file(
        pool_names=args.pools,
        sample_info_file=args.sample_info,
        output_file=args.output
    )
    if success:
        log_print("Done!")
    else:
        log_print("Failed to create unified file")
        exit(1)

if __name__ == "__main__":
    main()