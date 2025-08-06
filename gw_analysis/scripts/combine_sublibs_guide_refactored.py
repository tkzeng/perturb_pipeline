import argparse
import anndata as ad
import scanpy as sc
import pandas as pd
from pathlib import Path
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# TODO: Remove dependency on config.py and use YAML config instead
from config import *
from scripts.pipeline_utils import *

# This script combines GEX and guide data into unified format

def get_pools_guide_libraries(pool_names, sample_info_file):
    """Get all guide libraries for a list of pools from sample_info.xlsx
    
    Args:
        pool_names: List of pool names to process (empty list means all pools)
        sample_info_file: Path to sample_info.xlsx file
        
    Returns:
        Dictionary of {sample_id: sample_id} for all guide samples in the specified pools
    """
    # Load sample information
    sample_df = pd.read_excel(sample_info_file)
    
    # If empty list, use all pools
    if not pool_names:
        pool_names = sample_df['pool'].unique().tolist()
        log_print(f"No pools specified, using all pools: {pool_names}")
    
    # Filter for guide samples in specified pools
    pool_guide = sample_df[
        (sample_df['pool'].isin(pool_names)) & 
        (sample_df['sample_type'] == 'guide')
    ]
    
    if pool_guide.empty:
        raise ValueError(f"No guide samples found for pools {pool_names}")
    
    # Create library dictionary with standardized IDs for GEX/guide pairing
    # Extract standardized ID from sample_id (remove gex/guide prefix)
    guide_libs = {}
    for _, row in pool_guide.iterrows():
        sample_id = row['sample_id']
        # Extract standardized ID: guide_1 -> 1, guide_batch_2 -> batch_2
        if sample_id.startswith('guide_'):
            standardized_id = sample_id[6:]  # Remove 'guide_' prefix
        elif sample_id.startswith('guide'):
            standardized_id = sample_id[5:]  # Remove 'guide' prefix  
        else:
            standardized_id = sample_id  # Fallback
        guide_libs[sample_id] = standardized_id
    
    log_print(f"Found {len(guide_libs)} guide libraries across pools {pool_names}: {list(guide_libs.keys())}")
    return guide_libs

def load_and_format_gex_data(gex_file):
    """Load and format unified GEX data.
    
    Args:
        gex_file: Path to unified GEX file
        
    Returns:
        Formatted AnnData object
        
    Raises:
        FileNotFoundError: If GEX file doesn't exist
    """
    if not Path(gex_file).exists():
        raise FileNotFoundError(f"GEX file not found: {gex_file}. Create it first using combine_sublibs_refactored.py")
    
    log_print(f"Loading unified GEX data from {gex_file}...")
    adata_gex = sc.read_h5ad(gex_file)
    log_print(f"GEX data shape: {adata_gex.shape}")
    log_print(f"Available layers: {list(adata_gex.layers.keys()) if adata_gex.layers else 'None'}")
    
    # Cell indices already standardized during individual library processing
    log_print("GEX cell indices already standardized during processing")
    
    return adata_gex

def create_and_format_guide_data(pool_names, sample_info_file, target_cells=None, output_file="combined_guide.h5ad"):
    """Create and format guide data for specified pools.
    
    Args:
        pool_names: List of pool names to process (empty list means all pools)
        sample_info_file: Path to sample_info.xlsx file
        target_cells: Set of target cell IDs to subset to (filtered GEX cells)
        output_file: Path for guide file
        
    Returns:
        Combined guide AnnData object
        
    Raises:
        RuntimeError: If no guide libraries are processed successfully
    """
    log_print(f"Processing guide libraries for pools {pool_names}...")
    
    # Get guide libraries for specified pools
    guide_libs = get_pools_guide_libraries(pool_names, sample_info_file)
    
    adata_guide = []
    
    for lib in guide_libs.keys():
        adata = process_guide_library(lib, guide_libs=guide_libs, target_cells=target_cells)
        adata_guide.append(adata)
    
    if not adata_guide:
        raise RuntimeError("No guide libraries processed successfully")
    
    log_print(f"Concatenating {len(adata_guide)} guide libraries...")
    adata_guide_combined = ad.concat(adata_guide)
    log_print(f"Combined guide data shape: {adata_guide_combined.shape}")
    
    # Cell IDs already standardized within process_guide_library() before subsetting
    
    # Write guide file
    log_print(f"Writing {output_file}...")
    adata_guide_combined.write_h5ad(output_file)
    
    # Clean up memory
    del adata_guide
    gc.collect()
    
    return adata_guide_combined

def main():
    parser = argparse.ArgumentParser(description='Create guide data file for specified pools')
    parser.add_argument('--pools', nargs='*', default=[], help='Pool names to process (empty means all pools)')
    parser.add_argument('--sample-info', default='../references/sample_info.xlsx', help='Path to sample_info.xlsx file')
    parser.add_argument('--input-gex', required=True, help='Input GEX file path (for cell filtering)')
    parser.add_argument('--output-guide', required=True, help='Output guide file path')
    parser.add_argument('--output-metadata', required=True, help='Output guide metadata file path')
    args = parser.parse_args()
    
    # Set output files from command line arguments
    input_gex = args.input_gex
    output_guide = args.output_guide
    output_metadata = args.output_metadata
    pool_names = args.pools
    sample_info_file = args.sample_info
    
    # Note: Snakemake handles dependency tracking, so we always regenerate when called
    if Path(output_guide).exists():
        log_print(f"{output_guide} exists - regenerating as requested by Snakemake...")
    
    log_print(f"=== Creating guide data file for pools {pool_names} ===")
    
    # Load GEX data only to get the filtered cell list
    log_print("Loading GEX data to get filtered cell list...")
    adata_gex = load_and_format_gex_data(input_gex)
    
    gex_cells = set(adata_gex.obs.index)
    log_print(f"Found {len(gex_cells)} cells in GEX data for filtering")
    
    # Clean up GEX data - we only needed the cell list
    del adata_gex
    gc.collect()
    
    # Create and format guide data (subset to GEX cells, write to final output)
    adata_guide = create_and_format_guide_data(
        pool_names=pool_names,
        sample_info_file=sample_info_file,
        target_cells=gex_cells, 
        output_file=output_guide
    )
    
    # Verify guide data structure
    log_print("Verifying guide data structure...")
    log_print(f"Guide data shape: {adata_guide.shape}")
    log_print(f"Guide features: {adata_guide.shape[1]} guides")
    
    # Mark all features as guides
    adata_guide.var["type"] = "guide"
    
    # Write outputs
    log_print(f"Final guide data shape: {adata_guide.shape}")
    log_print(f"Final cell ID format (example): {adata_guide.obs.index[0] if len(adata_guide.obs.index) > 0 else 'N/A'}")
    log_print(f"Final guide ID format (example): {adata_guide.var.index[0] if len(adata_guide.var.index) > 0 else 'N/A'}")
    log_print(f"Guide data already written to: {output_guide}")
    
    log_print(f"Writing guide metadata to {output_metadata}...")
    adata_guide.obs.to_csv(output_metadata, sep='\t')
    
    log_print(f"Successfully created guide data for pools {pool_names}!")
    log_print(f"Outputs:")
    log_print(f"  - {output_guide}")
    log_print(f"  - {output_metadata}")
    log_print(f"Guide metadata columns: {list(adata_guide.obs.columns)}")

if __name__ == "__main__":
    main()