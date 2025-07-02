#!/usr/bin/env python3
"""
Calculate read statistics from kallisto-bustools output files.

Calculates:
- % of mapped reads
- % of mapped reads with correctable barcodes  
- % of total reads with correctable barcodes
- PCR duplication rate from count matrix
"""

import json
import os
import sys
import argparse
from pathlib import Path
import pandas as pd
import scanpy as sc
import numpy as np
import yaml


def load_json_file(filepath):
    """Load JSON file and return contents."""
    with open(filepath, 'r') as f:
        return json.load(f)


def find_inspect_files(kb_dir):
    """Find all inspect JSON files in the kb output directory."""
    inspect_files = {}
    
    # Check main directory
    for file in Path(kb_dir).glob("*_inspect.json"):
        inspect_files[file.stem.replace("_inspect", "")] = file
    
    # Check tmp directory
    tmp_dir = Path(kb_dir) / "tmp"
    if tmp_dir.exists():
        for file in tmp_dir.glob("*_inspect.json"):
            inspect_files[file.stem.replace("_inspect", "")] = file
    
    return inspect_files


def calculate_statistics(kb_dir, sample_id=None, config_file=None, cell_barcodes_file=None):
    """Calculate read statistics from kallisto-bustools output."""
    
    # Load run_info.json
    run_info_path = Path(kb_dir) / "run_info.json"
    if not run_info_path.exists():
        raise FileNotFoundError(f"run_info.json not found in {kb_dir}")
    
    run_info = load_json_file(run_info_path)
    
    # Find all inspect files
    inspect_files = find_inspect_files(kb_dir)
    
    # Key files we need
    required_files = {
        'output': 'Initial BUS file',
        'output.unfiltered': 'Corrected BUS file',
        'output_modified.unfiltered': 'Final modified BUS file'
    }
    
    # Load inspect data
    inspect_data = {}
    for key, desc in required_files.items():
        if key not in inspect_files:
            raise FileNotFoundError(f"{key}_inspect.json not found ({desc})")
        inspect_data[key] = load_json_file(inspect_files[key])
    
    # Extract statistics from run_info.json
    total_reads = run_info['n_processed']
    mapped_reads = run_info['n_pseudoaligned']
    mapping_rate = run_info['p_pseudoaligned']
    
    # Extract statistics from inspect files
    initial_bus_reads = inspect_data['output']['numReads']
    initial_bus_records = inspect_data['output']['numRecords']
    
    # Use the final modified BUS file statistics
    if 'output_modified.unfiltered' not in inspect_data:
        raise FileNotFoundError("output_modified.unfiltered_inspect.json not found - this is required")
    
    final_stats = inspect_data['output_modified.unfiltered']
    corrected_reads = final_stats['numReads']
    corrected_records = final_stats['numRecords']
    
    # Calculate percentages
    pct_mapped = (mapped_reads / total_reads) * 100
    pct_mapped_with_valid_bc = (corrected_reads / mapped_reads) * 100 if mapped_reads > 0 else 0
    pct_usable_reads = (corrected_reads / total_reads) * 100  # % usable reads
    
    # Calculate duplication rate from count matrix if available
    duplication_rate = 'N/A'
    total_umis = 'N/A'
    reads_per_umi = 'N/A'
    total_umis_in_cells = 'N/A'
    fraction_umis_in_cells = 'N/A'
    reads_per_umi_in_cells = 'N/A'
    n_cells = 'N/A'
    
    # Check for unfiltered count matrix
    h5ad_path = Path(kb_dir) / "counts_unfiltered" / "adata.h5ad"
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Count matrix not found at {h5ad_path}")
    
    # Load count matrix
    adata = sc.read_h5ad(h5ad_path)
    
    # Sum all layers to get total counts
    if 'layers' in adata.uns_keys() or hasattr(adata, 'layers'):
        if 'mature' in adata.layers and 'nascent' in adata.layers and 'ambiguous' in adata.layers:
            total_counts_matrix = adata.layers['mature'] + adata.layers['nascent'] + adata.layers['ambiguous']
        else:
            total_counts_matrix = adata.X
    else:
        total_counts_matrix = adata.X
    
    # Calculate total UMIs
    total_umis = int(total_counts_matrix.sum())
    
    # Calculate duplication rate
    if corrected_reads > 0:
        duplication_rate = (1 - (total_umis / corrected_reads)) * 100
    
    # Calculate reads per UMI
    if total_umis > 0:
        reads_per_umi = total_reads / total_umis
    
    
    print(f"Loaded count matrix: {adata.shape}")
    print(f"Total UMIs: {total_umis:,}")
    
    # Calculate cell-specific metrics if cell barcodes provided
    if cell_barcodes_file and os.path.exists(cell_barcodes_file):
        # Load cell barcodes
        with open(cell_barcodes_file, 'r') as f:
            cell_barcodes = [line.strip() for line in f]
        
        # Filter to called cells
        cell_mask = adata.obs.index.isin(cell_barcodes)
        n_cells = sum(cell_mask)
        
        if n_cells > 0:
            # Calculate UMIs in cells only
            total_umis_in_cells = int(total_counts_matrix[cell_mask, :].sum())
            
            # Fraction of UMIs in cells
            if total_umis > 0:
                fraction_umis_in_cells = (total_umis_in_cells / total_umis) * 100
            
            # Reads per UMI in cells (assuming all reads from valid barcodes)
            if total_umis_in_cells > 0 and fraction_umis_in_cells > 0:
                # Estimate reads from cells based on UMI fraction
                reads_from_cells = corrected_reads * (fraction_umis_in_cells / 100)
                reads_per_umi_in_cells = reads_from_cells / total_umis_in_cells
            
            print(f"Cells analyzed: {n_cells:,}")
            print(f"Total UMIs in cells: {total_umis_in_cells:,}")
            print(f"Fraction UMIs in cells: {fraction_umis_in_cells:.1f}%")
    
    # Load sample info if provided
    expected_cells = 'N/A'
    min_umi_threshold = 'N/A'
    
    if sample_id and config_file:
        # Load config
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # Load sample info
        sample_info_file = config.get('sample_info_file', 'sample_info.xlsx')
        if not os.path.exists(sample_info_file):
            raise FileNotFoundError(f"Sample info file not found: {sample_info_file}")
            
        if sample_info_file.endswith('.xlsx'):
            sample_df = pd.read_excel(sample_info_file)
        else:
            sample_df = pd.read_csv(sample_info_file, sep='\t')
        
        sample_row = sample_df[sample_df['sample_id'] == sample_id]
        if sample_row.empty:
            raise ValueError(f"Sample ID {sample_id} not found in sample info file")
            
        row = sample_row.iloc[0]
        expected_cells = row.get('expected_cells', 'N/A')
        
        # Get min_umi_threshold from config or sample info
        if 'min_umi_threshold' in sample_df.columns:
            min_umi_threshold = row.get('min_umi_threshold', 'N/A')
        else:
            # Fall back to config default
            min_umi_threshold = config.get('cell_calling', {}).get('min_umi_threshold', 500)
    
    # Create results dictionary
    results = {
        'Sample ID': sample_id if sample_id else 'N/A',
        'total_reads': total_reads,
        'mapped_reads': mapped_reads,
        'reads_with_correctable_barcodes': corrected_reads,
        'pct_mapped_reads': pct_mapped,
        'pct_mapped_reads_with_correctable_barcodes': pct_mapped_with_valid_bc,
        'pct_usable_reads': pct_usable_reads,
        'total_umis': total_umis,
        'pcr_duplication_rate_pct': duplication_rate,
        'reads_per_umi': reads_per_umi,
        'expected_cells': expected_cells,
        'min_umi_threshold': min_umi_threshold
    }
    
    return results


def main():
    parser = argparse.ArgumentParser(description='Calculate read statistics from kallisto-bustools output')
    parser.add_argument('kb_dir', help='Path to kallisto-bustools output directory')
    parser.add_argument('--output', '-o', help='Output TSV file (default: stdout)')
    parser.add_argument('--sample-id', help='Full sample ID (format: pool:sample) for loading sample-specific parameters')
    parser.add_argument('--config', help='Config YAML file for loading sample info')
    parser.add_argument('--cell-barcodes', help='File containing cell barcodes for cell-specific metrics')
    
    args = parser.parse_args()
    
    try:
        stats = calculate_statistics(args.kb_dir, args.sample_id, args.config, args.cell_barcodes)
        
        # Create TSV output
        df = pd.DataFrame([stats])
        output = df.to_csv(sep='\t', index=False)
        
        # Write output
        if args.output:
            with open(args.output, 'w') as f:
                f.write(output)
            print(f"Statistics written to {args.output}")
        else:
            print(output)
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()