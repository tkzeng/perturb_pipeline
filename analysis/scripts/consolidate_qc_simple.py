#!/usr/bin/env python3
"""
Consolidate QC metrics across all samples using simple concatenation.

All merging of read statistics and guide statistics is already done at the 
individual sample level by calculate_qc_metrics_by_biological_sample.

This script produces three consolidated output files:

1. all_metrics.tsv: Sample-level metrics
   - One row per sequencing sample (e.g., gex_1, gex_2)
   - Includes pool, sample_id, all read statistics (GEX and guide)
   - For each cell calling method: n_cells, UMI/gene/mito metrics, guides per cell
   - Most comprehensive view for sequencing QC
   
2. by_biological_sample.tsv: Biological replicate view
   - One row per biological sample per sequencing run
   - Only QC metrics (no read statistics)
   - Shows variation in cell quality across technical replicates
   - Same biological sample may appear multiple times if in different pools/runs
   
3. by_well.tsv: Well-level view
   - One row per well per sequencing run  
   - Only QC metrics (no read statistics)
   - Useful for identifying systematic plate position effects
   - Can spot problematic wells affecting biological samples

Works for any QC metrics files.
"""

import pandas as pd
import sys
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Simple consolidation of QC metrics files")
    parser.add_argument("--input-files", nargs='+', required=True,
                        help="List of TSV files to consolidate")
    parser.add_argument("--output", required=True,
                        help="Output consolidated TSV file")
    parser.add_argument("--add-pool", action='store_true',
                        help="Add pool column based on file path")
    
    args = parser.parse_args()
    
    # Simple consolidation - just concatenate all input files
    print(f"Consolidating {len(args.input_files)} files...")
    dfs = []
    
    for file in args.input_files:
        df = pd.read_csv(file, sep='\t')
        
        # Optionally add pool from file path
        if args.add_pool and 'pool' not in df.columns:
            # Extract pool from file path
            parts = file.split('/')
            for i, part in enumerate(parts):
                if part.startswith('pool'):
                    df['pool'] = part
                    break
        
        dfs.append(df)
    
    # Concatenate all dataframes
    consolidated_df = pd.concat(dfs, ignore_index=True)
    
    # Save consolidated file
    consolidated_df.to_csv(args.output, sep='\t', index=False)
    print(f"Saved consolidated metrics to {args.output}")
    print(f"Total rows: {len(consolidated_df)}")
    print(f"Total columns: {len(consolidated_df.columns)}")


if __name__ == "__main__":
    main()