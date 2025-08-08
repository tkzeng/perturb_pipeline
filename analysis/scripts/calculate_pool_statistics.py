#!/usr/bin/env python3
"""
Calculate pool-level statistics including undetermined read fraction.

This script reads count files for all samples in a pool plus the undetermined count,
then calculates the fraction of reads that were undetermined.
"""

import argparse
import pandas as pd
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description='Calculate pool-level statistics')
    parser.add_argument('--pool', required=True, help='Pool name')
    parser.add_argument('--sample-counts', nargs='+', required=True, 
                        help='Paths to sample count files')
    parser.add_argument('--undetermined-counts', required=True, 
                        help='Path to undetermined count file')
    parser.add_argument('--output', required=True, help='Output TSV file')
    
    args = parser.parse_args()
    
    print(f"Calculating pool statistics for {args.pool}...")
    
    # Read undetermined count
    undetermined_df = pd.read_csv(args.undetermined_counts)
    undetermined_count = undetermined_df.iloc[0]['R1_reads']
    print(f"Undetermined reads: {undetermined_count:,}")
    
    # Sum demultiplexed reads from all samples
    demux_total = 0
    for counts_file in args.sample_counts:
        df = pd.read_csv(counts_file)
        sample_count = df.iloc[0]['R1_reads']
        demux_total += sample_count
        print(f"  Added {sample_count:,} reads from {Path(counts_file).name}")
    
    print(f"Total demultiplexed reads: {demux_total:,}")
    
    # Calculate undetermined fraction
    total_reads = demux_total + undetermined_count
    if total_reads > 0:
        undetermined_fraction = (undetermined_count / total_reads) * 100
    else:
        undetermined_fraction = 0.0
    
    # Create output dataframe
    results = pd.DataFrame([{
        'pool': args.pool,
        'total_demultiplexed_reads': demux_total,
        'undetermined_reads': undetermined_count,
        'undetermined_fraction_percent': round(undetermined_fraction, 2)
    }])
    
    # Write output
    results.to_csv(args.output, sep='\t', index=False)
    print(f"\n✅ Pool statistics calculation completed")
    print(results.to_string(index=False))


if __name__ == '__main__':
    main()