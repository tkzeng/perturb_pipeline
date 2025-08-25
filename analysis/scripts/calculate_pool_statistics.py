#!/usr/bin/env python3
"""
Calculate pool-level statistics including undetermined read fraction.

This script counts reads directly from raw lane files for all samples in a pool
plus raw undetermined lane files, then calculates per-lane and pool-level statistics.
"""

import argparse
import pandas as pd
import yaml
import subprocess
import os
import sys
from pathlib import Path

# Add scripts directory to path to import helper functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from snakemake_helpers import load_sample_info, get_main_lane_files, get_undetermined_lane_files


def count_reads_from_fastq(fastq_file):
    """Count reads in a FASTQ file using pigz and wc"""
    cmd = f"pigz -cd {fastq_file} | echo $((`wc -l`/4))"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    assert result.returncode == 0, f"Failed to count reads in {fastq_file}: {result.stderr}"
    assert result.stdout.strip(), f"No output from read counting command for {fastq_file}"
    
    return int(result.stdout.strip())


def main():
    parser = argparse.ArgumentParser(description='Calculate pool-level statistics')
    parser.add_argument('--pool', required=True, help='Pool name')
    parser.add_argument('--sample-counts', nargs='+', required=True, 
                        help='Paths to sample count files (for processed sample totals)')
    parser.add_argument('--config', required=True, help='Path to config file')
    parser.add_argument('--output', required=True, help='Output TSV file')
    
    args = parser.parse_args()
    
    print(f"Calculating pool statistics for {args.pool}...")
    
    # Load config
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    # Get all samples in this pool
    sample_df = load_sample_info(config)
    pool_samples = sample_df[sample_df['pool'] == args.pool]
    
    assert not pool_samples.empty, f"No samples found for pool {args.pool}"
    
    # Count reads from raw undetermined lane files
    print(f"\nCounting raw undetermined lane files...")
    undetermined_total = 0
    undetermined_lanes = []
    
    undetermined_files = get_undetermined_lane_files(args.pool, config)
    for run_id, lane_base, r1_path, r2_path in undetermined_files:
        r1_count = count_reads_from_fastq(r1_path)
        undetermined_total += r1_count
        undetermined_lanes.append({
            'run_id': run_id,
            'lane': lane_base,
            'r1_path': r1_path,
            'reads': r1_count
        })
        print(f"  {lane_base}: {r1_count:,} reads")
    
    print(f"Total undetermined reads: {undetermined_total:,}")
    
    # Count reads from processed sample files (for comparison)
    print(f"\nProcessed sample totals:")
    demux_total = 0
    for counts_file in args.sample_counts:
        df = pd.read_csv(counts_file)
        sample_count = df.iloc[0]['R1_reads']
        demux_total += sample_count
        print(f"  {Path(counts_file).stem}: {sample_count:,} reads")
    
    print(f"Total demultiplexed reads: {demux_total:,}")
    
    # Calculate statistics
    total_reads = demux_total + undetermined_total
    assert total_reads > 0, f"No reads found for pool {args.pool}"
    
    undetermined_fraction = (undetermined_total / total_reads) * 100
    
    # Create output dataframe
    results = pd.DataFrame([{
        'pool': args.pool,
        'total_demultiplexed_reads': demux_total,
        'undetermined_reads': undetermined_total,
        'total_raw_reads': total_reads,
        'undetermined_fraction_percent': round(undetermined_fraction, 2),
        'num_undetermined_lanes': len(undetermined_lanes)
    }])
    
    # Write output
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    results.to_csv(args.output, sep='\t', index=False)
    print(f"\n✅ Pool statistics calculation completed")
    print(results.to_string(index=False))


if __name__ == '__main__':
    main()