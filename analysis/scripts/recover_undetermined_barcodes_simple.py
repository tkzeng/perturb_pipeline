#!/usr/bin/env python3
"""
Simple undetermined read recovery - outputs exactly which indices matched.
"""

import argparse
import pandas as pd
import pysam
import gzip
import os
import sys
from collections import defaultdict
import time

def reverse_complement(seq):
    """Get reverse complement of sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in seq[::-1])

def load_indices(indices_file):
    """Load all indices into separate dictionaries for direct and RC."""
    print(f"Loading indices from {indices_file}")
    
    df = pd.read_excel(indices_file, header=None)
    
    # Separate dictionaries for direct and reverse complement
    index_direct = {}  # Direct sequences
    index_rc = {}      # Reverse complement sequences
    
    for _, row in df.iterrows():
        name = str(row.iloc[1]).strip()  # Column 1: name (Split P5_A1, etc.)
        index_seq = str(row.iloc[3]).strip().upper()  # Column 3: index sequence
        
        # Crash if invalid data
        assert pd.notna(name), f"Name is NaN in row {_}"
        assert pd.notna(index_seq), f"Index sequence is NaN in row {_}"
        
        index_direct[index_seq] = name
        index_rc[reverse_complement(index_seq)] = name
    
    print(f"Loaded {len(index_direct)} direct indices and {len(index_rc)} RC indices")
    return index_direct, index_rc

def test_matching(index_direct, index_rc):
    """Test the matching with user examples."""
    print("\nRunning tests...")
    print("="*60)
    
    test_cases = [
        ("GAACTGAGCG", "GGGGGGGGGG"),
        ("GGGGGGGGGG", "CGCTCCACGA")
    ]
    
    for idx1, idx2 in test_cases:
        # Check both direct and RC
        name1 = "UNKNOWN"
        name2 = "UNKNOWN"
        
        if idx1 in index_direct:
            name1 = index_direct[idx1]
        elif idx1 in index_rc:
            name1 = index_rc[idx1] + "_RC"
            
        if idx2 in index_direct:
            name2 = index_direct[idx2]
        elif idx2 in index_rc:
            name2 = index_rc[idx2] + "_RC"
        
        print(f"\nTest: {idx1} + {idx2}")
        print(f"  Index1: {idx1} -> {name1}")
        print(f"  Index2: {idx2} -> {name2}")
        
        if name1 == "UNKNOWN" and name2 == "UNKNOWN":
            print("  ✗ FAIL: No indices matched")
        else:
            print("  ✓ PASS: At least one index matched")
    
    print("="*60)

def process_reads(r1_file, r2_file, index_direct, index_rc, output_dir):
    """Process reads and output based on index matches."""
    
    # Open input files
    fq_r1 = pysam.FastxFile(r1_file)
    fq_r2 = pysam.FastxFile(r2_file)
    
    # Track all unique index combinations seen
    seen_combinations = set()
    
    # Dictionary to store file handles
    output_handles = {}
    
    # Dictionary to store metadata for each output
    metadata = {}
    
    # Stats
    total = 0
    recovered = 0
    match_counts = defaultdict(int)
    
    # Buffer for combinations that haven't hit 100 reads yet
    pending_reads = defaultdict(list)  # output_name -> [(r1_data, r2_data), ...]
    
    # Timing
    start_time = time.time()
    last_log_time = time.time()
    
    print("\nProcessing reads...")
    
    # Process reads
    for r1, r2 in zip(fq_r1, fq_r2):
        total += 1
        
        
        # Extract indices from comment field - crash if invalid
        parts = r1.comment.split(':')
        assert len(parts) >= 4, f"Invalid comment format: {r1.comment}"
        assert '+' in parts[-1], f"No '+' in indices: {parts[-1]}"
        
        indices = parts[-1].split('+')
        assert len(indices) == 2, f"Expected 2 indices, got {len(indices)}: {indices}"
        
        idx1, idx2 = indices[0], indices[1]
        
        # Look up index names - check both direct and RC
        name1 = "UNKNOWN"
        name2 = "UNKNOWN"
        idx1_orientation = "unknown"
        idx2_orientation = "unknown"
        
        if idx1 in index_direct:
            name1 = index_direct[idx1]
            idx1_orientation = "direct"
        elif idx1 in index_rc:
            name1 = index_rc[idx1]
            idx1_orientation = "rc"
            
        if idx2 in index_direct:
            name2 = index_direct[idx2]
            idx2_orientation = "direct"
        elif idx2 in index_rc:
            name2 = index_rc[idx2]
            idx2_orientation = "rc"
        
        # Skip if both unknown
        if name1 == "UNKNOWN" and name2 == "UNKNOWN":
            continue
        
        # Create output file name - use "UNKNOWN" for unmatched indices
        idx1_for_filename = idx1 if name1 != "UNKNOWN" else "UNKNOWN"
        idx2_for_filename = idx2 if name2 != "UNKNOWN" else "UNKNOWN"
        output_name = f"{idx1_for_filename}+{idx2_for_filename}_{idx1_orientation}_{idx2_orientation}"
        recovered += 1
        match_counts[output_name] += 1
        
        # Store metadata for later
        if output_name not in metadata:
            metadata[output_name] = {
                'idx1_seq': idx1,
                'idx2_seq': idx2,
                'idx1_name': name1,
                'idx2_name': name2,
                'idx1_orientation': idx1_orientation,
                'idx2_orientation': idx2_orientation
            }
        
        # If files don't exist yet, check if we should create them
        if output_name not in output_handles:
            # Buffer reads until we hit 1000
            if match_counts[output_name] < 1000:
                # Store read data for later
                r1_data = f"@{r1.name} {r1.comment}\n{r1.sequence}\n+\n{r1.quality}\n"
                r2_data = f"@{r2.name} {r2.comment}\n{r2.sequence}\n+\n{r2.quality}\n"
                pending_reads[output_name].append((r1_data, r2_data))
                continue  # Don't write yet
            else:
                # Hit 1000 reads! Create files and write all pending reads
                seen_combinations.add(output_name)
                r1_path = os.path.join(output_dir, f"{output_name}_R1.fastq.gz")
                r2_path = os.path.join(output_dir, f"{output_name}_R2.fastq.gz")
                
                # Open gzip files with fast compression (level 1)
                r1_out = gzip.open(r1_path, 'wt', compresslevel=1)
                r2_out = gzip.open(r2_path, 'wt', compresslevel=1)
                
                output_handles[output_name] = (r1_out, r2_out)
                
                # Write all pending reads
                r1_out, r2_out = output_handles[output_name]
                for r1_data, r2_data in pending_reads[output_name]:
                    r1_out.write(r1_data)
                    r2_out.write(r2_data)
                
                # Clear the buffer
                del pending_reads[output_name]
                
                # Write the current (1000th) read
                r1_out.write(f"@{r1.name} {r1.comment}\n{r1.sequence}\n+\n{r1.quality}\n")
                r2_out.write(f"@{r2.name} {r2.comment}\n{r2.sequence}\n+\n{r2.quality}\n")
                
                # Flush immediately after file creation to ensure 1000 reads are visible
                r1_out.flush()
                r2_out.flush()
        else:
            # Files already exist, write directly (no buffering after 1000)
            r1_out, r2_out = output_handles[output_name]
            r1_out.write(f"@{r1.name} {r1.comment}\n{r1.sequence}\n+\n{r1.quality}\n")
            r2_out.write(f"@{r2.name} {r2.comment}\n{r2.sequence}\n+\n{r2.quality}\n")
            
            # Periodic flush every 1000 reads
            if match_counts[output_name] % 1000 == 0:
                r1_out.flush()
                r2_out.flush()
        
        # Progress report every 5 seconds
        current_time = time.time()
        if current_time - last_log_time >= 5:
            elapsed = current_time - start_time
            rate = total / elapsed
            recovery_rate = recovered / total * 100 if total > 0 else 0
            print(f"Processed {total:,} reads ({rate:.0f} reads/sec) - Recovered {recovered:,} ({recovery_rate:.1f}%)")
            sys.stdout.flush()  # Ensure immediate output
            last_log_time = current_time
    
    # Close everything
    fq_r1.close()
    fq_r2.close()
    
    for name, (r1_out, r2_out) in output_handles.items():
        r1_out.close()
        r2_out.close()
    
    elapsed = time.time() - start_time
    
    # Print stats
    print(f"\nCompleted in {elapsed/60:.1f} minutes ({total/elapsed:.0f} reads/sec)")
    print(f"\nRecovery Statistics:")
    print(f"Total reads: {total:,}")
    print(f"Successfully recovered: {recovered:,} ({recovered/total*100:.2f}%)")
    print(f"\nTop 20 index combinations:")
    for combo, count in sorted(match_counts.items(), key=lambda x: x[1], reverse=True)[:20]:
        print(f"  {combo}: {count:,}")
    
    # Write summary
    with open(os.path.join(output_dir, "recovery_summary.txt"), 'w') as f:
        f.write(f"Total reads: {total}\n")
        f.write(f"Recovered: {recovered}\n")
        f.write(f"Recovery rate: {recovered/total*100:.2f}%\n\n")
        f.write("Index combinations:\n")
        for combo, count in sorted(match_counts.items(), key=lambda x: x[1], reverse=True):
            f.write(f"{combo}\t{count}\n")
    
    # Write detailed metadata file
    with open(os.path.join(output_dir, "undetermined_recovery_metadata.tsv"), 'w') as f:
        # Write header
        f.write("filename\tidx1_seq\tidx2_seq\tidx1_name\tidx2_name\tidx1_orientation\tidx2_orientation\tread_count\n")
        
        # Write data for each file created
        for output_name, count in sorted(match_counts.items(), key=lambda x: x[1], reverse=True):
            if output_name in metadata:
                meta = metadata[output_name]
                # Write entries for both R1 and R2 files
                for read in ['R1', 'R2']:
                    filename = f"{output_name}_{read}.fastq.gz"
                    f.write(f"{filename}\t{meta['idx1_seq']}\t{meta['idx2_seq']}\t")
                    f.write(f"{meta['idx1_name']}\t{meta['idx2_name']}\t")
                    f.write(f"{meta['idx1_orientation']}\t{meta['idx2_orientation']}\t{count}\n")

def main():
    parser = argparse.ArgumentParser(description='Simple undetermined recovery')
    parser.add_argument('--fastq-r1', required=True, help='R1 FASTQ file')
    parser.add_argument('--fastq-r2', required=True, help='R2 FASTQ file')
    parser.add_argument('--indices', required=True, help='Index file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load indices
    index_direct, index_rc = load_indices(args.indices)
    
    # Run tests
    test_matching(index_direct, index_rc)
    
    # Process reads
    process_reads(args.fastq_r1, args.fastq_r2, index_direct, index_rc, args.output_dir)

if __name__ == "__main__":
    main()