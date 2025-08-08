#!/usr/bin/env python3
"""
Create undetermined FASTQ files for a single sample by concatenating matching recovery files.
"""

import argparse
import pandas as pd
import os
import glob
import subprocess
import gzip

def reverse_complement(seq):
    """Get reverse complement of sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in seq[::-1])

def load_index_mappings(primer_file):
    """Load primer name to sequence mappings from Split LP primer.xlsx"""
    primer_df = pd.read_excel(primer_file, header=None)
    seq_to_name = {}
    name_to_seq = {}
    
    for _, row in primer_df.iterrows():
        name = str(row.iloc[1]).strip()
        seq = str(row.iloc[3]).strip()
        seq_to_name[seq] = name
        name_to_seq[name] = seq
    
    return seq_to_name, name_to_seq

def find_matching_files(recovery_dir, i7_seq, i5_seq, seq_to_name):
    """Find all files matching the given i7/i5 sequences in any combination"""
    matches = []
    
    # Get all R1 files (to avoid counting twice)
    all_files = glob.glob(os.path.join(recovery_dir, "*_R1.fastq.gz"))
    
    for filepath in all_files:
        filename = os.path.basename(filepath)
        base_name = filename.replace("_R1.fastq.gz", "")
        
        # Parse filename
        parts = base_name.split("+")
        if len(parts) != 2:
            continue
            
        part1, remainder = parts[0], parts[1]
        
        # Handle orientation info
        if "_" in remainder:
            orient_parts = remainder.split("_")
            part2 = orient_parts[0]
            orientation = "_".join(orient_parts[1:])
        else:
            part2 = remainder
            orientation = "unknown"
        
        # Check if sequences match
        idx1_match = None
        idx2_match = None
        
        # Direct sequence match
        if part1 == i7_seq:
            idx1_match = seq_to_name.get(i7_seq, i7_seq)
        elif part1 == i5_seq:
            idx1_match = seq_to_name.get(i5_seq, i5_seq)
        elif part1 == "UNKNOWN":
            idx1_match = "UNKNOWN"
        
        # Same for part2
        if part2 == i7_seq:
            idx2_match = seq_to_name.get(i7_seq, i7_seq)
        elif part2 == i5_seq:
            idx2_match = seq_to_name.get(i5_seq, i5_seq)
        elif part2 == "UNKNOWN":
            idx2_match = "UNKNOWN"
        
        # Determine match type
        matched = False
        matched_by = None
        index_order = "normal"
        
        # Check normal order (i7+i5)
        if (idx1_match and idx1_match != "UNKNOWN" and 
            idx1_match == seq_to_name.get(i7_seq, i7_seq)):
            if (idx2_match and idx2_match != "UNKNOWN" and 
                idx2_match == seq_to_name.get(i5_seq, i5_seq)):
                matched_by = "both"
                matched = True
            elif idx2_match == "UNKNOWN":
                matched_by = "i7_only"
                matched = True
        elif idx1_match == "UNKNOWN":
            if (idx2_match and idx2_match != "UNKNOWN" and 
                idx2_match == seq_to_name.get(i5_seq, i5_seq)):
                matched_by = "i5_only"
                matched = True
        
        # Check swapped order (i5+i7)
        if not matched:
            if (idx1_match and idx1_match != "UNKNOWN" and 
                idx1_match == seq_to_name.get(i5_seq, i5_seq)):
                if (idx2_match and idx2_match != "UNKNOWN" and 
                    idx2_match == seq_to_name.get(i7_seq, i7_seq)):
                    matched_by = "both"
                    matched = True
                    index_order = "swapped"
                elif idx2_match == "UNKNOWN":
                    matched_by = "i5_only"
                    matched = True
                    index_order = "swapped"
            elif idx1_match == "UNKNOWN":
                if (idx2_match and idx2_match != "UNKNOWN" and 
                    idx2_match == seq_to_name.get(i7_seq, i7_seq)):
                    matched_by = "i7_only"
                    matched = True
                    index_order = "swapped"
        
        if matched:
            matches.append({
                'filename_base': base_name,
                'matched_by': matched_by,
                'i7_match': part1 if index_order == "normal" else part2,
                'i5_match': part2 if index_order == "normal" else part1,
                'orientation': orientation,
                'index_order': index_order
            })
    
    return matches

def create_dummy_fastq(output_file):
    """Create a dummy FASTQ file with reads that won't map"""
    with gzip.open(output_file, 'wt') as f:
        # Write a single dummy read (100bp)
        f.write("@DUMMY:READ:NO:MATCH 1:N:0:NNNNNNNN+NNNNNNNN\n")
        f.write("N" * 100 + "\n")
        f.write("+\n")
        f.write("!" * 100 + "\n")

def main():
    parser = argparse.ArgumentParser(description='Create undetermined FASTQ files for a single sample')
    parser.add_argument('--sample-info', required=True,
                        help='Path to sample_info.xlsx file')
    parser.add_argument('--primer-info', required=True,
                        help='Path to Split LP primer.xlsx file')
    parser.add_argument('--recovery-dir', required=True,
                        help='Path to undetermined recovery directory')
    parser.add_argument('--sample-id', required=True,
                        help='Full sample ID to process (format: pool:sample)')
    parser.add_argument('--output-r1', required=True,
                        help='Output R1 file path')
    parser.add_argument('--output-r2', required=True,
                        help='Output R2 file path')
    parser.add_argument('--min-reads', type=int, default=1000,
                        help='Minimum reads to create file (default: 1000)')
    
    args = parser.parse_args()
    
    # Load sample info
    print(f"Loading sample info for {args.sample_id}")
    sample_df = pd.read_excel(args.sample_info)
    
    # Find the sample
    sample_row = sample_df[sample_df['sample_id'] == args.sample_id]
    if sample_row.empty:
        print(f"ERROR: Sample {args.sample_id} not found in sample info")
        return 1
    
    # Get sample details
    row = sample_row.iloc[0]
    pool = row['pool']
    sample_type = row['sample_type']
    i7_seq = str(row['i7_barcode'])
    i5_seq = str(row['i5_barcode'])
    
    print(f"Sample: {args.sample_id} (pool: {pool}, type: {sample_type})")
    print(f"Indices: i7={i7_seq}, i5={i5_seq}")
    
    # Load primer mappings
    seq_to_name, name_to_seq = load_index_mappings(args.primer_info)
    
    # Load read counts from metadata if available
    metadata_file = os.path.join(args.recovery_dir, 'undetermined_recovery_metadata.tsv')
    read_counts = {}
    if os.path.exists(metadata_file):
        metadata_df = pd.read_csv(metadata_file, sep='\t')
        for _, mrow in metadata_df.iterrows():
            if '_R1.fastq.gz' in mrow['filename']:
                base = mrow['filename'].replace('_R1.fastq.gz', '')
                read_counts[base] = mrow['read_count']
    
    # Find matching files
    matches = find_matching_files(args.recovery_dir, i7_seq, i5_seq, seq_to_name)
    
    # Filter for valid matches
    valid_matches = []
    total_reads = 0
    
    for match in matches:
        # Valid if: has direct orientation for at least one index AND normal index order
        orientation = match['orientation']
        valid = ('direct' in orientation and match['index_order'] == 'normal')
        
        if valid:
            match['read_count'] = read_counts.get(match['filename_base'], 0)
            if match['read_count'] >= args.min_reads:
                valid_matches.append(match)
                total_reads += match['read_count']
    
    # Create output files
    if valid_matches:
        print(f"\nFound {len(valid_matches)} valid matches with {total_reads:,} total reads")
        
        # Collect input files
        input_files_r1 = []
        input_files_r2 = []
        
        for match in valid_matches:
            base = match['filename_base']
            r1_file = os.path.join(args.recovery_dir, f"{base}_R1.fastq.gz")
            r2_file = os.path.join(args.recovery_dir, f"{base}_R2.fastq.gz")
            
            if os.path.exists(r1_file) and os.path.exists(r2_file):
                input_files_r1.append(r1_file)
                input_files_r2.append(r2_file)
                print(f"  Including: {base} ({match['read_count']:,} reads)")
        
        # Concatenate files
        print(f"\nConcatenating {len(input_files_r1)} file pairs...")
        
        # For R1
        cmd_r1 = f"cat {' '.join(input_files_r1)} > {args.output_r1}"
        subprocess.run(cmd_r1, shell=True, check=True)
        
        # For R2
        cmd_r2 = f"cat {' '.join(input_files_r2)} > {args.output_r2}"
        subprocess.run(cmd_r2, shell=True, check=True)
        
        print(f"Created: {args.output_r1}")
        print(f"Created: {args.output_r2}")
        print(f"Total reads: {total_reads:,}")
    else:
        # No valid matches - create dummy files
        print(f"\nNo valid matches found for {args.sample_id}")
        print("Creating dummy placeholder files...")
        
        create_dummy_fastq(args.output_r1)
        create_dummy_fastq(args.output_r2)
        
        print(f"Created placeholder: {args.output_r1}")
        print(f"Created placeholder: {args.output_r2}")
    
    return 0

if __name__ == "__main__":
    exit(main())