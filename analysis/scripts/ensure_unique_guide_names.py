#!/usr/bin/env python3
"""
Ensure unique guide names in a guide TSV file.
If duplicates exist, append a counter suffix.
"""

import sys
import pandas as pd
from collections import defaultdict

def main(input_file, output_file):
    # Read guide file
    df = pd.read_csv(input_file, sep='\t', header=None, names=['sequence', 'name'])
    
    # Track name occurrences
    name_counts = defaultdict(int)
    unique_names = []
    n_duplicates = 0
    
    for original_name in df['name']:
        name_counts[original_name] += 1
        if name_counts[original_name] == 1:
            unique_names.append(original_name)
        else:
            # Append counter for duplicates
            unique_names.append(f"{original_name}_{name_counts[original_name]}")
            n_duplicates += 1
    
    # Update dataframe
    df['name'] = unique_names
    
    # Check for duplicate sequences
    dup_seqs = df[df.duplicated(subset='sequence', keep=False)]
    if not dup_seqs.empty:
        print("WARNING: Duplicate sequences found:", file=sys.stderr)
        print(dup_seqs[['sequence', 'name']].to_string(), file=sys.stderr)
    
    # Report changes
    if len(set(df['name'])) < len(df):
        print(f"ERROR: Still have duplicate names after renaming!", file=sys.stderr)
        sys.exit(1)
    
    # Count unique original names that had duplicates
    n_unique_names_with_dups = sum(1 for count in name_counts.values() if count > 1)
    print(f"Found {n_unique_names_with_dups} guide names with duplicates")
    print(f"Renamed {n_duplicates} duplicate occurrences")
    
    # Save
    df.to_csv(output_file, sep='\t', header=False, index=False)
    print(f"Saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python ensure_unique_guide_names.py input.tsv output.tsv")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2])