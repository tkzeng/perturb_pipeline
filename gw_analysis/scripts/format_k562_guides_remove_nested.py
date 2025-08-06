#!/usr/bin/env python3
"""
Format K562 guides file and remove nested guides.
When a 20bp guide is nested within a 21bp guide (first 20bp match), remove the 21bp guide.
Outputs in GUIDE_REFERENCE_SPECIFICATION format with ID and gene columns.
"""

import pandas as pd
import sys

# Hardcoded paths for this simple utility script
input_file = "/oak/stanford/groups/engreitz/Users/tonyzeng/GW_PERTURB/references/guides/all_guides_duplicates_filtered.k562_yann.txt"
output_file = "/oak/stanford/groups/engreitz/Users/tonyzeng/GW_PERTURB/references/guides/k562_yann_guides_qc_reference.txt"

# Read with header
df = pd.read_csv(input_file, sep='\t')

print(f"Starting with {len(df)} guides")

# Create a dictionary to track guides by their first 20bp
guides_by_prefix = {}

for idx, row in df.iterrows():
    seq = row['GuideSequence'].upper().strip()
    guide_id = row['guideID']
    guide_len = len(seq)
    
    # Get the first 20bp as the key
    prefix_20 = seq[:20]
    
    if prefix_20 not in guides_by_prefix:
        guides_by_prefix[prefix_20] = []
    
    guides_by_prefix[prefix_20].append({
        'seq': seq,
        'id': guide_id,
        'len': guide_len
    })

# Now process each group and keep only the 20bp guide if both 20bp and 21bp exist
final_guides = []

for prefix_20, guide_list in guides_by_prefix.items():
    if len(guide_list) == 1:
        # Only one guide with this prefix, keep it
        final_guides.append(guide_list[0])
    else:
        # Multiple guides with same 20bp prefix
        # Check if we have both 20bp and 21bp versions
        lengths = [g['len'] for g in guide_list]
        
        if 20 in lengths and 21 in lengths:
            # Keep only the 20bp version(s)
            for guide in guide_list:
                if guide['len'] == 20:
                    final_guides.append(guide)
                else:
                    print(f"Removing 21bp guide: {guide['id']} ({guide['seq']}) - has 20bp version")
        else:
            # No conflict, keep all
            final_guides.extend(guide_list)

print(f"After removing nested guides: {len(final_guides)} guides")

# Check for remaining duplicates
seen_seqs = {}
duplicates = []
for guide in final_guides:
    seq = guide['seq']
    if seq in seen_seqs:
        duplicates.append((guide['id'], seen_seqs[seq], seq))
    else:
        seen_seqs[seq] = guide['id']

if duplicates:
    print(f"\nWarning: Found {len(duplicates)} duplicate sequences:")
    for dup_id, orig_id, seq in duplicates[:5]:
        print(f"  {dup_id} and {orig_id}: {seq}")
    if len(duplicates) > 5:
        print(f"  ... and {len(duplicates) - 5} more")

# Write output in GUIDE_REFERENCE_SPECIFICATION format
# Required columns: ID and gene (using ID as placeholder for gene as requested)
with open(output_file, 'w') as f:
    # Write header
    f.write("ID\tgene\n")
    # Write data - using guide ID as placeholder for gene column
    for guide in final_guides:
        f.write(f"{guide['id']}\t{guide['id']}\n")

print(f"\nWrote {len(final_guides)} guides to {output_file}")
print("Output format: GUIDE_REFERENCE_SPECIFICATION (ID and gene columns)")
print("\nFirst 5 entries:")
for guide in final_guides[:5]:
    print(f"  {guide['id']}\t{guide['id']}")