#!/usr/bin/env python3
"""Convert Yann's plate layout TSV to plate map format"""

import pandas as pd
import sys

# Read the TSV file
input_file = "/oak/stanford/groups/engreitz/Users/tonyzeng/GW_PERTURB/references/k562_yann/cell thaw split seq jun 18 2025 - platelayout.tsv"
output_file = "/oak/stanford/groups/engreitz/Users/tonyzeng/GW_PERTURB/references/k562_yann/plate_map.xlsx"

# Read the plate layout (skip the metadata rows at the bottom)
df = pd.read_csv(input_file, sep='\t', nrows=8, index_col=0)

# Create well positions and samples list
wells = []
samples = []
reps = []

# Iterate through the plate
for row_idx, row_label in enumerate(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']):
    for col in range(1, 13):  # Columns 1-12
        well = f"{row_label}{col}"
        rep_num = df.iloc[row_idx, col-1]  # The value is the replicate number
        sample = f"rep{rep_num}"
        wells.append(well)
        samples.append(sample)
        reps.append(f"rep{rep_num}")

# Create output dataframe
output_df = pd.DataFrame({
    'Well Position': wells,
    'Sample': samples,
    'Rep': reps
})

# Write to Excel
with pd.ExcelWriter(output_file) as writer:
    output_df.to_excel(writer, sheet_name='plate1', index=False)

print(f"Converted plate layout saved to: {output_file}")
print(f"Total wells: {len(wells)}")
print("\nFirst 10 entries:")
print(output_df.head(10))