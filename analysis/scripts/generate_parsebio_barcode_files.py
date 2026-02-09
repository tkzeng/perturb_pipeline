#!/usr/bin/env python3
"""
Generate barcodes.XX.txt and replace.XX.txt files from Parse Bio CSV files.

This script converts Parse Bio's bc_data CSV files into the format needed by the pipeline:
- barcodes.XX.txt: Three-column whitelist (BC3, BC2, BC1)
- replace.XX.txt: Maps random hexamer (R-type) to poly-T (T-type) barcodes

Supported kits:
- Chemistry v1: WT_mini (12-well), WT (48-well), WT_mega (96-well)
- Chemistry v2: WT_mini (12-well), WT (48-well), WT_mega (96-well)
- Chemistry v3: WT_mini (12-well), WT (48-well), WT_mega (96-well)

Kit definitions are based on:
/oak/stanford/groups/engreitz/Users/tonyzeng/GW_PERTURB/references/ParseBiosciences-Pipeline.1.6.0/splitpipe/kits.py
"""

import pandas as pd
import sys
import os
from pathlib import Path

# Kit configuration mapping
KIT_CONFIGS = {
    'v1': {
        'WT_mini': {
            'bc1_csv': 'bc_data_n24_v4.csv',
            'bc2_csv': 'bc_data_v1.csv',
            'bc3_csv': 'bc_data_v1.csv',
            'output_suffix': '12_v1',
        },
        'WT': {
            'bc1_csv': 'bc_data_v2.csv',
            'bc2_csv': 'bc_data_v1.csv',
            'bc3_csv': 'bc_data_v1.csv',
            'output_suffix': '48_v1',
        },
        'WT_mega': {
            'bc1_csv': 'bc_data_n198_v5.csv',
            'bc2_csv': 'bc_data_v1.csv',
            'bc3_csv': 'bc_data_v1.csv',
            'output_suffix': '96_v1',
        }
    },
    'v2': {
        'WT_mini': {
            'bc1_csv': 'bc_data_n24_v4.csv',
            'bc2_csv': 'bc_data_v1.csv',
            'bc3_csv': 'bc_data_v1.csv',
            'output_suffix': '12_v2',
        },
        'WT': {
            'bc1_csv': 'bc_data_n99_v5.csv',
            'bc2_csv': 'bc_data_v1.csv',
            'bc3_csv': 'bc_data_v1.csv',
            'output_suffix': '48_v2',
        },
        'WT_mega': {
            'bc1_csv': 'bc_data_n198_v5.csv',
            'bc2_csv': 'bc_data_v1.csv',
            'bc3_csv': 'bc_data_v1.csv',
            'output_suffix': '96_v2',
        }
    },
    'v3': {
        'WT_mini': {
            'bc1_csv': 'bc_data_n37_R1_v3_6.csv',
            'bc2_csv': 'bc_data_v1.csv',
            'bc3_csv': 'bc_data_R3_v3.csv',
            'output_suffix': '12_v3',
        },
        'WT': {
            'bc1_csv': 'bc_data_n141_R1_v3_6.csv',
            'bc2_csv': 'bc_data_v1.csv',
            'bc3_csv': 'bc_data_R3_v3.csv',
            'output_suffix': '48_v3',
        },
        'WT_mega': {
            'bc1_csv': 'bc_data_n299_R1_v3_6.csv',
            'bc2_csv': 'bc_data_v1.csv',
            'bc3_csv': 'bc_data_R3_v3.csv',
            'output_suffix': '96_v3',
        }
    }
}

def read_barcode_csv(csv_path, barcode_type=None):
    """
    Read a Parse Bio barcode CSV file.
    
    Args:
        csv_path: Path to CSV file
        barcode_type: Filter by stype if specified ('T', 'R', 'L', etc.)
    
    Returns:
        DataFrame with barcode information
    """
    df = pd.read_csv(csv_path)
    if barcode_type:
        df = df[df['stype'] == barcode_type].copy()
    return df

def process_kit(chemistry, kit_name, parsebio_dir, output_dir):
    """
    Process a specific kit configuration to generate barcode files.
    
    Args:
        chemistry: Chemistry version ('v2' or 'v3')
        kit_name: Kit name ('WT_mini', 'WT', or 'WT_mega')
        parsebio_dir: Path to Parse Bio pipeline directory with CSV files
        output_dir: Output directory for generated files
    """
    if chemistry not in KIT_CONFIGS:
        raise ValueError(f"Chemistry {chemistry} not supported")
    
    if kit_name not in KIT_CONFIGS[chemistry]:
        raise ValueError(f"Kit {kit_name} not supported for chemistry {chemistry}")
    
    config = KIT_CONFIGS[chemistry][kit_name]
    barcode_dir = Path(parsebio_dir) / 'splitpipe' / 'barcodes'
    
    # Read BC1 (Round 1) barcodes - contains both T and R types
    bc1_csv_path = barcode_dir / config['bc1_csv']
    if not bc1_csv_path.exists():
        print(f"Error: BC1 CSV file not found: {bc1_csv_path}")
        return False
    
    df_bc1 = pd.read_csv(bc1_csv_path)
    df_bc1_t = df_bc1[df_bc1['stype'] == 'T'].copy()
    df_bc1_r = df_bc1[df_bc1['stype'] == 'R'].copy() if 'R' in df_bc1['stype'].values else pd.DataFrame()
    
    # Read BC2 (Round 2) barcodes
    bc2_csv_path = barcode_dir / config['bc2_csv']
    if not bc2_csv_path.exists():
        print(f"Error: BC2 CSV file not found: {bc2_csv_path}")
        return False
    
    df_bc2 = read_barcode_csv(bc2_csv_path, 'L')  # Ligation barcodes
    
    # Read BC3 (Round 3) barcodes
    bc3_csv_path = barcode_dir / config['bc3_csv']
    if not bc3_csv_path.exists():
        print(f"Error: BC3 CSV file not found: {bc3_csv_path}")
        return False
    
    df_bc3 = read_barcode_csv(bc3_csv_path, 'L')  # Ligation barcodes
    
    print(f"\nProcessing {kit_name} (chemistry {chemistry}):")
    print(f"  BC1 file: {config['bc1_csv']}")
    print(f"    T-type barcodes: {len(df_bc1_t)}")
    print(f"    R-type barcodes: {len(df_bc1_r)}")
    print(f"  BC2 file: {config['bc2_csv']} ({len(df_bc2)} barcodes)")
    print(f"  BC3 file: {config['bc3_csv']} ({len(df_bc3)} barcodes)")
    
    # Generate barcodes file
    # The three columns (BC3, BC2, BC1) are independent whitelists, one per
    # barcoding round. Each round always uses all of its barcodes regardless
    # of kit size (e.g. R2 and R3 always have 96 barcodes, even for 12-well
    # kits). Columns shorter than the longest are padded with '-' (wildcard).
    suffix = config['output_suffix']
    barcodes_file = Path(output_dir) / f"barcodes.{suffix}.txt"

    bc3_list = list(df_bc3['sequence'])
    bc2_list = list(df_bc2['sequence'])
    bc1_list = list(df_bc1_t['sequence']) + list(df_bc1_r['sequence'])

    max_len = max(len(bc3_list), len(bc2_list), len(bc1_list))
    bc3_list += ['-'] * (max_len - len(bc3_list))
    bc2_list += ['-'] * (max_len - len(bc2_list))
    bc1_list += ['-'] * (max_len - len(bc1_list))

    with open(barcodes_file, 'w') as f:
        for bc3, bc2, bc1 in zip(bc3_list, bc2_list, bc1_list):
            f.write(f"{bc3} {bc2} {bc1}\n")
    
    print(f"  Generated: {barcodes_file}")
    
    # Generate replace file (maps R-type to T-type for same well)
    replace_file = Path(output_dir) / f"replace.{suffix}.txt"
    
    replacements = []
    if not df_bc1_r.empty:
        # Create mapping from well to T-type barcode
        well_to_t = {}
        for idx, row in df_bc1_t.iterrows():
            well_to_t[row['well']] = row['sequence']
        
        # Map each R-type to corresponding T-type by well
        for idx, row in df_bc1_r.iterrows():
            r_barcode = row['sequence']
            well = row['well']
            if well in well_to_t:
                t_barcode = well_to_t[well]
                replacements.append((r_barcode, t_barcode))
    
    with open(replace_file, 'w') as f:
        for r_bc, t_bc in replacements:
            f.write(f"{r_bc} *{t_bc}\n")
    
    print(f"  Generated: {replace_file} ({len(replacements)} R->T mappings)")
    
    return True

def main():
    if len(sys.argv) < 2:
        print("Usage: python generate_parsebio_barcode_files.py <parsebio_pipeline_dir> [output_dir]")
        print("\nExample:")
        print("  python generate_parsebio_barcode_files.py ParseBiosciences-Pipeline.1.2.1/ ./")
        print("\nSupported kits:")
        print("  Chemistry v2: WT_mini (12-well), WT (48-well), WT_mega (96-well)")
        print("  Chemistry v3: WT_mini (12-well), WT_mega (96-well)")
        print("\nThis script generates:")
        print("  - barcodes.XX.txt: Three-column barcode whitelist (BC3, BC2, BC1)")
        print("  - replace.XX.txt: Maps random hexamer (R) to poly-T (T) barcodes")
        sys.exit(1)
    
    parsebio_dir = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "."
    
    # Check if Parse Bio directory exists
    if not Path(parsebio_dir).exists():
        print(f"Error: Parse Bio directory not found: {parsebio_dir}")
        sys.exit(1)
    
    # Create output directory if needed
    Path(output_dir).mkdir(exist_ok=True)
    
    # Process all supported kits
    kits_to_process = [
        ('v1', 'WT_mini'),
        ('v1', 'WT'),
        ('v1', 'WT_mega'),
        ('v2', 'WT_mini'),
        ('v2', 'WT'),
        ('v2', 'WT_mega'),
        ('v3', 'WT_mini'),
        ('v3', 'WT'),
        ('v3', 'WT_mega'),
    ]
    
    for chemistry, kit_name in kits_to_process:
        process_kit(chemistry, kit_name, parsebio_dir, output_dir)

if __name__ == "__main__":
    main()