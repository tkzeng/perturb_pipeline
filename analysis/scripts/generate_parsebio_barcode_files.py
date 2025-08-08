#!/usr/bin/env python3
"""
Generate barcodes.XX.txt and replace.XX.txt files from Parse Bio CSV files.

This script converts Parse Bio's bc_data CSV files into the format needed by the pipeline:
- barcodes.XX.txt: Three-column whitelist (BC3, BC2, BC1)
- replace.XX.txt: Maps random hexamer (R-type) to poly-T (T-type) barcodes

Supported kits:
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
    'v2': {
        'WT_mini': {
            'bc1_csv': 'bc_data_n24_v4.csv',
            'bc2_csv': 'bc_data_v1.csv',
            'bc3_csv': 'bc_data_v1.csv',
            'output_suffix': '12_v2',
            'num_wells': 12
        },
        'WT': {
            'bc1_csv': 'bc_data_n99_v5.csv',
            'bc2_csv': 'bc_data_v1.csv',
            'bc3_csv': 'bc_data_v1.csv',
            'output_suffix': '48_v2',
            'num_wells': 48
        },
        'WT_mega': {
            'bc1_csv': 'bc_data_n198_v5.csv',
            'bc2_csv': 'bc_data_v1.csv',
            'bc3_csv': 'bc_data_v1.csv',
            'output_suffix': '96_v2',
            'num_wells': 96
        }
    },
    'v3': {
        'WT_mini': {
            'bc1_csv': 'bc_data_n37_R1_v3_6.csv',
            'bc2_csv': 'bc_data_v1.csv',
            'bc3_csv': 'bc_data_R3_v3.csv',
            'output_suffix': '12_v3',
            'num_wells': 12
        },
        'WT': {
            'bc1_csv': 'bc_data_n141_R1_v3_6.csv',
            'bc2_csv': 'bc_data_v1.csv',
            'bc3_csv': 'bc_data_R3_v3.csv',
            'output_suffix': '48_v3',
            'num_wells': 48
        },
        'WT_mega': {
            'bc1_csv': 'bc_data_n299_R1_v3_6.csv',
            'bc2_csv': 'bc_data_v1.csv',
            'bc3_csv': 'bc_data_R3_v3.csv',
            'output_suffix': '96_v3',
            'num_wells': 96
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
    
    # Get the number of wells for this kit
    num_wells = config['num_wells']
    
    # Generate barcodes file
    suffix = config['output_suffix']
    barcodes_file = Path(output_dir) / f"barcodes.{suffix}.txt"
    
    with open(barcodes_file, 'w') as f:
        # Write T-type barcodes (poly-T primed)
        # Format: BC3 BC2 BC1
        for idx, row in df_bc1_t.iterrows():
            bc1 = row['sequence']
            # Use the barcode index to get corresponding BC2/BC3
            bci = int(row['bci']) - 1  # Convert to 0-based index
            
            # Get BC2 and BC3 for this index (cycling through if needed)
            bc2_idx = bci % len(df_bc2)
            bc3_idx = bci % len(df_bc3)
            
            bc2 = df_bc2.iloc[bc2_idx]['sequence']
            bc3 = df_bc3.iloc[bc3_idx]['sequence']
            
            f.write(f"{bc3} {bc2} {bc1}\n")
        
        # Write R-type barcodes (random hexamer primed)
        # Format: - - BC1
        for idx, row in df_bc1_r.iterrows():
            bc1 = row['sequence']
            f.write(f"- - {bc1}\n")
    
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
        ('v2', 'WT_mini'),
        ('v2', 'WT'),
        ('v2', 'WT_mega'),
        ('v3', 'WT_mini'),
        ('v3', 'WT'),
        ('v3', 'WT_mega'),
    ]
    
    success_count = 0
    for chemistry, kit_name in kits_to_process:
        try:
            if process_kit(chemistry, kit_name, parsebio_dir, output_dir):
                success_count += 1
        except Exception as e:
            print(f"Error processing {kit_name} (chemistry {chemistry}): {e}")
    
    print(f"\nSuccessfully generated files for {success_count}/{len(kits_to_process)} kits")

if __name__ == "__main__":
    main()