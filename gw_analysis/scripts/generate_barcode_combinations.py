#!/usr/bin/env python3

import sys
import itertools

def generate_barcode_combinations(input_file, output_file):
    """
    Generate all possible combinations of barcodes from a 3-column whitelist file.
    
    Args:
        input_file: Path to the input file with 3 columns of barcodes
        output_file: Path to the output file with concatenated barcode combinations
    """
    # Read the three columns
    col1_barcodes = set()
    col2_barcodes = set()
    col3_barcodes = set()
    
    with open(input_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                # Skip entries with '-' as they represent missing barcodes
                if parts[0] != '-':
                    col1_barcodes.add(parts[0])
                if parts[1] != '-':
                    col2_barcodes.add(parts[1])
                if parts[2] != '-':
                    col3_barcodes.add(parts[2])
    
    print(f"Column 1: {len(col1_barcodes)} unique barcodes")
    print(f"Column 2: {len(col2_barcodes)} unique barcodes")
    print(f"Column 3: {len(col3_barcodes)} unique barcodes")
    
    # Generate all combinations
    all_combinations = []
    for bc1, bc2, bc3 in itertools.product(col1_barcodes, col2_barcodes, col3_barcodes):
        combined = bc1 + bc2 + bc3
        all_combinations.append(combined)
    
    print(f"Total combinations: {len(all_combinations)}")
    
    # Write to output file
    with open(output_file, 'w') as f:
        for barcode in sorted(all_combinations):
            f.write(barcode + '\n')
    
    print(f"Written {len(all_combinations)} barcodes to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_barcode_combinations.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    generate_barcode_combinations(input_file, output_file)