#!/usr/bin/env python3
"""
Validate guide file format before kallisto kite index building.

Requirements:
- File must have no headers
- File must be tab delimited
- No duplicate guides (first column)
- No duplicate names (second column)
- No spaces in names (second column)
"""

import sys
import argparse
import pandas as pd
import re


def validate_guide_file(guide_file_path):
    """
    Validate guide file format.
    
    Parameters
    ----------
    guide_file_path : str
        Path to the guide file to validate
    
    Raises
    ------
    ValueError
        If any validation check fails
    """
    df = pd.read_csv(guide_file_path, sep='\t', header=None, dtype=str)
    
    if df.shape[1] < 2:
        raise ValueError("File must have at least 2 tab-delimited columns")
    
    guide_sequences = df[0]
    guide_names = df[1]
    
    # Check if first row looks like a header (first column not only ACTG)
    if not bool(re.match('^[ACTGactg]+$', guide_sequences.iloc[0])):
        raise ValueError(f"First row appears to be a header (guide sequence contains non-ACTG characters): '{guide_sequences.iloc[0]}'")
    
    # Check for duplicate guide sequences
    duplicated_guides = guide_sequences[guide_sequences.duplicated()]
    if not duplicated_guides.empty:
        raise ValueError(f"Duplicate guide sequences found: {duplicated_guides.unique().tolist()}")
    
    # Check for duplicate guide names
    duplicated_names = guide_names[guide_names.duplicated()]
    if not duplicated_names.empty:
        raise ValueError(f"Duplicate guide names found: {duplicated_names.unique().tolist()}")
    
    # Check for spaces in guide names
    names_with_spaces = guide_names[guide_names.str.contains(' ', na=False)]
    if not names_with_spaces.empty:
        raise ValueError(f"Guide names with spaces found: {names_with_spaces.tolist()}")
    
    if len(guide_sequences) == 0:
        raise ValueError("File is empty or contains no valid guide entries")
    
    print(f"✓ Guide file validation passed: {len(guide_sequences)} unique guides found")


def main():
    parser = argparse.ArgumentParser(description='Validate guide file format for kallisto kite')
    parser.add_argument('guide_file', help='Path to guide file to validate')
    
    args = parser.parse_args()
    
    try:
        validate_guide_file(args.guide_file)
        sys.exit(0)
    except ValueError as e:
        print(f"ERROR: Guide file validation failed: {e}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print(f"ERROR: Guide file not found: {args.guide_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Unexpected error during validation: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()