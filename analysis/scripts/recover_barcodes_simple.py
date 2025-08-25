#!/usr/bin/env python3
"""
Simple barcode recovery using direct whitelist matching.
Searches ±3 positions around expected barcode locations.
"""

import argparse
import gzip
import sys
from pathlib import Path
import logging
import pysam
import time
import subprocess

# Constants
SPACER1 = "GTGGCCGATGTTTCGCATCGGCGTACGACT"
SPACER2 = "ATCCACGTGCTTGAGACTGTGG"

# Expected barcode positions in read structure:
# [UMI:10][BC3:8][SPACER1:30][BC2:8][SPACER2:22][BC1:8]
# Positions: 0-9, 10-17, 18-47, 48-55, 56-77, 78-85


def load_barcode_whitelists(barcode_file):
    """Load barcode whitelists from file, handling '-' entries."""
    bc1_set = set()
    bc2_set = set()
    bc3_set = set()
    
    with open(barcode_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                # Map columns to their position in the read structure
                # Column 1 -> BC3 (leftmost, next to UMI)
                # Column 2 -> BC2 (middle)
                # Column 3 -> BC1 (rightmost)
                bc3_col, bc2_col, bc1_col = parts[0], parts[1], parts[2]
                if bc1_col != '-':
                    bc1_set.add(bc1_col)
                if bc2_col != '-':
                    bc2_set.add(bc2_col)
                if bc3_col != '-':
                    bc3_set.add(bc3_col)
    
    logging.info(f"Loaded whitelists: BC1={len(bc1_set)}, BC2={len(bc2_set)}, BC3={len(bc3_set)}")
    return bc1_set, bc2_set, bc3_set


def validate_spacer(sequence, barcode_pos, expected_spacer_end, barcode_name):
    """
    Validate spacer sequence before a barcode.
    Returns True if spacer is valid, False otherwise.
    """
    spacer_end_pos = barcode_pos - 5
    seq_len = len(sequence)
    
    if spacer_end_pos >= 0 and spacer_end_pos + 5 <= seq_len:
        spacer_end = sequence[spacer_end_pos:spacer_end_pos + 5]
        if spacer_end != expected_spacer_end:
            logging.debug(f"{barcode_name} spacer mismatch: expected '{expected_spacer_end}', got '{spacer_end}'")
            return False
    return True


def find_barcodes_simple(sequence, bc1_set, bc2_set, bc3_set, max_shift=4):
    """
    Find shifted barcodes by checking BC1 and BC2 at various positions.
    BC3 is always taken from position 10-18 without checking.
    Always validates spacer sequences to confirm correct barcode positions.
    Returns: (shift1, shift2, bc1, bc2, bc3, spacer_valid) or None if no shifts needed.
    """
    seq_len = len(sequence)
    
    # BC3 is always at fixed position 10-18 (no checking needed)
    bc3 = sequence[10:18]
    
    # Fast path: check if both BC2 and BC1 are at expected positions
    bc2_default = sequence[48:56]
    bc1_default = sequence[78:86]
    
    if bc2_default in bc2_set and bc1_default in bc1_set:
        # Both found at expected positions - no shift needed, skip this read
        return None
    
    # Find BC2 position (check in order: 0, -1, 1, -2, 2, -3, 3)
    best_bc2 = bc2_default
    best_shift2 = 0
    
    # Check position 0 first (if not already validated in fast path)
    if bc2_default not in bc2_set:
        # Check other positions in order of increasing shift
        for abs_shift in range(1, max_shift + 1):
            for shift2 in [-abs_shift, abs_shift]:
                bc2_start = 48 + shift2
                bc2_end = 56 + shift2
                if bc2_start >= 0 and bc2_end <= seq_len:
                    bc2 = sequence[bc2_start:bc2_end]
                    if bc2 in bc2_set:
                        # Check spacer validation for this position
                        if validate_spacer(sequence, 48 + shift2, "CGACT", f"BC2 at shift {shift2}"):
                            best_shift2 = shift2
                            best_bc2 = bc2
                            break
            if best_shift2 != 0:
                break
    
    # Find BC1 position (check in order: 0, -1, 1, -2, 2, -3, 3)
    best_bc1 = bc1_default
    best_shift1 = 0
    
    # Check position 0 first (if not already validated in fast path)
    if bc1_default not in bc1_set:
        # Check other positions in order of increasing shift
        for abs_shift in range(1, max_shift + 1):
            for shift1 in [-abs_shift, abs_shift]:
                bc1_start = 78 + shift1
                bc1_end = 86 + shift1
                if bc1_start >= 0 and bc1_end <= seq_len:
                    bc1 = sequence[bc1_start:bc1_end]
                    if bc1 in bc1_set:
                        # Check spacer validation for this position
                        if validate_spacer(sequence, 78 + shift1, "TGTGG", f"BC1 at shift {shift1}"):
                            best_shift1 = shift1
                            best_bc1 = bc1
                            break
            if best_shift1 != 0:
                break
    
    # Only return if at least one barcode needs shifting
    if best_shift1 != 0 or best_shift2 != 0:
        # Need to validate spacers for any barcode at default position
        spacer_valid = True
        
        # If BC2 is at default position (shift=0), validate its spacer
        if best_shift2 == 0 and bc2_default in bc2_set:
            if not validate_spacer(sequence, 48, "CGACT", "BC2 at default position"):
                spacer_valid = False
        
        # If BC1 is at default position (shift=0), validate its spacer  
        if best_shift1 == 0 and bc1_default in bc1_set:
            if not validate_spacer(sequence, 78, "TGTGG", "BC1 at default position"):
                spacer_valid = False
        
        return best_shift1, best_shift2, best_bc1, best_bc2, bc3, spacer_valid
    else:
        return None


def process_reads_simple(r1_file, r2_file, output_prefix, barcode_file, max_reads=None, max_shift=4, require_valid_spacers=True):
    """Process R1/R2 reads using simple whitelist matching."""
    
    logging.info("Using simple barcode recovery with whitelist matching")
    if require_valid_spacers:
        logging.info("Requiring valid spacer sequences for recovery")
    else:
        logging.info("Recovering all reads with valid barcodes (ignoring spacer validation)")
    
    # Load barcode whitelists
    bc1_set, bc2_set, bc3_set = load_barcode_whitelists(barcode_file)
    
    if max_reads:
        logging.info(f"Processing up to {max_reads:,} reads")
    
    # Open input files
    f1 = pysam.FastxFile(r1_file)
    f2 = pysam.FastxFile(r2_file)
    
    # Open output files with pigz compression
    r1_pigz = subprocess.Popen(['pigz', '-c', '-p', '4'], 
                               stdin=subprocess.PIPE,
                               stdout=open(f"{output_prefix}_recovered_R1.fastq.gz", 'wb'),
                               text=True)
    r2_pigz = subprocess.Popen(['pigz', '-c', '-p', '4'],
                               stdin=subprocess.PIPE,
                               stdout=open(f"{output_prefix}_recovered_R2.fastq.gz", 'wb'),
                               text=True)
    
    r1_out = r1_pigz.stdin
    r2_out = r2_pigz.stdin
    
    # Statistics
    stats = {
        'total': 0,
        'found': 0,
        'found_with_valid_spacers': 0,
        'found_with_invalid_spacers': 0,
        'written': 0,  # Actually written to output
        'shift_dist': {},
        'bc1_shifts': {},
        'bc2_shifts': {},
        'spacer_mismatches_by_shift': {}
    }
    
    processed = 0
    last_log_time = time.time()
    start_time = time.time()
    
    # Process reads
    for r1, r2 in zip(f1, f2):
        processed += 1
        stats['total'] += 1
        
        # Find shifted barcodes
        result = find_barcodes_simple(
            r2.sequence, bc1_set, bc2_set, bc3_set, max_shift
        )
        
        if result is not None:
            shift1, shift2, bc1, bc2, bc3, spacer_valid = result
            # Found valid combination
            stats['found'] += 1
            
            # Track spacer validation
            if spacer_valid:
                stats['found_with_valid_spacers'] += 1
            else:
                stats['found_with_invalid_spacers'] += 1
                # Track which shifts have spacer mismatches
                shift_combo = (shift1, shift2)
                stats['spacer_mismatches_by_shift'][shift_combo] = stats['spacer_mismatches_by_shift'].get(shift_combo, 0) + 1
            
            # Track shift combinations
            shift_combo = (shift1, shift2)
            stats['shift_dist'][shift_combo] = stats['shift_dist'].get(shift_combo, 0) + 1
            
            # Track individual barcode shifts
            stats['bc1_shifts'][shift1] = stats['bc1_shifts'].get(shift1, 0) + 1
            stats['bc2_shifts'][shift2] = stats['bc2_shifts'].get(shift2, 0) + 1
            
            # Only write reads if spacers are valid (or if we're not requiring valid spacers)
            if spacer_valid or not require_valid_spacers:
                stats['written'] += 1
                
                # Extract UMI (always at fixed position)
                umi = r2.sequence[0:10]
                
                # Build standardized R2 sequence
                new_r2_seq = umi + bc3 + SPACER1 + bc2 + SPACER2 + bc1
                
                # R2 must be at least 86bp - fail if shorter
                if len(r2.quality) < 86:
                    logging.error(f"R2 read {r2.name} has quality string length {len(r2.quality)}, expected at least 86bp")
                    sys.exit(1)
                
                new_r2_qual = r2.quality[:86]
                
                # Write recovered reads
                r2_out.write(f"@{r2.name}\n{new_r2_seq}\n+\n{new_r2_qual}\n")
                r1_out.write(f"@{r1.name}\n{r1.sequence}\n+\n{r1.quality}\n")
        
        # Progress report every 5 seconds
        current_time = time.time()
        if current_time - last_log_time >= 5:
            elapsed = current_time - start_time
            rate = processed / elapsed
            logging.info(f"Processed {processed:,} reads - {rate:.0f} reads/sec")
            last_log_time = current_time
        
        # Check max reads
        if max_reads and processed >= max_reads:
            break
    
    f1.close()
    f2.close()
    
    # Close output files
    r1_out.close()
    r2_out.close()
    r1_pigz.wait()
    r2_pigz.wait()
    
    elapsed = time.time() - start_time
    
    # Print statistics
    print(f"\nCompleted in {elapsed/60:.1f} minutes ({stats['total']/elapsed:.0f} reads/sec)")
    print(f"\nRecovery Statistics:")
    print(f"Total reads: {stats['total']:,}")
    print(f"Found with valid barcodes: {stats['found']:,} ({stats['found']/stats['total']*100:.2f}%)")
    print(f"  - With valid spacers: {stats['found_with_valid_spacers']:,} ({stats['found_with_valid_spacers']/stats['found']*100:.2f}% of found)")
    print(f"  - With invalid spacers: {stats['found_with_invalid_spacers']:,} ({stats['found_with_invalid_spacers']/stats['found']*100:.2f}% of found)")
    print(f"Written to output: {stats['written']:,} ({stats['written']/stats['total']*100:.2f}% of total)")
    
    print(f"\nBC1 shift distribution (of all found reads):")
    for shift in sorted(stats['bc1_shifts'].keys()):
        count = stats['bc1_shifts'][shift]
        print(f"  Shift {shift:+2d}: {count:,} ({count/stats['found']*100:.2f}% of found)")
    
    print(f"\nBC2 shift distribution (of all found reads):")
    for shift in sorted(stats['bc2_shifts'].keys()):
        count = stats['bc2_shifts'][shift]
        print(f"  Shift {shift:+2d}: {count:,} ({count/stats['found']*100:.2f}% of found)")
    
    print(f"\nShift combinations (BC1, BC2) - of all found reads:")
    for shift_combo in sorted(stats['shift_dist'].keys()):
        count = stats['shift_dist'][shift_combo]
        print(f"  Shifts {shift_combo}: {count:,} ({count/stats['total']*100:.2f}% of total)")
    
    if stats['spacer_mismatches_by_shift']:
        print(f"\nSpacers mismatches by shift combination:")
        for shift_combo in sorted(stats['spacer_mismatches_by_shift'].keys()):
            count = stats['spacer_mismatches_by_shift'][shift_combo]
            print(f"  Shifts {shift_combo}: {count:,} mismatches")
    
    # Write stats to file
    with open(f"{output_prefix}_stats.txt", 'w') as f:
        f.write(f"Total reads\t{stats['total']}\n")
        f.write(f"Found with valid barcodes\t{stats['found']}\n")
        f.write(f"Found rate\t{stats['found']/stats['total']*100:.2f}%\n")
        f.write(f"Valid spacers\t{stats['found_with_valid_spacers']}\n")
        f.write(f"Invalid spacers\t{stats['found_with_invalid_spacers']}\n")
        f.write(f"Written to output\t{stats['written']}\n")
        f.write(f"Output rate\t{stats['written']/stats['total']*100:.2f}%\n")
        f.write(f"\nBC1 shifts:\n")
        f.write(f"Shift\tCount\tPercentage\n")
        for shift in sorted(stats['bc1_shifts'].keys()):
            count = stats['bc1_shifts'][shift]
            f.write(f"{shift}\t{count}\t{count/stats['found']*100:.2f}% of found\n")
        
        f.write(f"\nBC2 shifts:\n")
        f.write(f"Shift\tCount\tPercentage\n")
        for shift in sorted(stats['bc2_shifts'].keys()):
            count = stats['bc2_shifts'][shift]
            f.write(f"{shift}\t{count}\t{count/stats['found']*100:.2f}% of found\n")
        
        f.write(f"\nShift combinations:\n")
        f.write(f"BC1_shift\tBC2_shift\tCount\tPercentage\n")
        for shift_combo in sorted(stats['shift_dist'].keys()):
            count = stats['shift_dist'][shift_combo]
            f.write(f"{shift_combo[0]}\t{shift_combo[1]}\t{count}\t{count/stats['total']*100:.2f}% of total\n")
        
        if stats['spacer_mismatches_by_shift']:
            f.write(f"\nSpacers mismatches by shift:\n")
            f.write(f"BC1_shift\tBC2_shift\tMismatches\n")
            for shift_combo in sorted(stats['spacer_mismatches_by_shift'].keys()):
                count = stats['spacer_mismatches_by_shift'][shift_combo]
                f.write(f"{shift_combo[0]}\t{shift_combo[1]}\t{count}\n")
    
    return stats


def run_self_tests():
    """Run self-tests to verify barcode recovery logic."""
    print("\n" + "="*60)
    print("RUNNING SELF-TESTS FOR BARCODE RECOVERY")
    print("="*60)
    
    # Create test whitelists
    # Using simple patterns to make it easy to identify: AAAA1111, BBBB2222, etc.
    bc1_set = {"AAAA1111", "BBBB1111", "CCCC1111", "DDDD1111"}
    bc2_set = {"AAAA2222", "BBBB2222", "CCCC2222", "DDDD2222"}  
    bc3_set = {"AAAA3333", "BBBB3333", "CCCC3333", "DDDD3333"}
    
    # Test cases with visual structure:
    # Format: UMI(10) + BC3(8) + SPACER1(30) + BC2(8) + SPACER2(22) + BC1(8)
    #         0-9      10-17    18-47           48-55    56-77          78-85
    
    test_cases = [
        {
            "name": "Test 1: All barcodes at correct positions",
            "sequence": (
                "UUUUUUUUUU" +  # UMI (10bp)
                "AAAA3333" +    # BC3 (8bp)
                "GTGGCCGATGTTTCGCATCGGCGTACGACT" +  # SPACER1 (30bp)
                "AAAA2222" +    # BC2 (8bp)
                "ATCCACGTGCTTGAGACTGTGG" +  # SPACER2 (22bp)
                "AAAA1111"      # BC1 (8bp)
            ),
            "expected": None,  # Should return None (no recovery needed)
            "description": "Barcodes at expected positions - should skip"
        },
        {
            "name": "Test 2: BC2 shifted +2 (insertion within SPACER1)",
            "sequence": (
                "UUUUUUUUUU" +  # UMI (10bp)
                "BBBB3333" +    # BC3 (8bp)
                "GTGGCCGATGTTTCGCATCGGCGTAXXCGACT" +  # SPACER1 with 2bp insertion before last 5bp
                "BBBB2222" +    # BC2 (8bp) - now at position 50-57
                "ATCCACGTGCTTGAGACTGTGG" +  # SPACER2 (22bp)
                "BBBB1111"      # BC1 (8bp) - now at position 80-87
            ),
            "expected": (2, 2, "BBBB1111", "BBBB2222", "BBBB3333", True),
            "description": "2bp insertion within SPACER1 (preserving last 5bp) shifts BC2 and BC1"
        },
        {
            "name": "Test 3: Invalid BC2 with valid BC1",
            "sequence": (
                "UUUUUUUUUU" +  # UMI (10bp)
                "CCCC3333" +    # BC3 (8bp)
                "GTGGCCGATGTTTCGCATCGGCGTACGACT" +  # SPACER1 (30bp)
                "XXXXXXXX" +    # Invalid BC2 (8bp)
                "ATCCACGTGCTTGAGACTGTGG" +  # SPACER2 (22bp)
                "CCCC1111"      # BC1 (8bp)
            ),
            "expected": None,  # No shifts needed, so returns None
            "description": "Invalid BC2 but no shifts - not a recovery case"
        },
        {
            "name": "Test 4: Indel between spacer and barcode - should fail",
            "sequence": (
                "UUUUUUUUUU" +  # UMI (10bp)
                "DDDD3333" +    # BC3 (8bp)
                "GTGGCCGATGTTTCGCATCGGCGTACGACT" +  # SPACER1 (30bp)
                "X" +           # Insertion between spacer and BC2
                "DDDD2222" +    # BC2 (8bp) at position 49
                "ATCCACGTGCTTGAGACTGTGG" +  # SPACER2 (22bp)
                "X" +           # Insertion between spacer and BC1
                "DDDD1111"      # BC1 (8bp) at position 80
            ),
            "expected": None,  # Should return None (spacers not adjacent to barcodes)
            "description": "Insertions between spacers and barcodes break spacer validation"
        },
        {
            "name": "Test 5: BC1 shifted -1 with BC2 substitution error",
            "sequence": (
                "UUUUUUUUUU" +  # UMI (10bp)
                "AAAA3333" +    # BC3 (8bp)
                "GTGGCCGATGTTTCGCATCGGCGTACGACT" +  # SPACER1 (30bp)
                "XAAA2222" +    # BC2 with substitution error (X instead of A, making it invalid)
                "ATCCACGTGCTTGAGATGTGG" +  # SPACER2 with 1bp deletion (removed C from GACT, now 21bp, keeping TGTGG at end)
                "AAAA1111"      # BC1 (8bp) - now at position 77-84 due to deletion
            ),
            "expected": (-1, 0, "AAAA1111", "XAAA2222", "AAAA3333", True),
            "description": "BC1 shifted left by deletion in SPACER2, BC2 has substitution error but no shift"
        }
    ]
    
    # Run tests
    passed = 0
    failed = 0
    
    for test in test_cases:
        print(f"\n{test['name']}")
        print(f"  Description: {test['description']}")
        print(f"  Sequence length: {len(test['sequence'])}")
        
        result = find_barcodes_simple(test['sequence'], bc1_set, bc2_set, bc3_set, max_shift=4)
        
        if result == test['expected']:
            print(f"  ✓ PASSED - Got expected result: {result}")
            passed += 1
        else:
            print(f"  ✗ FAILED")
            print(f"    Expected: {test['expected']}")
            print(f"    Got:      {result}")
            failed += 1
            
            # Debug info for failures
            if result:
                shift1, shift2, bc1, bc2, bc3, spacer_valid = result
                print(f"    Details: BC1='{bc1}' (shift={shift1}), BC2='{bc2}' (shift={shift2}), BC3='{bc3}', spacers_valid={spacer_valid}")
    
    # Summary
    print("\n" + "="*60)
    print(f"TEST SUMMARY: {passed} passed, {failed} failed out of {len(test_cases)} tests")
    if failed > 0:
        print("ERROR: Some tests failed! Please check the barcode recovery logic.")
        sys.exit(1)
    else:
        print("SUCCESS: All tests passed!")
    print("="*60 + "\n")


def main():
    parser = argparse.ArgumentParser(description='Simple barcode recovery using whitelist matching')
    parser.add_argument('r1', help='R1 FASTQ.gz file', nargs='?', default=None)
    parser.add_argument('r2', help='R2 FASTQ.gz file', nargs='?', default=None)
    parser.add_argument('output_prefix', help='Output prefix for recovered R1/R2 files', nargs='?', default=None)
    parser.add_argument('--barcode-file', required=False,
                        help='Barcode whitelist file')
    parser.add_argument('--chemistry', required=False,
                        help='kb_chemistry string from config')
    parser.add_argument('--max-reads', type=int, help='Maximum number of reads to process')
    parser.add_argument('--max-shift', type=int, default=4, 
                        help='Maximum shift to search (default: 4)')
    parser.add_argument('--allow-invalid-spacers', action='store_true',
                        help='Also recover reads with invalid spacer sequences (default: require valid spacers)')
    parser.add_argument('--log-level', default='INFO', help='Logging level')
    parser.add_argument('--test', action='store_true', help='Run self-tests and exit')
    
    args = parser.parse_args()
    
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Always run self-tests first to verify logic is correct
    print("Running self-tests before processing...")
    run_self_tests()
    print("Self-tests completed successfully.\n")
    
    # If only testing, exit now
    if args.test:
        print("Test mode - exiting without processing files.")
        sys.exit(0)
    
    # Check if required arguments are provided for actual processing
    if args.r1 is None:
        print("Error: No input files provided. Run with --help for usage or --test to run self-tests.")
        sys.exit(1)
    
    # Verify required arguments for actual processing
    if not args.barcode_file or not args.chemistry:
        logging.error("ERROR: --barcode-file and --chemistry are required for processing")
        sys.exit(1)
    
    # Verify chemistry is v2
    expected_chemistry = "1,10,18,1,48,56,1,78,86:1,0,10:0,0,0"
    if args.chemistry != expected_chemistry:
        logging.error(f"ERROR: This script only supports v2 chemistry")
        logging.error(f"Expected: {expected_chemistry}")
        logging.error(f"Got: {args.chemistry}")
        sys.exit(1)
    
    logging.info(f"Verified v2 chemistry: {args.chemistry}")
    
    try:
        import pysam
    except ImportError:
        logging.error("pysam not installed. Install with: pip install pysam")
        sys.exit(1)
    
    process_reads_simple(
        args.r1, args.r2, args.output_prefix, args.barcode_file,
        args.max_reads, args.max_shift, 
        require_valid_spacers=not args.allow_invalid_spacers
    )


if __name__ == '__main__':
    main()