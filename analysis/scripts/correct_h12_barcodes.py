#!/usr/bin/env python3
"""
Correct H12 BC1 barcode collision in Parse Biosciences R2 reads.

The random hexamer (R) plate BC1 barcodes are shifted by +1 well relative to
the polyT (T) plate. The replace file handles wells A1-H11, but the rotation
wraps around: physical R well A1 received ACATTTGG (the T plate H12 barcode),
creating a collision with genuine polyT H12 reads. Under the rotated mapping,
R well A1 should map to CCTCTGTG (R plate H12 barcode, never physically used).

This script:
1. Extracts BC1 from position 78:86 in R2 reads
2. ED-corrects BC1 against the whitelist (edit distance <= 1, substitutions only)
3. For reads where corrected BC1 == ACATTTGG, checks polyT signature
4. Replaces BC1 with CCTCTGTG for random hexamer reads (no polyT)
5. Outputs all reads (1:1 with R1) with corrected BC1 where applicable
"""

import argparse
import logging
import subprocess
import sys
import time

import pysam

# H12 collision constants
H12_POLYT_BC1 = "ACATTTGG"  # polyT plate H12 (also physically on R plate A1 due to shift)
H12_HEXAMER_BC1 = "CCTCTGTG"  # correct rotated barcode for R plate A1 (R-H12, never physically used)

# PolyT detection parameters
BC1_START = 78
BC1_END = 86
POLYT_LENGTH = 15
POLYT_MIN_TS = 11

BASES = "ACGT"


def load_bc1_whitelist(barcode_file):
    """Load BC1 whitelist from barcodes file (column 3)."""
    bc1_set = set()
    with open(barcode_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                bc1 = parts[2]
                if bc1 != '-':
                    bc1_set.add(bc1)
    return bc1_set


def build_ed_lookup(bc1_set):
    """Build dictionary mapping 1-substitution variants to their corrected BC1.

    For each BC1 in the whitelist, generate all single-substitution variants.
    If a variant maps to multiple different BC1s, mark as ambiguous (None).

    Returns dict: variant_seq -> (corrected_bc1, edit_distance) or None if ambiguous
    """
    lookup = {}
    for bc in bc1_set:
        lookup[bc] = (bc, 0)

    for bc in bc1_set:
        for i in range(len(bc)):
            for base in BASES:
                if base == bc[i]:
                    continue
                variant = bc[:i] + base + bc[i + 1:]
                if variant in lookup:
                    existing = lookup[variant]
                    if existing is None:
                        continue
                    existing_bc, existing_ed = existing
                    if existing_ed == 0:
                        continue
                    if existing_bc != bc:
                        lookup[variant] = None
                else:
                    lookup[variant] = (bc, 1)

    n_exact = sum(1 for v in lookup.values() if v is not None and v[1] == 0)
    n_ed1 = sum(1 for v in lookup.values() if v is not None and v[1] == 1)
    n_ambig = sum(1 for v in lookup.values() if v is None)
    logging.info(f"ED lookup: {n_exact} exact, {n_ed1} ed1, {n_ambig} ambiguous")

    return lookup


def is_polyt_read(r2_sequence, bc1_end_pos):
    """Check if read has polyT signature after BC1 position.

    Returns True (polyT), False (random hexamer), or None (read too short to determine).
    """
    primer_start = bc1_end_pos
    primer_end = primer_start + POLYT_LENGTH
    if len(r2_sequence) < primer_end:
        return None  # Read too short to check primer region
    primer_region = r2_sequence[primer_start:primer_end].upper()
    return primer_region.count('T') >= POLYT_MIN_TS


def process_r2(r2_file, barcode_file, output_path, threads=4, max_reads=None):
    """Process R2 FASTQ, correcting H12 BC1 collision."""

    bc1_set = load_bc1_whitelist(barcode_file)
    logging.info(f"Loaded {len(bc1_set)} BC1 barcodes")

    lookup = build_ed_lookup(bc1_set)

    f2 = pysam.FastxFile(r2_file)

    pigz = subprocess.Popen(
        ['pigz', '-c', '-p', str(threads)],
        stdin=subprocess.PIPE,
        stdout=open(output_path, 'wb'),
        text=True
    )
    out = pigz.stdin

    stats = {
        'total': 0,
        'bc1_matched': 0,
        'bc1_is_h12': 0,
        'h12_polyt': 0,
        'h12_hexamer_replaced': 0,
        'h12_too_short': 0,
        'bc1_ambiguous': 0,
        'bc1_no_match': 0,
        'too_short': 0,
        'ed_dist': {0: 0, 1: 0},
    }

    start_time = time.time()
    last_log_time = start_time

    for r2 in f2:
        stats['total'] += 1
        seq = r2.sequence
        qual = r2.quality

        if len(seq) < BC1_END:
            stats['too_short'] += 1
            out.write(f"@{r2.name}\n{seq}\n+\n{qual}\n")
            continue

        raw_bc1 = seq[BC1_START:BC1_END]
        result = lookup.get(raw_bc1)

        if result is None:
            if raw_bc1 in lookup:
                stats['bc1_ambiguous'] += 1
            else:
                stats['bc1_no_match'] += 1
            out.write(f"@{r2.name}\n{seq}\n+\n{qual}\n")
            continue

        corrected_bc1, ed = result
        stats['bc1_matched'] += 1
        stats['ed_dist'][ed] += 1

        # Check if this is the H12 collision barcode
        if corrected_bc1 == H12_POLYT_BC1:
            stats['bc1_is_h12'] += 1

            polyt_result = is_polyt_read(seq, BC1_END)
            if polyt_result is None:
                stats['h12_too_short'] += 1
                # Read too short to determine — leave unchanged
                if ed > 0:
                    seq = seq[:BC1_START] + corrected_bc1 + seq[BC1_END:]
            elif polyt_result:
                stats['h12_polyt'] += 1
                # PolyT read — keep as ACATTTGG (may need ED correction though)
                if ed > 0:
                    seq = seq[:BC1_START] + corrected_bc1 + seq[BC1_END:]
            else:
                stats['h12_hexamer_replaced'] += 1
                # Random hexamer read — replace with correct A1 barcode
                seq = seq[:BC1_START] + H12_HEXAMER_BC1 + seq[BC1_END:]
        elif ed > 0:
            # Non-H12 barcode that needed ED correction — apply correction
            seq = seq[:BC1_START] + corrected_bc1 + seq[BC1_END:]

        out.write(f"@{r2.name}\n{seq}\n+\n{qual}\n")

        current_time = time.time()
        if current_time - last_log_time >= 5:
            elapsed = current_time - start_time
            rate = stats['total'] / elapsed
            logging.info(f"Processed {stats['total']:,} reads — {rate:.0f} reads/sec")
            last_log_time = current_time

        if max_reads and stats['total'] >= max_reads:
            break

    f2.close()
    out.close()
    pigz.wait()

    elapsed = time.time() - start_time

    print(f"\nCompleted in {elapsed / 60:.1f} minutes ({stats['total'] / elapsed:.0f} reads/sec)")
    print(f"\nH12 Correction Statistics:")
    print(f"Total reads: {stats['total']:,}")
    print(f"Too short (<{BC1_END}bp): {stats['too_short']:,}")
    print(f"BC1 matched: {stats['bc1_matched']:,} ({stats['bc1_matched'] / stats['total'] * 100:.2f}%)")
    print(f"  ED=0: {stats['ed_dist'][0]:,}")
    print(f"  ED=1: {stats['ed_dist'][1]:,}")
    print(f"BC1 ambiguous: {stats['bc1_ambiguous']:,}")
    print(f"BC1 no match: {stats['bc1_no_match']:,}")
    print(f"BC1 == ACATTTGG (H12): {stats['bc1_is_h12']:,}")
    print(f"  PolyT (kept as H12): {stats['h12_polyt']:,}")
    print(f"  Random hexamer (replaced → CCTCTGTG): {stats['h12_hexamer_replaced']:,}")
    if stats['h12_too_short'] > 0:
        print(f"  Too short for polyT detection: {stats['h12_too_short']:,}")
        print(f"  WARNING: R2 reads are too short (<{BC1_END + POLYT_LENGTH}bp) for polyT detection.")
        print(f"  These reads were left unchanged.")

    return stats


def main():
    parser = argparse.ArgumentParser(
        description='Correct H12 BC1 barcode collision in R2 FASTQ'
    )
    parser.add_argument('r2', help='R2 FASTQ.gz file')
    parser.add_argument('--barcode-file', required=True, help='Barcode whitelist file')
    parser.add_argument('--output', '-o', required=True, help='Output corrected R2 FASTQ.gz')
    parser.add_argument('--threads', type=int, default=4, help='Threads for pigz (default: 4)')
    parser.add_argument('--max-reads', type=int, help='Max reads to process')
    parser.add_argument('--log-level', default='INFO', help='Logging level')

    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    process_r2(
        args.r2, args.barcode_file, args.output,
        threads=args.threads,
        max_reads=args.max_reads
    )


if __name__ == '__main__':
    main()
