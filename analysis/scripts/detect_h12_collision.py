#!/usr/bin/env python3
"""
Detect whether an R2 FASTQ has the H12 BC1 barcode collision.

Samples reads, checks for exact BC1=ACATTTGG match, and classifies by polyT signature.
Reports whether:
- Only random hexamer reads are present (A1 only) → use replace file
- Only polyT reads are present (H12 only) → no correction needed
- Both are present → use correct_h12_barcodes.py for polyT detection
"""

import argparse
import logging

import pysam

from correct_h12_barcodes import (
    is_polyt_read,
    H12_POLYT_BC1,
    BC1_START,
    BC1_END,
)


def detect_collision(r2_file, sample_size=1_000_000):
    """Sample reads and classify H12 BC1 reads as polyT or random hexamer."""

    f2 = pysam.FastxFile(r2_file)

    n_total = 0
    n_h12 = 0
    n_polyt = 0
    n_hexamer = 0
    n_too_short = 0

    for r2 in f2:
        n_total += 1
        seq = r2.sequence

        if len(seq) < BC1_END:
            continue

        raw_bc1 = seq[BC1_START:BC1_END]
        if raw_bc1 != H12_POLYT_BC1:
            continue

        n_h12 += 1
        result = is_polyt_read(seq, BC1_END)
        if result is None:
            n_too_short += 1
        elif result:
            n_polyt += 1
        else:
            n_hexamer += 1

        if n_total >= sample_size:
            break

    f2.close()
    return n_total, n_h12, n_polyt, n_hexamer, n_too_short


def main():
    parser = argparse.ArgumentParser(
        description='Detect H12 BC1 barcode collision in R2 FASTQ'
    )
    parser.add_argument('r2', help='R2 FASTQ.gz file')
    parser.add_argument('--sample-size', type=int, default=1_000_000,
                        help='Number of reads to sample (default: 1M)')
    parser.add_argument('--log-level', default='INFO', help='Logging level')

    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    n_total, n_h12, n_polyt, n_hexamer, n_too_short = detect_collision(
        args.r2,
        sample_size=args.sample_size,
    )

    print(f"\nH12 Collision Detection Report")
    print(f"{'=' * 50}")
    print(f"Reads sampled: {n_total:,}")
    print(f"Reads with BC1=ACATTTGG (exact match): {n_h12:,} ({n_h12 / n_total * 100:.3f}%)")

    if n_h12 == 0:
        print(f"\nNo reads with BC1=ACATTTGG found.")
        print(f"Recommendation: No H12 correction needed.")
        return

    if n_too_short > 0:
        print(f"  Too short for polyT detection: {n_too_short:,} ({n_too_short / n_h12 * 100:.1f}% of H12 reads)")

    if n_too_short == n_h12:
        print(f"\nERROR: All R2 reads are too short (<{BC1_END + 15}bp) for polyT detection.")
        print(f"Cannot distinguish polyT H12 from random hexamer R-A1.")
        print(f"If this library uses well A1, use the replace file approach:")
        print(f"  ACATTTGG *CCTCTGTG")
        return

    print(f"  PolyT (genuine H12): {n_polyt:,} ({n_polyt / n_h12 * 100:.1f}% of H12 reads)")
    print(f"  Random hexamer (misassigned R-A1): {n_hexamer:,} ({n_hexamer / n_h12 * 100:.1f}% of H12 reads)")

    print(f"\n{'=' * 50}")
    if n_polyt == 0 and n_hexamer > 0:
        print(f"Recommendation: R-A1 ONLY — add this line to your replace file:")
        print(f"  ACATTTGG *CCTCTGTG")
    elif n_polyt > 0 and n_hexamer == 0:
        print(f"Recommendation: H12 ONLY — no correction needed.")
    elif n_polyt > 0 and n_hexamer > 0:
        print(f"Recommendation: BOTH R-A1+H12 — use correct_h12_barcodes.py with polyT detection.")
        print(f"  python scripts/correct_h12_barcodes.py <R2.fastq.gz> \\")
        print(f"      --barcode-file <barcodes_file> \\")
        print(f"      --output <corrected_R2.fastq.gz>")


if __name__ == '__main__':
    main()
