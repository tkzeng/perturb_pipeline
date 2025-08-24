#!/usr/bin/env python3
"""
Wrapper script for kb count to reduce duplication in Snakefile.
Handles recovered reads concatenation and provides sensible defaults.
"""

import argparse
import subprocess
import os
import sys
import logging

def run_kb_count(args):
    """Run kb count with appropriate parameters."""
    
    # Build the basic command with all the common parameters
    cmd = [
        "kb", "count",
        "-i", args.index,
        "-g", args.t2g,
        "-o", args.output,
        "-x", args.tech,
        "-t", str(args.threads),
        "--workflow", args.workflow,
        "--strand", args.strand,
        "--overwrite",
        "--h5ad",
        "--sum", "total",
        "--verbose",
        "-w", args.barcodes,
        "-r", args.replace
    ]
    
    # Add workflow-specific parameters
    if args.workflow == "nac":
        if args.cdna and args.nascent:
            cmd.extend(["-c1", args.cdna, "-c2", args.nascent])
        if args.mm:
            cmd.append("--mm")
    
    # Add subsampling if specified
    if args.max_reads:
        cmd.extend(["-N", str(args.max_reads)])
    
    # Add the input files
    cmd.extend([args.r1, args.r2])
    
    # Execute the command
    logging.info(f"Running kb count for {args.output}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"kb count failed with exit code {e.returncode}")
        print(e.stdout)
        print(e.stderr, file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(description="Wrapper for kb count with sensible defaults")
    
    # Positional arguments that always change
    parser.add_argument("workflow", choices=["nac", "kite"], help="Workflow type")
    parser.add_argument("strand", choices=["unstranded", "forward"], help="Strandedness")
    parser.add_argument("r1", help="R1 FASTQ file")
    parser.add_argument("r2", help="R2 FASTQ file") 
    parser.add_argument("output", help="Output directory")
    parser.add_argument("threads", type=int, help="Number of threads")
    
    # Common file paths from Snakemake
    parser.add_argument("index", help="Index file")
    parser.add_argument("t2g", help="T2G file")
    parser.add_argument("barcodes", help="Barcodes file")
    parser.add_argument("replace", help="Replace file")
    parser.add_argument("tech", help="Technology string for -x parameter")
    
    # Optional workflow-specific arguments
    parser.add_argument("--cdna", help="cDNA file (for nac workflow)")
    parser.add_argument("--nascent", help="Nascent file (for nac workflow)")
    parser.add_argument("--mm", action="store_true", help="Use --mm flag for nac")
    parser.add_argument("--max-reads", type=int, help="Maximum reads for subsampling")
    
    args = parser.parse_args()
    
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(message)s'
    )
    
    # Run kb count
    success = run_kb_count(args)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()