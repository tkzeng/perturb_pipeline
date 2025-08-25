#!/usr/bin/env bash

# Submit FASTQ preparation pipeline to SLURM
# Usage: ./submit_prepare_fastq.sh [snakemake args]
#
# This script runs prepare_fastq.smk which handles:
#   - Lane concatenation for multi-lane samples
#   - Undetermined read recovery and per-sample splitting
#   - Barcode recovery on main and undetermined reads
#   - Cross-pool concatenation for multi-pool samples
#   - Final merging of all FASTQ sources
#
# Examples:
#   CONFIG=config.gw_pilot.yaml ./submit_prepare_fastq.sh --dry-run
#   CONFIG=config.gw_pilot.yaml ./submit_prepare_fastq.sh --jobs 50

# Require config file via environment variable
if [ -z "$CONFIG" ]; then
    echo "Error: CONFIG environment variable must be set"
    echo "Usage: CONFIG=config.gw_pilot.yaml ./submit_prepare_fastq.sh [args]"
    exit 1
fi
CONFIG_FILE="$CONFIG"

# Extract SLURM settings from config file
SLURM_PARTITION=$(python -c "import yaml; config = yaml.safe_load(open('$CONFIG_FILE')); print(config['slurm']['partition'])")
SLURM_ACCOUNT=$(python -c "import yaml; config = yaml.safe_load(open('$CONFIG_FILE')); print(config['slurm']['account'])")


# Default arguments that can be overridden
# Resource defaults:
#   - threads: 4 (good default for most rules, specific rules override as needed)
#   - mem_mb: 32768 (32GB, sufficient for most operations)
#   - runtime: 2880 (48 hours in minutes)
DEFAULT_ARGS=(
    --executor slurm
    --keep-going
    --retries 0
    --keep-incomplete
    --rerun-incomplete
    --printshellcmds
    --rerun-triggers mtime
    --max-jobs-per-timespan 50/60s
    --default-resources \
        threads=4 \
        mem_mb=32768 \
        runtime=2880 \
        slurm_partition="$SLURM_PARTITION" \
        slurm_account="$SLURM_ACCOUNT" \
        slurm_job_name="{rule}"
    # List of allowed rules - comment out any rules you want to skip
    --allowed-rules \
        all \
        # STAGE 1: Lane concatenation, read counting, and QC \
        concatenate_lanes \
        count_reads \
        fastp_qc_lane \
        # STAGE 2: Undetermined recovery (parallelized per lane) \
        recover_undetermined_per_lane \
        create_undetermined_fastq \
        # STAGE 3: Barcode recovery (parallelized per lane) \
        recover_barcodes_main_per_lane \
        merge_barcode_recovery_main \
        recover_barcodes_undetermined \
        # STAGE 4: Final merging and pool statistics \
        merge_all_fastqs \
        calculate_pool_statistics
)

# If user didn't specify --jobs, add default
if [[ ! " $@ " =~ " --jobs " ]] && [[ ! " $@ " =~ " -j " ]]; then
    DEFAULT_ARGS+=(--jobs 100)
fi

# Add config file and snakefile arguments
DEFAULT_ARGS+=(--configfile "$CONFIG_FILE")
DEFAULT_ARGS+=(-s prepare_fastq.smk)

echo "Using config file: $CONFIG_FILE"
echo "Running FASTQ preparation pipeline: prepare_fastq.smk"
echo "EXECUTING: snakemake ${DEFAULT_ARGS[@]} $@"

snakemake "${DEFAULT_ARGS[@]}" "$@"
