#!/usr/bin/env bash

# Submit Snakemake pipeline to SLURM
# Usage: ./submit.sh [snakemake args]
# Example: ./submit.sh --dry-run
#          ./submit.sh --jobs 50
#          ./submit.sh --dry-run --rerun-triggers mtime

# Ensure we're in the kb conda environment
if [[ "$CONDA_DEFAULT_ENV" != "kb" ]]; then
    echo "Error: Must be in 'kb' conda environment"
    echo "Run: conda activate kb"
    exit 1
fi

# Default arguments that can be overridden
DEFAULT_ARGS=(
    --executor slurm
    --keep-going
    --retries 0
    --keep-incomplete
    --rerun-incomplete
    --printshellcmds
    --max-jobs-per-timespan 50/60s
    --default-resources mem_mb=65536 runtime=2880 slurm_partition="pritch,engreitz,hns" slurm_account="pritch" slurm_job_name="{rule}"
)

# If user didn't specify --jobs, add default
if [[ ! " $@ " =~ " --jobs " ]] && [[ ! " $@ " =~ " -j " ]]; then
    DEFAULT_ARGS+=(--jobs 100)
fi

echo "EXECUTING: snakemake ${DEFAULT_ARGS[@]} $@"

snakemake "${DEFAULT_ARGS[@]}" "$@"