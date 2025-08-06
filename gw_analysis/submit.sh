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
    --max-jobs-per-timespan 50/60s
    --default-resources \
        threads=4 \
        mem_mb=32768 \
        runtime=2880 \
        slurm_partition="pritch,engreitz,hns" \
        slurm_account="pritch" \
        slurm_job_name="{rule}"
    # List of allowed rules - comment out any rules you want to skip
    --allowed-rules \
        all \
        #count_reads \
        #kallisto_gex \
        #kallisto_guide \
        #kallisto_gex_subsampled \
        #kallisto_guide_subsampled \
        #generate_combined_whitelist \
        #inspect_bus_files \
        #merge_all_fastqs \
        #filter_and_annotate_sublibrary \
        #fastp_qc \
        #calculate_read_statistics \
        #cell_calling_analysis \
        #cell_calling_plots \
        calculate_qc_metrics_stratified \
        #umi_saturation_analysis \
        #umi_saturation_analysis_guide \
        #check_undetermined_barcodes \
        #create_undetermined_fastq \
        #barcode_recovery \
        consolidate_qc_metrics \
        visualize_consolidated_qc \
        process_pool_metrics \
        generate_qc_report \
        #calculate_pool_statistics
)

# If user didn't specify --jobs, add default
if [[ ! " $@ " =~ " --jobs " ]] && [[ ! " $@ " =~ " -j " ]]; then
    DEFAULT_ARGS+=(--jobs 100)
fi

echo "EXECUTING: snakemake ${DEFAULT_ARGS[@]} $@"

snakemake "${DEFAULT_ARGS[@]}" "$@"
