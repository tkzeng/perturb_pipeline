#!/usr/bin/env bash

# Submit Snakemake pipeline to SLURM
# Usage: ./submit.sh -s <snakefile> [snakemake args]
#
# Available Snakefiles:
#   preprocessing.smk - Initial processing, QC, and cell calling
#   downstream.smk    - Downstream analysis (UMAP, clustering, DE)
#
# Examples:
#   CONFIG=config.yann_k562.yaml ./submit.sh -s preprocessing.smk --dry-run
#   CONFIG=config.yann_k562.yaml ./submit.sh -s downstream.smk --jobs 50

# Require config file via environment variable
if [ -z "$CONFIG" ]; then
    echo "Error: CONFIG environment variable must be set"
    echo "Usage: CONFIG=config.yann_k562.yaml ./submit.sh -s <snakefile> [args]"
    exit 1
fi
CONFIG_FILE="$CONFIG"

# Check if -s is specified
if [[ ! " $@ " =~ " -s " ]]; then
    echo "Error: Must specify Snakefile with -s"
    echo "Available Snakefiles:"
    echo "  -s preprocessing.smk  # For initial processing and QC"
    echo "  -s downstream.smk     # For downstream analysis"
    echo ""
    echo "Example: CONFIG=$CONFIG_FILE ./submit.sh -s preprocessing.smk --dry-run"
    exit 1
fi

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
        # STAGE 1: Input processing and counting \
        #count_reads \
        #calculate_pool_statistics \
        # STAGE 2: Undetermined recovery and barcode correction \
        #check_undetermined_barcodes \
        #create_undetermined_fastq \
        #barcode_recovery \
        #merge_all_fastqs \
        # STAGE 3: Reference preparation \
        #kite_index \
        #generate_combined_whitelist \
        # STAGE 4: Alignment \
        #kallisto_gex \
        #kallisto_guide \
        #kallisto_gex_subsampled \
        #kallisto_guide_subsampled \
        # STAGE 5: Post-alignment processing \
        #inspect_bus_files \
        #filter_and_annotate_sublibrary \
        #calculate_read_statistics \
        # STAGE 6: QC and analysis \
        #fastp_qc \
        cell_calling_analysis \
        cell_calling_plots \
        generate_qc_cell_lists \
        calculate_qc_metrics_stratified \
        #umi_saturation_analysis \
        #umi_saturation_analysis_guide \
        # STAGE 7: Consolidation and visualization \
        consolidate_qc_metrics \
        visualize_consolidated_qc \
        #process_pool_metrics \
        # STAGE 8: Final outputs \
        generate_qc_report \
        # DOWNSTREAM PIPELINE (from downstream.smk) \
        # Note: To run these, use: snakemake -s downstream.smk ... \
        combine_sublibraries \
        standard_analyses \
        plot_standard_analyses_umap \
        generate_final_report
)

# If user didn't specify --jobs, add default
if [[ ! " $@ " =~ " --jobs " ]] && [[ ! " $@ " =~ " -j " ]]; then
    DEFAULT_ARGS+=(--jobs 100)
fi

# Add config file argument
DEFAULT_ARGS+=(--configfile "$CONFIG_FILE")

echo "Using config file: $CONFIG_FILE"
echo "EXECUTING: snakemake ${DEFAULT_ARGS[@]} $@"

snakemake "${DEFAULT_ARGS[@]}" "$@"
