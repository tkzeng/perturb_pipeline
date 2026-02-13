#!/usr/bin/env bash

# Submit Snakemake pipeline to SLURM
# Usage: ./submit.sh [snakemake args]
# Example: ./submit.sh --dry-run
#          ./submit.sh --jobs 50
#          ./submit.sh --dry-run --rerun-triggers mtime
#          CONFIG=config.yann_k562.yaml ./submit.sh --dry-run

# Require config file via environment variable
if [ -z "$CONFIG" ]; then
    echo "Error: CONFIG environment variable must be set"
    echo "Usage: CONFIG=config.yann_k562.yaml ./submit.sh [args]"
    exit 1
fi
CONFIG_FILE="$CONFIG"

# Extract SLURM settings from config file
SLURM_PARTITION=$(python -c "import yaml; config = yaml.safe_load(open('$CONFIG_FILE')); print(config['slurm']['partition'])")
SLURM_ACCOUNT=$(python -c "import yaml; config = yaml.safe_load(open('$CONFIG_FILE')); print(config['slurm']['account'])")
analysis_name=$(python -c "import yaml; config = yaml.safe_load(open('$CONFIG_FILE')); print(config['analysis_name'])")

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
        count_reads \
        calculate_pool_statistics \
        # STAGE 2: Undetermined recovery and barcode correction \
        #check_undetermined_barcodes \
        #create_undetermined_fastq \
        #barcode_recovery \
        #merge_all_fastqs \
        # STAGE 3: Reference preparation \
        kite_index \
        generate_combined_whitelist \
        # STAGE 4: Alignment \
        kallisto_gex \
        kallisto_guide \
        kallisto_gex_subsampled \
        kallisto_guide_subsampled \
        # STAGE 5: Post-alignment processing \
        inspect_bus_files \
        filter_and_annotate_sublibrary \
        calculate_read_statistics \
        # STAGE 6: QC and analysis \
        fastp_qc \
        cell_calling_analysis \
        cell_calling_plots \
        calculate_qc_metrics_stratified \
        umi_saturation_analysis \
        umi_saturation_analysis_guide \
        # STAGE 7: Consolidation and visualization \
        consolidate_qc_metrics \
        visualize_consolidated_qc \
        process_pool_metrics \
        # STAGE 8: Final outputs \
        generate_qc_report \
        #combine_sublibraries
)

# If user didn't specify --jobs, add default
if [[ ! " $@ " =~ " --jobs " ]] && [[ ! " $@ " =~ " -j " ]]; then
    DEFAULT_ARGS+=(--jobs 100)
fi

# Add config file argument
DEFAULT_ARGS+=(--configfile "$CONFIG_FILE")

echo "Using config file: $CONFIG_FILE"
echo "Analysis name: ${analysis_name}"
echo "EXECUTING: snakemake ${DEFAULT_ARGS[@]} $@"

snakemake "${DEFAULT_ARGS[@]}" "$@"
