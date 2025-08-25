# PREPROCESSING PIPELINE
# ======================
# This pipeline handles initial data processing, QC, and cell calling.
# 
# NOTE: This pipeline should be run from the 'kb' conda environment
# All tools (kallisto-bustools, fastqc, python packages) are available in this environment
#
# To run:
#   snakemake -s preprocessing.smk --configfile config.yaml
#
# For downstream analysis (UMAP, clustering, DE), use:
#   snakemake -s downstream.smk --configfile config.yaml

# Config file must be specified via --configfile when running snakemake

import os
import sys
from datetime import datetime
from scripts.pipeline_utils import get_guide_gex_pairings
from scripts.snakemake_helpers import (
    # Path helper functions
    get_scratch_path, get_results_path, get_logs_path,
    get_reference_path, get_raw_data_path, print_path_configuration,
    # Sample management functions
    load_sample_info, get_sample_ids, get_samples_by_type,
    extract_pool, extract_sample, get_sample_pool,
    # File finding functions
    find_processed_fastq, print_sample_summary,
    # Output generation functions
    get_preprocessing_outputs
)

# Define base paths using helper functions
SCRATCH = get_scratch_path(config=config)
RESULTS = get_results_path(config=config)
LOGS = get_logs_path(config=config)
REFERENCES = get_reference_path(config=config)

# Print configuration at startup
print_path_configuration(config)

# Create logs directory once at startup
os.makedirs(LOGS, exist_ok=True)

# Load sample information and IDs
IDS = get_sample_ids(config)
GEX_IDS = get_samples_by_type('gex', config)
GUIDE_IDS = get_samples_by_type('guide', config)
# Get guide-GEX pairings from the centralized function
guide_to_gex, gex_to_guide = get_guide_gex_pairings(config['sample_info_file'])

# Print sample summary
print_sample_summary(config)

# Wrapper function for find_processed_fastq to pass SCRATCH
def find_fastq_file_wrapper(sample_id, read, source, processing):
    return find_processed_fastq(sample_id, read, source, processing, config, SCRATCH)

# Wrapper function for load_sample_info  
def load_sample_info_wrapper():
    return load_sample_info(config)

rule all:
    input:
        # The QC report DONE file now depends on ALL pipeline outputs via get_all_outputs()
        # This ensures all plots, metrics, and analyses are complete before packaging
        f"{RESULTS}/preprocessing_report/DONE.txt",
        # Additional outputs not covered by get_all_outputs:
        # Ensure kite index is built
        f"kite_references/{os.path.splitext(os.path.basename(config['guide_file']))[0]}_index/mismatch.idx",
        # Ensure gene annotation cache is built
        f"{REFERENCES}/gene_annotations.tsv" 



# =============================================================================
# STAGE 0: REFERENCE PREPARATION
# =============================================================================
rule generate_gene_annotation_table:
    """Generate pre-computed gene annotation table for fast lookups."""
    input:
        gene_database=config['input_paths']['gene_database'],
        ribosomal_genes=config['input_paths']['ribosomal_genes'],
        cell_cycle_genes=config['input_paths']['cell_cycle_genes'],
        script="scripts/generate_gene_annotation_table.py"
    output:
        table_tsv=f"{REFERENCES}/gene_annotations.tsv"
    log:
        f"{LOGS}/generate_gene_annotation_table.log"
    threads: 1
    resources:
        mem_mb=config["resources"]["light"]["mem_mb"]
    shell:
        """
        python3 {input.script} \
            --gene-database {input.gene_database} \
            --ribosomal-genes {input.ribosomal_genes} \
            --cell-cycle-genes {input.cell_cycle_genes} \
            --output {output.table_tsv} \
            &> {log}
        """


# =============================================================================
# STAGE 1: INPUT PROCESSING AND COUNTING
# =============================================================================








# =============================================================================
# STAGE 2: CROSS-RUN MERGING
# =============================================================================
# NOTE: Per-run processing (lane-level recovery + per-run concatenation) is handled by per_run_processing.smk
# This stage merges samples that appear in multiple sequencing runs




# =============================================================================
# STAGE 3: REFERENCE PREPARATION
# =============================================================================
rule kite_index:
    """Build kite index from guides file using kb ref"""
    input:
        guides=lambda wildcards: config['guide_file']
    output:
        index_idx=f"kite_references/{os.path.splitext(os.path.basename(config['guide_file']))[0]}_index/mismatch.idx",
        t2g=f"kite_references/{os.path.splitext(os.path.basename(config['guide_file']))[0]}_index/t2g.txt",
        fa=f"kite_references/{os.path.splitext(os.path.basename(config['guide_file']))[0]}_index/mismatch.fa"
    log:
        f"{LOGS}/kite_index.log"
    params:
        outdir=f"kite_references/{os.path.splitext(os.path.basename(config['guide_file']))[0]}_index"
    threads: config["resources"]["alignment"]["threads"]
    resources:
        mem_mb=config["resources"]["alignment"]["mem_mb"],
    shell:
        """
        # Validate guide file format
        python scripts/validate_guide_file.py {input.guides} || exit 1
        
        # Check if the installed version has the modified collision detection
        if ! python -c "import kb_python.ref; import inspect; exit(0 if 'remove_collisions' in dir(kb_python.ref) else 1)" 2>/dev/null; then
            echo "ERROR: The installed kb_python does not have the modified remove_collisions function"
            echo "Please install the modified version from the kb_python directory"
            echo "Run: cd /path/to/modified/kb_python && pip install -e ."
            exit 1
        fi
        
        # Remove existing index files to ensure clean build
        rm -rf {params.outdir}/* tmp
        
        # Build kite index with kb ref
        kb ref -i {output.index_idx} -f1 {output.fa} -g {output.t2g} \
            --workflow kite {input.guides} -t {threads} &> {log}
        """


rule generate_combined_whitelist:
    input:
        barcodes=config["input_paths"]["barcodes_file"],
        script="scripts/generate_barcode_combinations.py"
    output:
        combined=f"{REFERENCES}/barcodes_combined.txt"
    log:
        f"{LOGS}/generate_combined_whitelist.log"
    shell:
        """
        python {input.script} {input.barcodes} {output.combined} &> {log}
        """



# =============================================================================
# STAGE 4: ALIGNMENT (KALLISTO/BUSTOOLS)
# =============================================================================
rule kallisto_gex:
    input:
        fq1=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R1", wildcards.source, wildcards.processing),
        fq2=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R2", wildcards.source, wildcards.processing),
        replace=config["input_paths"]["replace_file"],
        barcodes=config["input_paths"]["barcodes_file"],
        index_idx=f"{REFERENCES}/nascent_all/index.idx",
        t2g=f"{REFERENCES}/nascent_all/t2g.txt",
        cdna=f"{REFERENCES}/nascent_all/cdna.txt",
        nascent=f"{REFERENCES}/nascent_all/nascent.txt",
        # No dependencies needed - files are already created by prepare_fastq.smk
        # Script dependency
        script="scripts/run_kb_count.py"
    output:
        kb_info=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}/kb_info.json",
        h5ad=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}/counts_unfiltered/adata.h5ad",
        # BUS files that will be inspected downstream
        bus_output=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}/output.bus",
        bus_unfiltered=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}/output.unfiltered.bus",
        bus_modified=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}/output_modified.unfiltered.bus"
    log:
        f"{LOGS}/kallisto_gex_{{sample_id}}_{{type}}_{{source}}_{{processing}}.log"
    wildcard_constraints:
        sample_id="|".join(GEX_IDS),  # Only match GEX sample IDs
        type="all",  # Only match "all" - prevents conflict with guide rule
        source="main|undetermined|all",
        processing="raw|recovered|merged"
    params:
        out=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}",
        strand=config["kallisto"]["strand"],
        kb_chemistry=config["kallisto"]["kb_chemistry"]
    resources:
        mem_mb=config["resources"]["alignment"]["mem_mb"],
    threads:
        config["resources"]["alignment"]["threads"]
    shell:
        """
        # Clean up any existing tmp directory to prevent conflicts
        rm -rf {params.out}/tmp
        
        # Note: Using --mm for fractional assignment of multimapping reads (even splitting)
        # Alternative: --em uses expectation maximization algorithm (hidden option)
        # Warning: --em can drop reads when genes have no unique reads (s=0 condition in bustools)
        
        # Run with configured strand option
        python scripts/run_kb_count.py \
            nac {params.strand} \
            {input.fq1} {input.fq2} {params.out} {threads} \
            {input.index_idx} {input.t2g} {input.barcodes} {input.replace} \
            {params.kb_chemistry} \
            --cdna {input.cdna} --nascent {input.nascent} --mm &> {log}
        """


rule kallisto_guide:
    input:
        fq1=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R1", wildcards.source, wildcards.processing),
        fq2=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R2", wildcards.source, wildcards.processing),
        barcodes=config["input_paths"]["barcodes_file"],
        replace=config["input_paths"]["replace_file"],
        index_idx=f"kite_references/{os.path.splitext(os.path.basename(config['guide_file']))[0]}_index/mismatch.idx",
        t2g=f"kite_references/{os.path.splitext(os.path.basename(config['guide_file']))[0]}_index/t2g.txt",
        # No dependencies needed - files are already created by prepare_fastq.smk
        # Script dependency
        script="scripts/run_kb_count.py"
    output:
        kb_info=f"{SCRATCH}/{{sample_id}}/kb_guide_{{source}}_{{processing}}/kb_info.json",
        h5ad=f"{SCRATCH}/{{sample_id}}/kb_guide_{{source}}_{{processing}}/counts_unfiltered/adata.h5ad",
        # BUS files that will be inspected downstream
        bus_output=f"{SCRATCH}/{{sample_id}}/kb_guide_{{source}}_{{processing}}/output.bus",
        bus_unfiltered=f"{SCRATCH}/{{sample_id}}/kb_guide_{{source}}_{{processing}}/output.unfiltered.bus",
        bus_modified=f"{SCRATCH}/{{sample_id}}/kb_guide_{{source}}_{{processing}}/output_modified.unfiltered.bus"
    log:
        f"{LOGS}/kallisto_guide_{{sample_id}}_{{source}}_{{processing}}.log"
    wildcard_constraints:
        sample_id="|".join(GUIDE_IDS),  # Only match guide sample IDs
        source="main|undetermined|all",
        processing="raw|recovered|merged"
    params:
        out=f"{SCRATCH}/{{sample_id}}/kb_guide_{{source}}_{{processing}}",
        strand=config["kallisto"]["strand"],
        kb_chemistry=config["kallisto"]["kb_chemistry"]
    resources:
        mem_mb=config["resources"]["alignment"]["mem_mb"],
    threads:
        config["resources"]["alignment"]["threads"]
    shell:
        """
        # Clean up any existing tmp directory to prevent conflicts
        rm -rf {params.out}/tmp
        
        # Run with configured strand option
        python scripts/run_kb_count.py \
            kite {params.strand} \
            {input.fq1} {input.fq2} {params.out} {threads} \
            {input.index_idx} {input.t2g} {input.barcodes} {input.replace} \
            {params.kb_chemistry} &> {log}
        """


rule kallisto_gex_subsampled:
    input:
        fq1=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R1", "main", "raw"),
        fq2=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R2", "main", "raw"),
        read_counts=f"{RESULTS}/{{sample_id}}/counts.txt",
        replace=config["input_paths"]["replace_file"],
        barcodes=config["input_paths"]["barcodes_file"],
        index_idx=f"{REFERENCES}/nascent_all/index.idx",
        t2g=f"{REFERENCES}/nascent_all/t2g.txt",
        cdna=f"{REFERENCES}/nascent_all/cdna.txt",
        nascent=f"{REFERENCES}/nascent_all/nascent.txt",
        script="scripts/run_kb_count.py"
    output:
        kb_info=f"{SCRATCH}/tmp/umi_sat_{{sample_id}}-{{fraction}}/kb_all/kb_info.json",
        h5ad=f"{SCRATCH}/tmp/umi_sat_{{sample_id}}-{{fraction}}/kb_all/counts_unfiltered/adata.h5ad"
    log:
        f"{LOGS}/kallisto_gex_subsampled_{{sample_id}}_{{fraction}}.log"
    wildcard_constraints:
        sample_id="|".join(GEX_IDS),  # Only match GEX sample IDs
        fraction="0\\.1|0\\.25|0\\.5|0\\.75"
    params:
        out=f"{SCRATCH}/tmp/umi_sat_{{sample_id}}-{{fraction}}/kb_all",
        strand=config["kallisto"]["strand"],
        kb_chemistry=config["kallisto"]["kb_chemistry"]
    threads:
        config["resources"]["alignment"]["threads"]
    resources:
        mem_mb=config["resources"]["alignment"]["mem_mb"],
    shell:
        """
        # Clean up any existing tmp directory to prevent conflicts
        rm -rf {params.out}/tmp
        
        # Read total read count and calculate target reads
        R1_COUNT=$(tail -n 1 {input.read_counts} | cut -d',' -f1)
        TARGET_READS=$(echo "$R1_COUNT * {wildcards.fraction}" | bc | cut -d'.' -f1)
        
        # Run with configured strand option
        python scripts/run_kb_count.py \
            nac {params.strand} \
            {input.fq1} {input.fq2} {params.out} {threads} \
            {input.index_idx} {input.t2g} {input.barcodes} {input.replace} \
            {params.kb_chemistry} \
            --cdna {input.cdna} --nascent {input.nascent} --mm \
            --max-reads $TARGET_READS &> {log}
        """


rule kallisto_guide_subsampled:
    input:
        fq1=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R1", "main", "raw"),
        fq2=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R2", "main", "raw"),
        read_counts=f"{RESULTS}/{{sample_id}}/counts.txt",
        replace=config["input_paths"]["replace_file"],
        barcodes=config["input_paths"]["barcodes_file"],
        index_idx=f"kite_references/{os.path.splitext(os.path.basename(config['guide_file']))[0]}_index/mismatch.idx",
        t2g=f"kite_references/{os.path.splitext(os.path.basename(config['guide_file']))[0]}_index/t2g.txt",
        script="scripts/run_kb_count.py"
    output:
        kb_info=f"{SCRATCH}/tmp/umi_sat_{{sample_id}}-{{fraction}}/kb_guide/kb_info.json",
        h5ad=f"{SCRATCH}/tmp/umi_sat_{{sample_id}}-{{fraction}}/kb_guide/counts_unfiltered/adata.h5ad"
    log:
        f"{LOGS}/kallisto_guide_subsampled_{{sample_id}}_{{fraction}}.log"
    wildcard_constraints:
        sample_id="|".join(GUIDE_IDS),  # Only match guide sample IDs
        fraction="0\\.1|0\\.25|0\\.5|0\\.75"
    params:
        out=f"{SCRATCH}/tmp/umi_sat_{{sample_id}}-{{fraction}}/kb_guide",
        strand=config["kallisto"]["strand"],
        kb_chemistry=config["kallisto"]["kb_chemistry"]
    threads:
        config["resources"]["alignment"]["threads"]
    resources:
        mem_mb=config["resources"]["alignment"]["mem_mb"],
    shell:
        """
        # Clean up any existing tmp directory to prevent conflicts
        rm -rf {params.out}/tmp
        
        # Read total read count and calculate target reads
        R1_COUNT=$(tail -n 1 {input.read_counts} | cut -d',' -f1)
        TARGET_READS=$(echo "$R1_COUNT * {wildcards.fraction}" | bc | cut -d'.' -f1)
        
        # Run with configured strand option
        python scripts/run_kb_count.py \
            kite {params.strand} \
            {input.fq1} {input.fq2} {params.out} {threads} \
            {input.index_idx} {input.t2g} {input.barcodes} {input.replace} \
            {params.kb_chemistry} \
            --max-reads $TARGET_READS &> {log}
        """



# =============================================================================
# STAGE 5: POST-ALIGNMENT PROCESSING
# =============================================================================
rule inspect_bus_files:
    input:
        combined_whitelist=f"{REFERENCES}/barcodes_combined.txt",
        kb_info=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}/kb_info.json",
        # Explicit BUS file inputs
        bus_output=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}/output.bus",
        bus_unfiltered=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}/output.unfiltered.bus",
        bus_modified=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}/output_modified.unfiltered.bus"
    output:
        # The key inspect JSON files that calculate_read_statistics needs
        output_inspect=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}/output_inspect.json",
        output_unfiltered_inspect=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}/output.unfiltered_inspect.json",
        output_modified_unfiltered_inspect=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}/output_modified.unfiltered_inspect.json"
    log:
        f"{LOGS}/inspect_bus_{{sample_id}}_{{type}}_{{source}}_{{processing}}.log"
    shell:
        """
        # Inspect each BUS file and create corresponding inspect JSON
        bustools inspect -o {output.output_inspect} -w {input.combined_whitelist} {input.bus_output} 2>&1 | tee -a {log}
        bustools inspect -o {output.output_unfiltered_inspect} -w {input.combined_whitelist} {input.bus_unfiltered} 2>&1 | tee -a {log}
        bustools inspect -o {output.output_modified_unfiltered_inspect} -w {input.combined_whitelist} {input.bus_modified} 2>&1 | tee -a {log}
        """


rule filter_and_annotate_sublibrary:
    input:
        # GEX unfiltered kallisto output - specific files for dependency tracking
        gex_h5ad=f"{SCRATCH}/{{gex_sample_id}}/kb_all_{{source}}_{{processing}}/counts_unfiltered/adata.h5ad",
        # Guide unfiltered kallisto output - specific file for dependency tracking
        guide_h5ad=lambda wildcards: f"{SCRATCH}/{gex_to_guide[wildcards.gex_sample_id]}/kb_guide_{wildcards.source}_{wildcards.processing}/counts_unfiltered/adata.h5ad",
        # Gene annotation table
        gene_table=f"{REFERENCES}/gene_annotations.tsv",
        # Config file
        config_file=workflow.configfiles[0],
        script="scripts/sublibrary_annotation.py"
    output:
        # Filtered output directory containing h5ad and MTX files
        output_dir=directory(f"{SCRATCH}/{{gex_sample_id}}/kb_all_{{source}}_{{processing}}/counts_filtered"),
        # Annotated h5ad file
        annotated_h5ad=f"{SCRATCH}/{{gex_sample_id}}/kb_all_{{source}}_{{processing}}/counts_filtered/adata.h5ad",
        # MTX file as representative output
        mtx_file=f"{SCRATCH}/{{gex_sample_id}}/kb_all_{{source}}_{{processing}}/counts_filtered/cells_x_genes.total.mtx"
    log:
        f"{LOGS}/filter_and_annotate_sublibrary_{{gex_sample_id}}_{{source}}_{{processing}}.log"
    wildcard_constraints:
        gex_sample_id="|".join(GEX_IDS),  # Only match GEX sample IDs
        source="main|undetermined|all",
        processing="raw|recovered|merged"
    params:
        gex_sample=lambda wildcards: extract_sample(wildcards.gex_sample_id),  # Just the sample name
        guide_sample=lambda wildcards: extract_sample(gex_to_guide[wildcards.gex_sample_id]),  # Just the sample name
        gex_sample_id="{gex_sample_id}",  # Full sample_id for lookups
        guide_sample_id=lambda wildcards: gex_to_guide[wildcards.gex_sample_id],  # Full sample_id for lookups
        # Directory paths for the script
        gex_kb_dir=f"{SCRATCH}/{{gex_sample_id}}/kb_all_{{source}}_{{processing}}",
        guide_kb_dir=lambda wildcards: f"{SCRATCH}/{gex_to_guide[wildcards.gex_sample_id]}/kb_guide_{wildcards.source}_{wildcards.processing}"
    resources:
        mem_mb=config["resources"]["analysis"]["mem_mb"],
    shell:
        """
        # Build command
        CMD="python3 {input.script} \
            --gex-kb-dir {params.gex_kb_dir} \
            --guide-kb-dir {params.guide_kb_dir} \
            --output-dir {output.output_dir} \
            --config {input.config_file} \
            --sample-id {params.gex_sample_id} \
            --guide-sample-id {params.guide_sample_id} \
            --source {wildcards.source} \
            --processing {wildcards.processing} \
            --gene-annotation-table {input.gene_table}"
        
        # Execute
        $CMD &> {log}
        """


rule calculate_read_statistics:
    input:
        kb_info=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}/kb_info.json",
        # Get inspect JSON from inspect_bus_files rule output
        inspect_json=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}/output_modified.unfiltered_inspect.json",
        h5ad=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}/counts_unfiltered/adata.h5ad",
        script="scripts/calculate_read_statistics.py"
    output:
        stats=f"{RESULTS}/{{sample_id}}/qc/{{sample_id}}_{{type}}_{{source}}_{{processing}}_read_statistics.tsv"
    log:
        f"{LOGS}/calculate_read_statistics_{{sample_id}}_{{type}}_{{source}}_{{processing}}.log"
    params:
        kb_dir=f"{SCRATCH}/{{sample_id}}/kb_{{type}}_{{source}}_{{processing}}",
        sample_id="{sample_id}",
        config_file=workflow.configfiles[0]
    shell:
        """
        python3 {input.script} \
            {params.kb_dir} \
            --sample-id {params.sample_id} \
            --config {params.config_file} \
            --output {output.stats} &> {log}
        """



# =============================================================================
# STAGE 6: QC AND ANALYSIS
# =============================================================================


rule cell_calling_analysis:
    input:
        kb_info=f"{SCRATCH}/{{sample_id}}/kb_all_{{source}}_{{processing}}/kb_info.json",
        h5ad=f"{SCRATCH}/{{sample_id}}/kb_all_{{source}}_{{processing}}/counts_filtered/adata.h5ad",
        script="scripts/cell_calling_analysis.py"
    output:
        # Directory containing all cell calling outputs
        cell_calling_dir=directory(f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/cell_calling"),
        # Key summary files for downstream use
        summary=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/cell_calling/results.tsv",
        # Default method cell barcodes file (needed by generate_qc_cell_lists)
        default_barcodes=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/cell_calling/{{sample_id}}_{config['cell_calling']['default_method']}_cell_barcodes.txt"
    log:
        f"{LOGS}/cell_calling_{{sample_id}}_{{source}}_{{processing}}.log"
    wildcard_constraints:
        sample_id="|".join(GEX_IDS),  # Only match GEX sample IDs
        source="main|undetermined|all",
        processing="raw|recovered|merged"
    params:
        kb_dir=f"{SCRATCH}/{{sample_id}}/kb_all_{{source}}_{{processing}}",
        sample_id="{sample_id}"
    threads:
        config["resources"]["analysis"]["threads"]
    resources:
        mem_mb=config["resources"]["analysis"]["mem_mb"],
    shell:
        """
        # ml R/4.4  # R is now installed in the kb conda environment via bioconda::bioconductor-dropletutils
        # Rscript is called with --vanilla flag in the Python script to ensure only conda packages are used
        
        python3 {input.script} \
            --h5ad_file {input.h5ad} \
            --kb_dir {params.kb_dir} \
            --sample-id {params.sample_id} \
            --config {workflow.configfiles[0]} \
            --output_dir {output.cell_calling_dir} \
            --ncores {threads} &> {log}
        """


rule cell_calling_plots:
    input:
        h5ad=f"{SCRATCH}/{{sample_id}}/kb_all_{{source}}_{{processing}}/counts_filtered/adata.h5ad",
        cell_calling_dir=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/cell_calling",
        summary=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/cell_calling/results.tsv",
        read_stats=f"{RESULTS}/{{sample_id}}/qc/{{sample_id}}_all_{{source}}_{{processing}}_read_statistics.tsv",
        script="scripts/cell_calling_plots_no_json.py"
    output:
        plot_dir=directory(f"{RESULTS}/preprocessing_report/plots/cell_calling/{{source}}_{{processing}}/{{sample_id}}"),
        complete=f"{RESULTS}/preprocessing_report/plots/cell_calling/{{source}}_{{processing}}/{{sample_id}}.complete"
    params:
        plot_dir=f"{RESULTS}/preprocessing_report/plots",
        sample_id="{sample_id}"
    log:
        f"{LOGS}/cell_calling_plots_{{sample_id}}_{{source}}_{{processing}}.log"
    wildcard_constraints:
        sample_id="|".join(GEX_IDS),  # Only match GEX sample IDs
        source="main|undetermined|all",
        processing="raw|recovered|merged"
    shell:
        """
        python3 {input.script} \
            --h5ad_file {input.h5ad} \
            --cell_calling_dir {input.cell_calling_dir} \
            --read_stats {input.read_stats} \
            --sample-id {params.sample_id} \
            --plot_dir {params.plot_dir} \
            --source {wildcards.source} \
            --processing {wildcards.processing} &> {log}
        
        # Create sentinel file to indicate successful completion
        touch {output.complete}
        """


rule generate_qc_cell_lists:
    """Generate QC cell filtering decisions based on mitochondrial and expression metrics.
    
    Creates a TSV file with boolean filtering decisions for multiple QC methods
    and continuous GMM posterior probabilities for visualization.
    """
    input:
        h5ad=f"{SCRATCH}/{{sample_id}}/kb_all_{{source}}_{{processing}}/counts_filtered/adata.h5ad",
        # Add cell calling results as dependency
        cell_barcodes=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/cell_calling/{{sample_id}}_{config['cell_calling']['default_method']}_cell_barcodes.txt",
        script="scripts/generate_qc_cell_lists.py"
    output:
        qc_lists=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/qc_cell_lists.tsv",
        debug_plot_dir=directory(f"{RESULTS}/preprocessing_report/plots/cell_quality_qc/{{source}}_{{processing}}/{{sample_id}}")
    params:
        config_file=workflow.configfiles[0]
    log:
        f"{LOGS}/generate_qc_cell_lists_{{sample_id}}_{{source}}_{{processing}}.log"
    wildcard_constraints:
        sample_id="|".join(GEX_IDS),  # Only match GEX sample IDs
        source="main|undetermined|all",
        processing="raw|recovered|merged"
    shell:
        """
        python {input.script} \
            --h5ad {input.h5ad} \
            --cell-barcodes {input.cell_barcodes} \
            --output {output.qc_lists} \
            --plot-output {output.debug_plot_dir}/plot.png \
            --config {params.config_file} &> {log}
        """


rule calculate_gmm_thresholds:
    """Calculate GMM thresholds and generate mixture model plots"""
    input:
        h5ad=f"{SCRATCH}/{{sample_id}}/kb_all_{{source}}_{{processing}}/counts_filtered/adata.h5ad",
        cell_calling_dir=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/cell_calling",
        script="scripts/guide_analysis.py"
    output:
        threshold_table=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/qc_metrics/gmm_thresholds_per_cell.tsv"
    log: f"{LOGS}/gmm_thresholds_{{sample_id}}_{{source}}_{{processing}}.log"
    threads: config["resources"]["analysis"]["threads"]
    wildcard_constraints:
        sample_id="|".join(GEX_IDS),
        source="main|undetermined|all",
        processing="raw|recovered|merged"
    params:
        sample_id="{sample_id}",
        config_file=workflow.configfiles[0],
        plot_dir=f"{RESULTS}/preprocessing_report/plots/per_cell/{{source}}_{{processing}}/{{sample_id}}"
    shell:
        """
        python {input.script} \
            --h5ad {input.h5ad} \
            --cell-calling-dir {input.cell_calling_dir} \
            --sample-id {params.sample_id} \
            --config {params.config_file} \
            --output {output.threshold_table} \
            --plot-dir {params.plot_dir} \
            --threads {threads} \
            > {log} 2>&1
        """


rule calculate_qc_metrics_sample:
    """Calculate QC metrics at sample level"""
    input:
        h5ad=f"{SCRATCH}/{{sample_id}}/kb_all_{{source}}_{{processing}}/counts_filtered/adata.h5ad",
        cell_calling_dir=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/cell_calling",
        read_stats=f"{RESULTS}/{{sample_id}}/qc/{{sample_id}}_all_{{source}}_{{processing}}_read_statistics.tsv",
        guide_stats=lambda wildcards: f"{RESULTS}/{gex_to_guide[wildcards.sample_id]}/qc/{gex_to_guide[wildcards.sample_id]}_guide_{wildcards.source}_{wildcards.processing}_read_statistics.tsv",
        qc_lists=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/qc_cell_lists.tsv",
        gmm_thresholds=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/qc_metrics/gmm_thresholds_per_cell.tsv",
        script="scripts/calculate_qc_metrics_by_biological_sample.py"
    output:
        qc_metrics_sample=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/qc_metrics/by_sample.tsv"
    log: f"{LOGS}/qc_metrics_sample_{{sample_id}}_{{source}}_{{processing}}.log"
    threads: config["resources"]["analysis"]["threads"]
    wildcard_constraints:
        sample_id="|".join(GEX_IDS),
        source="main|undetermined|all",
        processing="raw|recovered|merged"
    params:
        sample_id="{sample_id}",
        pool=lambda wildcards: extract_pool(wildcards.sample_id),
        plot_dir=f"{RESULTS}/preprocessing_report/plots/per_cell/{{source}}_{{processing}}/{{sample_id}}"
    shell:
        """
        python3 {input.script} \
            --h5ad {input.h5ad} \
            --cell-calling-dir {input.cell_calling_dir} \
            --read-stats {input.read_stats} \
            --guide-read-stats {input.guide_stats} \
            --gmm-thresholds {input.gmm_thresholds} \
            --sample-id {params.sample_id} \
            --config {workflow.configfiles[0]} \
            --stratify-by sample \
            --output {output.qc_metrics_sample} \
            --plot-dir {params.plot_dir} \
            --pool {params.pool} \
            --source {wildcards.source} \
            --processing {wildcards.processing} \
            --qc-cell-lists {input.qc_lists} \
            --threads {threads} \
            --per-cell-plot-method {config['cell_calling']['default_method']} &> {log}
        """

rule calculate_qc_metrics_biological:
    """Calculate QC metrics at biological sample level"""
    input:
        h5ad=f"{SCRATCH}/{{sample_id}}/kb_all_{{source}}_{{processing}}/counts_filtered/adata.h5ad",
        cell_calling_dir=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/cell_calling",
        read_stats=f"{RESULTS}/{{sample_id}}/qc/{{sample_id}}_all_{{source}}_{{processing}}_read_statistics.tsv",
        guide_stats=lambda wildcards: f"{RESULTS}/{gex_to_guide[wildcards.sample_id]}/qc/{gex_to_guide[wildcards.sample_id]}_guide_{wildcards.source}_{wildcards.processing}_read_statistics.tsv",
        gmm_thresholds=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/qc_metrics/gmm_thresholds_per_cell.tsv",
        script="scripts/calculate_qc_metrics_by_biological_sample.py"
    output:
        qc_metrics_bio=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/qc_metrics/by_biological_sample.tsv"
    log: f"{LOGS}/qc_metrics_biological_{{sample_id}}_{{source}}_{{processing}}.log"
    threads: config["resources"]["analysis"]["threads"]
    wildcard_constraints:
        sample_id="|".join(GEX_IDS),
        source="main|undetermined|all",
        processing="raw|recovered|merged"
    params:
        sample_id="{sample_id}",
        pool=lambda wildcards: extract_pool(wildcards.sample_id),
        plot_dir=f"{RESULTS}/preprocessing_report/plots/per_cell/{{source}}_{{processing}}/{{sample_id}}"
    shell:
        """
        python3 {input.script} \
            --h5ad {input.h5ad} \
            --cell-calling-dir {input.cell_calling_dir} \
            --read-stats {input.read_stats} \
            --guide-read-stats {input.guide_stats} \
            --gmm-thresholds {input.gmm_thresholds} \
            --sample-id {params.sample_id} \
            --config {workflow.configfiles[0]} \
            --stratify-by biological_sample \
            --output {output.qc_metrics_bio} \
            --plot-dir {params.plot_dir} \
            --pool {params.pool} \
            --source {wildcards.source} \
            --processing {wildcards.processing} \
            --threads {threads} \
            --per-cell-plot-method {config['cell_calling']['default_method']} &> {log}
        """

rule calculate_qc_metrics_well:
    """Calculate QC metrics at well level"""
    input:
        h5ad=f"{SCRATCH}/{{sample_id}}/kb_all_{{source}}_{{processing}}/counts_filtered/adata.h5ad",
        cell_calling_dir=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/cell_calling",
        read_stats=f"{RESULTS}/{{sample_id}}/qc/{{sample_id}}_all_{{source}}_{{processing}}_read_statistics.tsv",
        guide_stats=lambda wildcards: f"{RESULTS}/{gex_to_guide[wildcards.sample_id]}/qc/{gex_to_guide[wildcards.sample_id]}_guide_{wildcards.source}_{wildcards.processing}_read_statistics.tsv",
        gmm_thresholds=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/qc_metrics/gmm_thresholds_per_cell.tsv",
        script="scripts/calculate_qc_metrics_by_biological_sample.py"
    output:
        qc_metrics_well=f"{RESULTS}/preprocessing_report/data/per_sample/{{source}}_{{processing}}/{{sample_id}}/qc_metrics/by_well.tsv"
    log: f"{LOGS}/qc_metrics_well_{{sample_id}}_{{source}}_{{processing}}.log"
    threads: config["resources"]["analysis"]["threads"]
    wildcard_constraints:
        sample_id="|".join(GEX_IDS),
        source="main|undetermined|all",
        processing="raw|recovered|merged"
    params:
        sample_id="{sample_id}",
        pool=lambda wildcards: extract_pool(wildcards.sample_id),
        plot_dir=f"{RESULTS}/preprocessing_report/plots/per_cell/{{source}}_{{processing}}/{{sample_id}}"
    shell:
        """
        python3 {input.script} \
            --h5ad {input.h5ad} \
            --cell-calling-dir {input.cell_calling_dir} \
            --read-stats {input.read_stats} \
            --guide-read-stats {input.guide_stats} \
            --gmm-thresholds {input.gmm_thresholds} \
            --sample-id {params.sample_id} \
            --config {workflow.configfiles[0]} \
            --stratify-by well \
            --output {output.qc_metrics_well} \
            --plot-dir {params.plot_dir} \
            --pool {params.pool} \
            --source {wildcards.source} \
            --processing {wildcards.processing} \
            --threads {threads} \
            --per-cell-plot-method {config['cell_calling']['default_method']} &> {log}
        """


rule umi_saturation_analysis:
    input:
        # Main data (100% point) - using the main raw output
        kb_info=f"{SCRATCH}/{{sample_id}}/kb_all_main_raw/kb_info.json",
        main_unfiltered=f"{SCRATCH}/{{sample_id}}/kb_all_main_raw/counts_unfiltered/adata.h5ad",
        cell_calling_results=f"{RESULTS}/preprocessing_report/data/per_sample/main_raw/{{sample_id}}/cell_calling/results.tsv",
        read_stats=f"{RESULTS}/{{sample_id}}/qc/{{sample_id}}_all_main_raw_read_statistics.tsv",
        # Subsampled matrices - Snakemake will automatically generate these (unfiltered)
        subsampled_01=f"{SCRATCH}/tmp/umi_sat_{{sample_id}}-0.1/kb_all/counts_unfiltered/adata.h5ad",
        subsampled_025=f"{SCRATCH}/tmp/umi_sat_{{sample_id}}-0.25/kb_all/counts_unfiltered/adata.h5ad",
        subsampled_05=f"{SCRATCH}/tmp/umi_sat_{{sample_id}}-0.5/kb_all/counts_unfiltered/adata.h5ad",
        subsampled_075=f"{SCRATCH}/tmp/umi_sat_{{sample_id}}-0.75/kb_all/counts_unfiltered/adata.h5ad",
        script="scripts/umi_saturation_analysis.py"
    output:
        saturation=f"{RESULTS}/preprocessing_report/data/per_sample/main_raw/{{sample_id}}/saturation/umi_saturation.tsv",
        plot_dir=directory(f"{RESULTS}/preprocessing_report/plots/saturation/main_raw/gex/{{sample_id}}")
    log:
        f"{LOGS}/umi_saturation_gex_{{sample_id}}.log"
    threads: config["resources"]["analysis"]["threads"]
    wildcard_constraints:
        sample_id="|".join(GEX_IDS)  # Only match GEX sample IDs
    params:
        cell_calling_dir=f"{RESULTS}/preprocessing_report/data/per_sample/main_raw/{{sample_id}}/cell_calling",
        sample_id="{sample_id}",
        scratch_dir=SCRATCH
    shell:
        """
        python3 {input.script} \
            --sample-id {params.sample_id} \
            --config {workflow.configfiles[0]} \
            --read_stats {input.read_stats} \
            --output_tsv {output.saturation} \
            --output_plot {output.plot_dir}/{wildcards.sample_id}_umi_saturation.png \
            --cell_calling_dir {params.cell_calling_dir} \
            --sample_type gex \
            --source main \
            --processing raw \
            --scratch-dir {params.scratch_dir} \
            --threads {threads} &> {log}
        """


rule umi_saturation_analysis_guide:
    input:
        # Main data (100% point) - using the main raw output
        kb_info=f"{SCRATCH}/{{sample_id}}/kb_guide_main_raw/kb_info.json",
        main_unfiltered=f"{SCRATCH}/{{sample_id}}/kb_guide_main_raw/counts_unfiltered/adata.h5ad",
        read_stats=f"{RESULTS}/{{sample_id}}/qc/{{sample_id}}_guide_main_raw_read_statistics.tsv",
        # Depend on the paired GEX cell calling results
        cell_calling_results=lambda wildcards: f"{RESULTS}/preprocessing_report/data/per_sample/main_raw/{guide_to_gex[wildcards.sample_id]}/cell_calling/results.tsv",
        # Subsampled guide matrices - Snakemake will automatically generate these (unfiltered)
        subsampled_01=f"{SCRATCH}/tmp/umi_sat_{{sample_id}}-0.1/kb_guide/counts_unfiltered/adata.h5ad",
        subsampled_025=f"{SCRATCH}/tmp/umi_sat_{{sample_id}}-0.25/kb_guide/counts_unfiltered/adata.h5ad",
        subsampled_05=f"{SCRATCH}/tmp/umi_sat_{{sample_id}}-0.5/kb_guide/counts_unfiltered/adata.h5ad",
        subsampled_075=f"{SCRATCH}/tmp/umi_sat_{{sample_id}}-0.75/kb_guide/counts_unfiltered/adata.h5ad",
        script="scripts/umi_saturation_analysis.py"
    output:
        saturation=f"{RESULTS}/preprocessing_report/data/per_sample/main_raw/{{sample_id}}/saturation/guide_umi_saturation.tsv",
        plot_dir=directory(f"{RESULTS}/preprocessing_report/plots/saturation/main_raw/guide/{{sample_id}}")
    log:
        f"{LOGS}/umi_saturation_guide_{{sample_id}}.log"
    threads: config["resources"]["analysis"]["threads"]
    wildcard_constraints:
        sample_id="|".join(GUIDE_IDS)  # Only match guide sample IDs
    params:
        # Cell calling directory is from the paired GEX sample
        cell_calling_dir=lambda wildcards: f"{RESULTS}/preprocessing_report/data/per_sample/main_raw/{guide_to_gex[wildcards.sample_id]}/cell_calling",
        sample_id="{sample_id}",
        scratch_dir=SCRATCH
    shell:
        """
        python3 {input.script} \
            --sample-id {params.sample_id} \
            --config {workflow.configfiles[0]} \
            --read_stats {input.read_stats} \
            --output_tsv {output.saturation} \
            --output_plot {output.plot_dir}/{wildcards.sample_id}_guide_umi_saturation.png \
            --cell_calling_dir {params.cell_calling_dir} \
            --sample_type guide \
            --source main \
            --processing raw \
            --scratch-dir {params.scratch_dir} \
            --threads {threads} &> {log}
        """



# =============================================================================
# STAGE 7: CONSOLIDATION AND VISUALIZATION
# =============================================================================
rule consolidate_qc_metrics:
    """Consolidate QC metrics across all samples using simple concatenation"""
    input:
        # QC metrics after cell calling for all GEX samples (already includes read stats and guide stats)
        qc_metrics=lambda wildcards: [f"{RESULTS}/preprocessing_report/data/per_sample/{wildcards.source}_{wildcards.processing}/{row['sample_id']}/qc_metrics/by_sample.tsv" 
                                      for _, row in load_sample_info(config).iterrows() 
                                      if row['sample_type'] == 'gex'],
        # Biological sample QC metrics
        bio_metrics=lambda wildcards: [f"{RESULTS}/preprocessing_report/data/per_sample/{wildcards.source}_{wildcards.processing}/{row['sample_id']}/qc_metrics/by_biological_sample.tsv" 
                                       for _, row in load_sample_info(config).iterrows() 
                                       if row['sample_type'] == 'gex'],
        # Well QC metrics
        well_metrics=lambda wildcards: [f"{RESULTS}/preprocessing_report/data/per_sample/{wildcards.source}_{wildcards.processing}/{row['sample_id']}/qc_metrics/by_well.tsv" 
                                        for _, row in load_sample_info(config).iterrows() 
                                        if row['sample_type'] == 'gex'],
        # Pool statistics (only for main/raw since undetermined doesn't apply to other combinations)
        pool_stats=lambda wildcards: expand(f"{RESULTS}/preprocessing_report/data/per_pool/main_raw/{{pool}}/pool_statistics.tsv",
                                          pool=load_sample_info(config)['pool'].unique()) if wildcards.source == 'main' and wildcards.processing == 'raw' else [],
        script="scripts/consolidate_qc_simple.py"
    output:
        combined=f"{RESULTS}/preprocessing_report/data/consolidated/{{source}}_{{processing}}/all_metrics.tsv",
        combined_bio=f"{RESULTS}/preprocessing_report/data/consolidated/{{source}}_{{processing}}/by_biological_sample.tsv",
        combined_well=f"{RESULTS}/preprocessing_report/data/consolidated/{{source}}_{{processing}}/by_well.tsv"
    log:
        f"{LOGS}/consolidate_qc_metrics_{{source}}_{{processing}}.log"
    wildcard_constraints:
        source="main|undetermined|all",
        processing="raw|recovered|merged"
    shell:
        """
        # Sample-level consolidation (simple concatenation)
        python {input.script} \
            --input-files {input.qc_metrics} \
            --output {output.combined} \
            --add-pool &> {log}
        
        # Biological sample-level consolidation (simple concatenation)
        python {input.script} \
            --input-files {input.bio_metrics} \
            --output {output.combined_bio} \
            --add-pool &>> {log}
        
        # Well-level consolidation (simple concatenation)
        python {input.script} \
            --input-files {input.well_metrics} \
            --output {output.combined_well} \
            --add-pool &>> {log}
        """


rule visualize_consolidated_qc:
    """Generate visualizations for consolidated QC metrics"""
    input:
        sample_metrics=f"{RESULTS}/preprocessing_report/data/consolidated/{{source}}_{{processing}}/all_metrics.tsv",
        bio_metrics=f"{RESULTS}/preprocessing_report/data/consolidated/{{source}}_{{processing}}/by_biological_sample.tsv",
        well_metrics=f"{RESULTS}/preprocessing_report/data/consolidated/{{source}}_{{processing}}/by_well.tsv",
        script="scripts/visualize_aggregated_qc_metrics_no_json.py"
    output:
        # Track both consolidated plot directories
        plot_dir_general=directory(f"{RESULTS}/preprocessing_report/plots/consolidated_general/{{source}}_{{processing}}"),
        plot_dir_cell_based=directory(f"{RESULTS}/preprocessing_report/plots/consolidated_cell_based/{{source}}_{{processing}}"),
        # Sentinel file to track completion
        complete=f"{RESULTS}/preprocessing_report/plots/consolidated_{{source}}_{{processing}}.complete"
    params:
        output_dir=f"{RESULTS}/preprocessing_report/plots"
    wildcard_constraints:
        source="main|undetermined|all",
        processing="raw|recovered|merged"
    log:
        f"{LOGS}/visualize_consolidated_qc_{{source}}_{{processing}}.log"
    shell:
        """
        echo "Generating QC metric visualizations for {wildcards.source} {wildcards.processing}..." | tee {log}
        
        # Sample-level metrics
        python {input.script} \
            --input {input.sample_metrics} \
            --output_dir {params.output_dir} \
            --stratify_by sample \
            --source {wildcards.source} \
            --processing {wildcards.processing} \
            --threads {threads} &> {log}
        
        # Biological sample metrics
        python {input.script} \
            --input {input.bio_metrics} \
            --output_dir {params.output_dir} \
            --stratify_by biological_sample \
            --source {wildcards.source} \
            --processing {wildcards.processing} \
            --threads {threads} &>> {log}
        
        # Well-level metrics (heatmaps)
        python {input.script} \
            --input {input.well_metrics} \
            --output_dir {params.output_dir} \
            --stratify_by well \
            --source {wildcards.source} \
            --processing {wildcards.processing} \
            --threads {threads} &>> {log}
        
        # Create sentinel file to indicate successful completion
        touch {output.complete}
        """


rule process_pool_metrics:
    """Consolidate and visualize pool-level metrics (only for main/raw)"""
    input:
        pool_stats=expand(f"{RESULTS}/preprocessing_report/data/per_pool/main_raw/{{pool}}/pool_statistics.tsv",
                         pool=load_sample_info(config)['pool'].unique()),
        consolidate_script="scripts/consolidate_qc_simple.py",
        visualize_script="scripts/visualize_aggregated_qc_metrics_no_json.py"
    output:
        combined_pool=f"{RESULTS}/preprocessing_report/data/consolidated/main_raw/by_pool.tsv"
    params:
        output_dir=f"{RESULTS}/preprocessing_report/plots"
    log:
        f"{LOGS}/process_pool_metrics.log"
    shell:
        """
        # Consolidate pool statistics
        python {input.consolidate_script} \
            --input-files {input.pool_stats} \
            --output {output.combined_pool} &> {log}
        
        # Generate visualizations
        python {input.visualize_script} \
            --input {output.combined_pool} \
            --output_dir {params.output_dir} \
            --stratify_by pool \
            --source main \
            --processing raw \
            --threads {threads} &>> {log}
        """



# =============================================================================
# STAGE 8: FINAL OUTPUTS AND PACKAGING
# =============================================================================
rule generate_preprocessing_report:
    """Package QC outputs for laptop viewing with Streamlit dashboard"""
    input:
        # Get ALL outputs from the pipeline (except the report itself to avoid circular dependency)
        # This includes:
        # - Read statistics
        # - Cell calling results and plots  
        # - QC metrics (sample, biological sample, well, pool)
        # - Saturation analysis data and plots (both GEX and guide)
        # - Consolidated metrics and visualizations
        # - All plot directories and completion markers
        all_outputs=get_preprocessing_outputs(config, report_dir="preprocessing_report"),
        # Additional explicit inputs needed by the packaging script
        sample_info=config['sample_info_file'],
        dashboard_script="scripts/streamlit_qc_dashboard_optimized.py",
        package_script="scripts/package_qc_for_laptop_fast_no_json.py"
    output:
        done=f"{RESULTS}/preprocessing_report/DONE.txt"
    params:
        timestamp=lambda wildcards: datetime.now().strftime("%Y%m%d_%H%M%S")
    log:
        f"{LOGS}/generate_preprocessing_report.log"
    shell:
        """
        # Run packaging script with QC report directory
        python -u {input.package_script} \
            --qc-report-dir {RESULTS}/preprocessing_report \
            --sample-info {input.sample_info} \
            --dashboard-script {input.dashboard_script} \
            --output-archive {RESULTS}/qc_dashboard_preprocessing_{params.timestamp}.tar.gz \
            --per-cell-method-filter {config['cell_calling']['default_method']} \
            --guide-cutoff-filter 1,2 \
            --threads {threads} &> {log}
        
        # Create done file
        touch {output.done}
        """



