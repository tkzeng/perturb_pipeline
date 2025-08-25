# FASTQ PREPARATION PIPELINE
# ==========================
# This pipeline handles all FASTQ preprocessing before downstream analysis:
# 1. Lane concatenation for multi-lane samples
# 2. Undetermined read recovery and per-sample splitting  
# 3. Barcode recovery on main and undetermined reads
# 4. Cross-pool concatenation for multi-pool samples
# 5. Final merging of all FASTQ sources
#
# To run:
#   snakemake -s prepare_fastq.smk --configfile config.yaml

import os
import glob
import pandas as pd
from scripts.snakemake_helpers import (
    get_scratch_path, get_results_path, get_logs_path, load_sample_info,
    extract_pool, extract_sample, find_processed_fastq, find_raw_fastq_files, get_raw_data_path, get_sample_ids,
    get_undetermined_lane_files, get_main_lane_files
)

# Wrapper function for find_processed_fastq
def find_fastq_file_wrapper(sample_id, read, source, processing):
    return find_processed_fastq(sample_id, read, source, processing, config, SCRATCH)

def get_lane_file_paths(pool, run_id, lane_base, read, config):
    """Get file path for specific undetermined lane file"""
    lane_files = get_undetermined_lane_files(pool, config)
    for r_id, l_base, r1_path, r2_path in lane_files:
        if r_id == run_id and l_base == lane_base:
            return r1_path if read == "R1" else r2_path
    raise ValueError(f"Could not find {read} file for {pool}/{run_id}/{lane_base}")

def get_main_lane_file_paths(sample_id, run_id, lane_base, read, config):
    """Get file path for specific main lane file"""
    lane_files = get_main_lane_files(sample_id, config)
    for r_id, l_base, r1_path, r2_path in lane_files:
        if r_id == run_id and l_base == lane_base:
            return r1_path if read == "R1" else r2_path
    raise ValueError(f"Could not find {read} file for {sample_id}/{run_id}/{lane_base}")

# Define base paths
SCRATCH = get_scratch_path(config=config)
RESULTS = get_results_path(config=config)
LOGS = get_logs_path(config=config)

# Create logs directory
os.makedirs(LOGS, exist_ok=True)

# Get unique pools from sample sheet
sample_df = load_sample_info(config)
POOLS = sample_df['pool'].unique().tolist()


rule all:
    input:
        # Raw concatenated files (main/raw) - from lane concatenation
        expand(f"{SCRATCH}/concatenated_lanes/{{sample_id}}_R1.fastq.gz", sample_id=get_sample_ids(config)),
        expand(f"{SCRATCH}/concatenated_lanes/{{sample_id}}_R2.fastq.gz", sample_id=get_sample_ids(config)),
        # Undetermined fastqs per sample (undetermined/raw)
        expand(f"{SCRATCH}/undetermined_fastqs/{{sample_id}}_R1.fastq.gz", sample_id=get_sample_ids(config)),
        expand(f"{SCRATCH}/undetermined_fastqs/{{sample_id}}_R2.fastq.gz", sample_id=get_sample_ids(config)),
        # Barcode recovery outputs (main/recovered and undetermined/recovered)
        expand(f"{SCRATCH}/barcode_recovery/main/{{sample_id}}_merged_recovered_R1.fastq.gz", sample_id=get_sample_ids(config)),
        expand(f"{SCRATCH}/barcode_recovery/main/{{sample_id}}_merged_recovered_R2.fastq.gz", sample_id=get_sample_ids(config)),
        expand(f"{SCRATCH}/barcode_recovery/undetermined/{{sample_id}}_recovered_R1.fastq.gz", sample_id=get_sample_ids(config)),
        expand(f"{SCRATCH}/barcode_recovery/undetermined/{{sample_id}}_recovered_R2.fastq.gz", sample_id=get_sample_ids(config)),
        # Final merged files (all/merged)
        expand(f"{SCRATCH}/merged_fastqs/{{sample_id}}_all_merged_R1.fastq.gz", sample_id=get_sample_ids(config)),
        expand(f"{SCRATCH}/merged_fastqs/{{sample_id}}_all_merged_R2.fastq.gz", sample_id=get_sample_ids(config)),
        # Read count files
        expand(f"{RESULTS}/{{sample_id}}/counts.txt", sample_id=get_sample_ids(config)),
        # Pool statistics
        expand(f"{RESULTS}/preprocessing_report/data/per_pool/main_raw/{{pool}}/pool_statistics.tsv", pool=load_sample_info(config)['pool'].unique()),
        # Lane-level FASTQ QC reports
        #[f"{RESULTS}/lane_qc/{sample_info['sample_id']}/{run_id}/{lane_base}/fastp_report.json"
        # for _, sample_info in load_sample_info(config).iterrows()
        # for run_id, lane_base, r1_path, r2_path in get_main_lane_files(sample_info['sample_id'], config)],
        #[f"{RESULTS}/lane_qc/{sample_info['sample_id']}/{run_id}/{lane_base}/fastp_report.html"
        # for _, sample_info in load_sample_info(config).iterrows()
        # for run_id, lane_base, r1_path, r2_path in get_main_lane_files(sample_info['sample_id'], config)]

rule concatenate_lanes:
    """Concatenate FASTQ files from multiple lanes for the same sample (or copy single-lane files)"""
    input:
        r1_files=lambda wildcards: find_raw_fastq_files(wildcards.sample_id, "R1", config),
        r2_files=lambda wildcards: find_raw_fastq_files(wildcards.sample_id, "R2", config)
    output:
        r1=f"{SCRATCH}/concatenated_lanes/{{sample_id}}_R1.fastq.gz",
        r2=f"{SCRATCH}/concatenated_lanes/{{sample_id}}_R2.fastq.gz"
    log:
        f"{LOGS}/concatenate_lanes_{{sample_id}}.log"
    threads: config["resources"]["single"]["threads"]
    resources:
        mem_mb=config["resources"]["single"]["mem_mb"]
    shell:
        """
        # Concatenate R1 files (already sorted by lane in input)
        cat {input.r1_files} > {output.r1} 2> {log}
        
        # Concatenate R2 files (already sorted by lane in input)
        cat {input.r2_files} > {output.r2} 2>> {log}
        
        echo "Concatenated $(echo '{input.r1_files}' | wc -w) lane files for sample {wildcards.sample_id}" >> {log}
        """



rule recover_undetermined_per_lane:
    """Process individual undetermined lane files (one job per lane)"""
    input:
        script="scripts/recover_undetermined_barcodes_simple.py"
    output:
        report=f"{SCRATCH}/undetermined_index/{{pool}}/{{run_id}}/{{lane_base}}/recovery_summary.txt"
    log:
        f"{LOGS}/undetermined_recovery_{{pool}}_{{run_id}}_{{lane_base}}.log"
    threads: config["resources"]["light"]["threads"]
    resources:
        mem_mb=config["resources"]["light"]["mem_mb"]
    params:
        r1_file=lambda wildcards: get_lane_file_paths(wildcards.pool, wildcards.run_id, wildcards.lane_base, "R1", config),
        r2_file=lambda wildcards: get_lane_file_paths(wildcards.pool, wildcards.run_id, wildcards.lane_base, "R2", config),
        outdir=f"{SCRATCH}/undetermined_index/{{pool}}/{{run_id}}/{{lane_base}}",
        max_open_files=config["undetermined_recovery"]["max_open_files"]
    shell:
        """
        mkdir -p {params.outdir}
        
        python3 -u {input.script} \
            --fastq-r1 "{params.r1_file}" \
            --fastq-r2 "{params.r2_file}" \
            --indices "{config[input_paths][primer_info_file]}" \
            --output-dir {params.outdir} &> {log}
        """


rule create_undetermined_fastq:
    input:
        lane_reports=lambda wildcards: [f"{SCRATCH}/undetermined_index/{extract_pool(wildcards.sample_id)}/{run_id}/{lane_base}/recovery_summary.txt"
                                       for run_id, lane_base, r1_path, r2_path in get_undetermined_lane_files(extract_pool(wildcards.sample_id), config)],
        sample_info=config["sample_info_file"],
        primer_info=config["input_paths"]["primer_info_file"],
        script="scripts/create_undetermined_single_sample.py"
    output:
        r1=f"{SCRATCH}/undetermined_fastqs/{{sample_id}}_R1.fastq.gz",
        r2=f"{SCRATCH}/undetermined_fastqs/{{sample_id}}_R2.fastq.gz"
    log:
        f"{LOGS}/create_undetermined_{{sample_id}}.log"
    threads: config["resources"]["single"]["threads"]
    resources:
        mem_mb=config["resources"]["single"]["mem_mb"]
    params:
        pool=lambda wildcards: extract_pool(wildcards.sample_id),
        lane_dirs=lambda wildcards: [f"{SCRATCH}/undetermined_index/{extract_pool(wildcards.sample_id)}/{run_id}/{lane_base}"
                                    for run_id, lane_base, r1_path, r2_path in get_undetermined_lane_files(extract_pool(wildcards.sample_id), config)],
        output_dir=f"{SCRATCH}/undetermined_fastqs",
        sample_id="{sample_id}"
    shell:
        """
        # Run the script with multiple recovery directories (no temp directory needed)
        python {input.script} \
            --sample-info {input.sample_info} \
            --primer-info {input.primer_info} \
            --recovery-dirs {params.lane_dirs} \
            --sample-id {params.sample_id} \
            --output-r1 {output.r1} \
            --output-r2 {output.r2} \
            --min-reads 1000 &> {log}
        """

rule recover_barcodes_main_per_lane:
    """Process individual main lane files (one job per lane) - for main reads only"""
    input:
        barcodes=config["input_paths"]["barcodes_file"],
        script="scripts/recover_barcodes_simple.py"
    output:
        r1_recovered=f"{SCRATCH}/barcode_recovery/main/{{sample_id}}/{{run_id}}/{{lane_base}}_recovered_R1.fastq.gz",
        r2_recovered=f"{SCRATCH}/barcode_recovery/main/{{sample_id}}/{{run_id}}/{{lane_base}}_recovered_R2.fastq.gz",
        stats=f"{SCRATCH}/barcode_recovery/main/{{sample_id}}/{{run_id}}/{{lane_base}}_stats.txt"
    log:
        f"{LOGS}/barcode_recovery_main_{{sample_id}}_{{run_id}}_{{lane_base}}.log"
    threads: config["resources"]["light"]["threads"]
    resources:
        mem_mb=config["resources"]["light"]["mem_mb"]
    params:
        r1_file=lambda wildcards: get_main_lane_file_paths(wildcards.sample_id, wildcards.run_id, wildcards.lane_base, "R1", config),
        r2_file=lambda wildcards: get_main_lane_file_paths(wildcards.sample_id, wildcards.run_id, wildcards.lane_base, "R2", config),
        outdir=f"{SCRATCH}/barcode_recovery/main/{{sample_id}}/{{run_id}}",
        output_prefix=f"{SCRATCH}/barcode_recovery/main/{{sample_id}}/{{run_id}}/{{lane_base}}",
        max_reads=config["barcode_recovery"]["max_reads_recovery"],
        max_shift=config["barcode_recovery"]["max_shift"],
        chemistry=config["kallisto"]["kb_chemistry"]
    shell:
        """
        mkdir -p {params.outdir}
        
        # Build command with optional parameters
        CMD="python {input.script} \
            {params.r1_file} \
            {params.r2_file} \
            {params.output_prefix} \
            --barcode-file {input.barcodes} \
            --chemistry {params.chemistry} \
            --max-shift {params.max_shift}"
        
        # Add max-reads if specified
        if [ -n "{params.max_reads}" ] && [ "{params.max_reads}" != "None" ]; then
            CMD="$CMD --max-reads {params.max_reads}"
        fi
        
        $CMD &> {log}
        """

rule merge_barcode_recovery_main:
    """Merge main barcode recovery results from all lanes for a sample"""
    input:
        lane_stats=lambda wildcards: [f"{SCRATCH}/barcode_recovery/main/{wildcards.sample_id}/{run_id}/{lane_base}_stats.txt"
                                     for run_id, lane_base, r1_path, r2_path in get_main_lane_files(wildcards.sample_id, config)],
        lane_r1=lambda wildcards: [f"{SCRATCH}/barcode_recovery/main/{wildcards.sample_id}/{run_id}/{lane_base}_recovered_R1.fastq.gz"
                                  for run_id, lane_base, r1_path, r2_path in get_main_lane_files(wildcards.sample_id, config)],
        lane_r2=lambda wildcards: [f"{SCRATCH}/barcode_recovery/main/{wildcards.sample_id}/{run_id}/{lane_base}_recovered_R2.fastq.gz"
                                  for run_id, lane_base, r1_path, r2_path in get_main_lane_files(wildcards.sample_id, config)]
    output:
        r1_recovered=f"{SCRATCH}/barcode_recovery/main/{{sample_id}}_merged_recovered_R1.fastq.gz",
        r2_recovered=f"{SCRATCH}/barcode_recovery/main/{{sample_id}}_merged_recovered_R2.fastq.gz",
        stats=f"{SCRATCH}/barcode_recovery/main/{{sample_id}}_merged_stats.txt"
    log:
        f"{LOGS}/merge_barcode_recovery_main_{{sample_id}}.log"
    threads: config["resources"]["single"]["threads"]
    resources:
        mem_mb=config["resources"]["single"]["mem_mb"]
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.r1_recovered})
        
        # Concatenate all lane R1 files
        cat {input.lane_r1} > {output.r1_recovered} 2> {log}
        
        # Concatenate all lane R2 files  
        cat {input.lane_r2} > {output.r2_recovered} 2>> {log}
        
        # Merge stats files
        echo "Merged main barcode recovery results from $(echo {input.lane_stats} | wc -w) lanes" > {output.stats}
        echo "Lanes processed:" >> {output.stats}
        
        for stats_file in {input.lane_stats}; do
            echo "  $stats_file" >> {output.stats}
            echo "--- Content from $stats_file ---" >> {output.stats}
            cat "$stats_file" >> {output.stats}
            echo "" >> {output.stats}
        done
        
        echo "Merge completed at $(date)" >> {output.stats}
        """

rule recover_barcodes_undetermined:
    """Process undetermined reads per sample (simple per-sample processing)"""
    input:
        fq1=f"{SCRATCH}/undetermined_fastqs/{{sample_id}}_R1.fastq.gz",
        fq2=f"{SCRATCH}/undetermined_fastqs/{{sample_id}}_R2.fastq.gz",
        barcodes=config["input_paths"]["barcodes_file"],
        script="scripts/recover_barcodes_simple.py"
    output:
        r1_recovered=f"{SCRATCH}/barcode_recovery/undetermined/{{sample_id}}_recovered_R1.fastq.gz",
        r2_recovered=f"{SCRATCH}/barcode_recovery/undetermined/{{sample_id}}_recovered_R2.fastq.gz",
        stats=f"{SCRATCH}/barcode_recovery/undetermined/{{sample_id}}_stats.txt"
    log:
        f"{LOGS}/barcode_recovery_undetermined_{{sample_id}}.log"
    threads: config["resources"]["light"]["threads"]
    resources:
        mem_mb=config["resources"]["light"]["mem_mb"]
    params:
        outdir=f"{SCRATCH}/barcode_recovery/undetermined",
        output_prefix=f"{SCRATCH}/barcode_recovery/undetermined/{{sample_id}}",
        max_reads=config["barcode_recovery"]["max_reads_recovery"],
        max_shift=config["barcode_recovery"]["max_shift"],
        chemistry=config["kallisto"]["kb_chemistry"]
    shell:
        """
        mkdir -p {params.outdir}
        
        # Build command with optional parameters
        CMD="python {input.script} \
            {input.fq1} \
            {input.fq2} \
            {params.output_prefix} \
            --barcode-file {input.barcodes} \
            --chemistry {params.chemistry} \
            --max-shift {params.max_shift}"
        
        # Add max-reads if specified
        if [ -n "{params.max_reads}" ] && [ "{params.max_reads}" != "None" ]; then
            CMD="$CMD --max-reads {params.max_reads}"
        fi
        
        $CMD &> {log}
        """


rule count_reads:
    input:
        fq1=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R1", "main", "raw"),
        fq2=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R2", "main", "raw")
    output:
        f"{RESULTS}/{{sample_id}}/counts.txt"
    log:
        f"{LOGS}/count_reads_{{sample_id}}.log"
    threads: config["resources"]["light"]["threads"]
    resources:
        mem_mb=config["resources"]["light"]["mem_mb"]
    shell:
        """
        echo "R1_reads,R2_reads" > {output}
        R1_COUNT=$(pigz -cd -p {threads} {input.fq1} | echo $((`wc -l`/4)))
        R2_COUNT=$(pigz -cd -p {threads} {input.fq2} | echo $((`wc -l`/4)))
        echo "$R1_COUNT,$R2_COUNT" >> {output}
        """


rule calculate_pool_statistics:
    """Calculate pool-level statistics including undetermined read fraction"""
    input:
        # All sample count files for the pool
        sample_counts=lambda wildcards: [f"{RESULTS}/{row['sample_id']}/counts.txt"
                                        for _, row in load_sample_info(config).iterrows() 
                                        if row['pool'] == wildcards.pool],
        script="scripts/calculate_pool_statistics.py"
    output:
        stats=f"{RESULTS}/preprocessing_report/data/per_pool/main_raw/{{pool}}/pool_statistics.tsv"
    log:
        f"{LOGS}/calculate_pool_statistics_{{pool}}.log"
    threads: config["resources"]["light"]["threads"]
    resources:
        mem_mb=config["resources"]["light"]["mem_mb"]
    params:
        pool="{pool}",
        config_file=workflow.configfiles[0]
    shell:
        """
        python3 {input.script} \
            --pool {params.pool} \
            --sample-counts {input.sample_counts} \
            --config {params.config_file} \
            --output {output.stats} &> {log}
        """


rule fastp_qc_lane:
    input:
        r1_file=lambda wildcards: get_main_lane_file_paths(wildcards.sample_id, wildcards.run_id, wildcards.lane_base, "R1", config),
        r2_file=lambda wildcards: get_main_lane_file_paths(wildcards.sample_id, wildcards.run_id, wildcards.lane_base, "R2", config)
    output:
        json=f"{RESULTS}/lane_qc/{{sample_id}}/{{run_id}}/{{lane_base}}/fastp_report.json",
        html=f"{RESULTS}/lane_qc/{{sample_id}}/{{run_id}}/{{lane_base}}/fastp_report.html"
    params:
        outdir=f"{RESULTS}/lane_qc/{{sample_id}}/{{run_id}}/{{lane_base}}",
        pool=lambda wildcards: extract_pool(wildcards.sample_id),
        sample=lambda wildcards: extract_sample(wildcards.sample_id)
    log:
        f"{LOGS}/fastp_qc_lane_{{sample_id}}_{{run_id}}_{{lane_base}}.log"
    threads: config["resources"]["single"]["threads"]
    resources:
        mem_mb=config["resources"]["single"]["mem_mb"]
    shell:
        """
        mkdir -p {params.outdir}
        
        fastp -i {input.r1_file} -I {input.r2_file} \
              --json {output.json} --html {output.html} \
              --thread {threads} \
              --dont_eval_duplication \
              --report_title "{params.pool}_{params.sample}_{wildcards.run_id}_{wildcards.lane_base}_Lane_QC"
        """


rule merge_all_fastqs:
    """Merge all FASTQs for a sample: main raw + main recovered + undetermined raw + undetermined recovered"""
    input:
        main_raw_r1=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R1", "main", "raw"),
        main_raw_r2=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R2", "main", "raw"),
        main_recovered_r1=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R1", "main", "recovered"),
        main_recovered_r2=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R2", "main", "recovered"),
        undetermined_raw_r1=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R1", "undetermined", "raw"),
        undetermined_raw_r2=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R2", "undetermined", "raw"),
        undetermined_recovered_r1=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R1", "undetermined", "recovered"),
        undetermined_recovered_r2=lambda wildcards: find_fastq_file_wrapper(wildcards.sample_id, "R2", "undetermined", "recovered")
    output:
        merged_r1=f"{SCRATCH}/merged_fastqs/{{sample_id}}_all_merged_R1.fastq.gz",
        merged_r2=f"{SCRATCH}/merged_fastqs/{{sample_id}}_all_merged_R2.fastq.gz"
    log:
        f"{LOGS}/merge_all_fastqs_{{sample_id}}.log"
    threads: config["resources"]["single"]["threads"]
    resources:
        mem_mb=config["resources"]["single"]["mem_mb"]
    shell:
        """
        # Merge all R1 files (main raw + main recovered + undetermined raw + undetermined recovered)
        cat {input.main_raw_r1} {input.main_recovered_r1} {input.undetermined_raw_r1} {input.undetermined_recovered_r1} > {output.merged_r1} 2> {log}
        
        # Merge all R2 files  
        cat {input.main_raw_r2} {input.main_recovered_r2} {input.undetermined_raw_r2} {input.undetermined_recovered_r2} > {output.merged_r2} 2>> {log}
        
        echo "Successfully merged files for {wildcards.sample_id}" >> {log}
        """

