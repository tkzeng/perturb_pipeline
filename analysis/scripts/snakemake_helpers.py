"""
Helper functions for Snakemake pipeline
Contains path builders, sample management, and utility functions
"""

import os
import sys
import glob
import pandas as pd
import subprocess


# =============================================================================
# Environment-aware path helper functions
# =============================================================================

def get_scratch_path(*parts, config=None):
    """Build scratch directory paths using $SCRATCH environment variable and analysis_name"""
    SCRATCH_BASE = os.environ['SCRATCH']
    if not config:
        raise ValueError("config is required for get_scratch_path")
    analysis_name = config['analysis_name']
    
    scratch_dir = os.path.join(SCRATCH_BASE, analysis_name)
    if parts and parts[0]:  # Only join if parts are provided and not empty
        path_parts = [str(p) for p in parts]
        return os.path.join(scratch_dir, *path_parts)
    return scratch_dir



def get_results_path(*parts, config=None):
    """Build results directory paths automatically from analysis_name"""
    if not config:
        raise ValueError("config is required for get_results_path")
    
    # Auto-generate from current directory and analysis_name
    analysis_name = config['analysis_name']
    results_base = os.path.join(os.getcwd(), f"results_{analysis_name}")
    
    if parts and parts[0]:
        path_parts = [str(p) for p in parts]
        return os.path.join(results_base, *path_parts)
    return results_base


def get_logs_path(*parts, config=None):
    """Build logs directory paths automatically from analysis_name"""
    if not config:
        raise ValueError("config is required for get_logs_path")
    
    # Auto-generate from current directory and analysis_name
    analysis_name = config['analysis_name']
    logs_base = os.path.join(os.getcwd(), f"logs_{analysis_name}")
    
    if parts and parts[0]:
        path_parts = [str(p) for p in parts]
        return os.path.join(logs_base, *path_parts)
    return logs_base


def get_reference_path(*parts, config=None):
    """Build reference directory paths"""
    if not config:
        raise ValueError("config is required for get_reference_path")
    ref_base = config['input_paths']['reference_base']
    
    if parts and parts[0]:
        path_parts = [str(p) for p in parts]
        return os.path.join(ref_base, *path_parts)
    return ref_base


def get_raw_data_path(pool, config):
    """Get raw data directory paths for a pool from sample_info TSV
    
    Returns list of fastq_dir paths (semicolon-separated values split into list)
    """
    sample_df = load_sample_info(config)
    pool_samples = sample_df[sample_df['pool'] == pool]
    if pool_samples.empty:
        raise ValueError(f"Could not find fastq_dir for pool {pool} in sample_info.tsv")
    
    fastq_dir = pool_samples.iloc[0]['fastq_dir']
    # Split semicolon-separated directories
    return [dir_path.strip() for dir_path in fastq_dir.split(';')]


def get_undetermined_lane_files(pool, config):
    """Get all undetermined lane files for a pool across all runs
    
    STRICT VALIDATION: ALL directories must contain undetermined files
    
    Returns list of tuples: (run_id, lane_file_base, full_r1_path, full_r2_path)
    Example: ("run0", "Undetermined_S0_L001", "/path/run1/Undetermined_S0_L001_R1_001.fastq.gz", "/path/run1/Undetermined_S0_L001_R2_001.fastq.gz")
    """
    fastq_dirs = get_raw_data_path(pool, config)
    lane_files = []
    
    for run_idx, fastq_dir in enumerate(fastq_dirs):
        # Find all undetermined R1 files in this directory
        r1_files = sorted(glob.glob(os.path.join(fastq_dir, "*Undetermined*R1*.fastq.gz")))
        
        if not r1_files:
            raise FileNotFoundError(f"No undetermined R1 files found in directory: {fastq_dir}\nExpected pattern: *Undetermined*R1*.fastq.gz")
        
        run_lane_count = 0
        for r1_file in r1_files:
            # Extract base name (e.g., "Undetermined_S0_L001")
            r1_basename = os.path.basename(r1_file)
            # Remove _R1_001.fastq.gz suffix to get base
            lane_base = r1_basename.replace("_R1_001.fastq.gz", "")
            
            # Construct corresponding R2 file path
            r2_file = r1_file.replace("_R1_", "_R2_")
            
            # Verify R2 file exists
            if os.path.exists(r2_file):
                run_id = f"run{run_idx}"
                lane_files.append((run_id, lane_base, r1_file, r2_file))
                run_lane_count += 1
            else:
                raise FileNotFoundError(f"Missing R2 file for undetermined R1: {r1_file}\nExpected R2 file: {r2_file}")
        
        if run_lane_count == 0:
            raise FileNotFoundError(f"No valid undetermined R1/R2 file pairs found in directory: {fastq_dir}")
    
    if not lane_files:
        raise FileNotFoundError(f"No undetermined lane files found for pool '{pool}' in any directory")
    
    return lane_files


def get_main_lane_files(sample_id, config):
    """Get all main (non-undetermined) lane files for a sample across all runs
    
    STRICT VALIDATION: Sample must exist in ALL specified directories
    
    Returns list of tuples: (run_id, lane_file_base, full_r1_path, full_r2_path)
    Example: ("run0", "GEX_S1_L001", "/path/run1/GEX_S1_L001_R1_001.fastq.gz", "/path/run1/GEX_S1_L001_R2_001.fastq.gz")
    """
    sample_df = load_sample_info(config)
    sample_row = sample_df[sample_df['sample_id'] == sample_id]
    if sample_row.empty:
        raise ValueError(f"Sample ID '{sample_id}' not found in sample info file")
    
    fastq_dir_str = sample_row.iloc[0]['fastq_dir']
    sample_name = extract_sample(sample_id)
    
    # Split semicolon-separated directories
    fastq_dirs = [dir_path.strip() for dir_path in fastq_dir_str.split(';')]
    
    lane_files = []
    for run_idx, fastq_dir in enumerate(fastq_dirs):
        # Find all R1 files for this sample in this directory
        r1_pattern = os.path.join(fastq_dir, f"{sample_name}*_S*_R1_001.fastq.gz")
        r1_files = sorted(glob.glob(r1_pattern))
        
        if not r1_files:
            raise FileNotFoundError(f"Sample '{sample_name}' R1 files not found in directory: {fastq_dir}\nExpected pattern: {sample_name}*_S*_R1_001.fastq.gz")
        
        run_lane_count = 0
        for r1_file in r1_files:
            # Extract base name (e.g., "GEX_S1_L001") 
            r1_basename = os.path.basename(r1_file)
            # Remove _R1_001.fastq.gz suffix to get base
            lane_base = r1_basename.replace("_R1_001.fastq.gz", "")
            
            # Construct corresponding R2 file path
            r2_file = r1_file.replace("_R1_", "_R2_")
            
            # Verify R2 file exists
            if os.path.exists(r2_file):
                run_id = f"run{run_idx}"
                lane_files.append((run_id, lane_base, r1_file, r2_file))
                run_lane_count += 1
            else:
                raise FileNotFoundError(f"Missing R2 file for sample R1: {r1_file}\nExpected R2 file: {r2_file}")
        
        if run_lane_count == 0:
            raise FileNotFoundError(f"No valid R1/R2 file pairs found for sample '{sample_name}' in directory: {fastq_dir}")
    
    if not lane_files:
        raise FileNotFoundError(f"No lane files found for sample '{sample_id}' in any directory")
    
    return lane_files


def print_path_configuration(config):
    """Print path configuration for debugging"""
    print("=" * 60, file=sys.stderr)
    print("PATH CONFIGURATION", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"Scratch directory: {get_scratch_path(config=config)}", file=sys.stderr)
    print(f"Results directory: {get_results_path(config=config)}", file=sys.stderr)
    print(f"Logs directory: {get_logs_path(config=config)}", file=sys.stderr)
    print(f"Reference base: {get_reference_path(config=config)}", file=sys.stderr)
    print("=" * 60, file=sys.stderr)


# =============================================================================
# Sample management functions
# =============================================================================

# Cache for sample info to avoid multiple file reads
_sample_info_cache = None


def load_sample_info(config):
    """Load sample information from TSV file"""
    global _sample_info_cache
    if _sample_info_cache is not None:
        return _sample_info_cache
    
    sample_info_file = config["sample_info_file"]
    if not os.path.exists(sample_info_file):
        raise FileNotFoundError(f"Sample info file not found: {sample_info_file}")
    
    # Read TSV file
    df = pd.read_csv(sample_info_file, sep='\t')
    
    # Filter out 'other' samples - we don't process these
    df = df[df['sample_type'].isin(['gex', 'guide'])]
    
    # Filter by sample_ids config if specified
    if config["sample_ids"]:
        df = df[df['sample_id'].isin(config["sample_ids"])]
    
    _sample_info_cache = df
    return df


def get_sample_ids(config):
    """Get sample IDs from sample info or config"""
    if config["sample_ids"]:
        return config["sample_ids"]
    
    sample_df = load_sample_info(config)
    if sample_df.empty:
        raise ValueError("No samples found in sample info file")
    return sample_df['sample_id'].tolist()


def get_samples_by_type(sample_type, config):
    """Get sample IDs filtered by type"""
    sample_df = load_sample_info(config)
    return sample_df[sample_df['sample_type'] == sample_type]['sample_id'].tolist()


def extract_pool(sample_id):
    """Extract pool from sample_id format 'pool:sample'"""
    if ':' not in sample_id:
        raise ValueError(f"Invalid sample_id format: {sample_id}. Expected 'pool:sample'")
    return sample_id.split(':')[0]


def extract_sample(sample_id):
    """Extract sample name from sample_id format 'pool:sample'"""
    if ':' not in sample_id:
        raise ValueError(f"Invalid sample_id format: {sample_id}. Expected 'pool:sample'")
    return sample_id.split(':')[1]


def get_sample_pool(sample_id):
    """Extract pool from sample_id (e.g., 'pool1:gex_1' -> 'pool1')"""
    return extract_pool(sample_id)


# =============================================================================
# File finding functions
# =============================================================================

def find_raw_fastq_files(sample_id, read, config):
    """Find raw FASTQ files in sequencing directories (used by prepare_fastq.smk)
    
    STRICT VALIDATION: Sample must exist in ALL specified directories
    """
    sample_df = load_sample_info(config)
    sample_row = sample_df[sample_df['sample_id'] == sample_id]
    if sample_row.empty:
        raise ValueError(f"Sample ID '{sample_id}' not found in sample info file")
    
    fastq_dir_str = sample_row.iloc[0]['fastq_dir']
    sample_name = extract_sample(sample_id)
    
    # Split semicolon-separated directories
    fastq_dirs = [dir_path.strip() for dir_path in fastq_dir_str.split(';')]
    
    all_matches = []
    for fastq_dir in fastq_dirs:
        # Try multi-lane pattern first
        pattern = os.path.join(fastq_dir, f"{sample_name}*_S*_L*_{read}_001.fastq.gz")
        matches = sorted(glob.glob(pattern))
        if matches:
            all_matches.extend(matches)
        else:
            # Fallback to single-lane pattern
            pattern = os.path.join(fastq_dir, f"{sample_name}*_S*_{read}_001.fastq.gz")
            matches = sorted(glob.glob(pattern))
            if matches:
                all_matches.extend(matches)
            else:
                raise FileNotFoundError(f"Sample '{sample_name}' {read} files not found in directory: {fastq_dir}\nTried patterns: *_S*_L*_{read}_001.fastq.gz and *_S*_{read}_001.fastq.gz")
    
    if not all_matches:
        raise FileNotFoundError(f"No {read} files found for sample '{sample_id}' in any directory")
    
    return all_matches


def find_processed_fastq(sample_id, read, source, processing, config, scratch_base):
    """Find processed FASTQ file in scratch directories (used by preprocessing.smk)"""
    # Define path patterns
    path_patterns = {
        ("main", "raw"): None,  # Special handling below
        ("main", "recovered"): os.path.join(scratch_base, "barcode_recovery", "main", f"{sample_id}_merged_recovered_{read}.fastq.gz"),
        ("undetermined", "raw"): os.path.join(scratch_base, "undetermined_fastqs", f"{sample_id}_{read}.fastq.gz"),
        ("undetermined", "recovered"): os.path.join(scratch_base, "barcode_recovery", "undetermined", f"{sample_id}_recovered_{read}.fastq.gz"),
        ("all", "merged"): os.path.join(scratch_base, "merged_fastqs", f"{sample_id}_all_merged_{read}.fastq.gz")
    }
    
    # Special handling for main/raw - always return concatenated path
    if source == "main" and processing == "raw":
        return os.path.join(scratch_base, "concatenated_lanes", f"{sample_id}_{read}.fastq.gz")
    
    # Use dictionary lookup for standard patterns
    key = (source, processing)
    if key in path_patterns:
        path = path_patterns[key]
        if path is None:
            raise FileNotFoundError(f"Could not find {read} file for sample {sample_id} with source={source}, processing={processing}")
        return path
    
    raise FileNotFoundError(f"Could not find {read} file for sample {sample_id} with source={source}, processing={processing}")


def print_sample_summary(config):
    """Print sample detection summary"""
    from scripts.pipeline_utils import get_guide_gex_pairings
    
    IDS = get_sample_ids(config)
    sample_df = load_sample_info(config)
    
    print("=" * 60, file=sys.stderr)
    print("SAMPLE DETECTION SUMMARY", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"Total samples detected: {len(IDS)}", file=sys.stderr)
    
    if not sample_df.empty:
        print("\nPOOL\tGEXSAMP\tGUIDESAMP", file=sys.stderr)
        print("-" * 40, file=sys.stderr)
        for pool in sorted(sample_df['pool'].unique()):
            pool_samples = sample_df[sample_df['pool'] == pool]
            gex_in_pool = pool_samples[pool_samples['sample_type'] == 'gex']['sample_id'].tolist()
            guide_in_pool = pool_samples[pool_samples['sample_type'] == 'guide']['sample_id'].tolist()
            
            # Print each GEX-Guide pair in the pool
            max_samples = max(len(gex_in_pool), len(guide_in_pool))
            for i in range(max_samples):
                # Extract just the sample name for display
                gex_sample = extract_sample(gex_in_pool[i]) if i < len(gex_in_pool) else ""
                guide_sample = extract_sample(guide_in_pool[i]) if i < len(guide_in_pool) else ""
                pool_name = pool if i == 0 else ""  # Only show pool name on first row
                print(f"{pool_name}\t{gex_sample}\t{guide_sample}", file=sys.stderr)
    
    print("=" * 60, file=sys.stderr)


# =============================================================================
# Output generation functions
# =============================================================================

def get_preprocessing_outputs(config, combinations=None, as_dict=False, report_dir="preprocessing_report"):
    """Generate all expected output files for the preprocessing pipeline
    
    This function is used by preprocessing.smk to define all outputs:
    - Read statistics
    - Cell calling results and plots
    - QC metrics (sample, biological sample, well, pool)
    - Saturation analysis data and plots
    - Consolidated metrics and visualizations
    
    Does NOT include:
    - UMAP plots (handled by downstream.smk)
    - Final report DONE file (to avoid circular dependency)
    """
    if combinations is None:
        # Get combinations from config (top level)
        combinations = config['combinations']
        # Convert list of lists to list of tuples
        combinations = [tuple(combo) for combo in combinations]
    
    # Initialize with preprocessing categories only
    outputs_dict = {
        'read_stats': [],
        'cell_calling': [],
        'qc_metrics': [],
        'saturation': [],
        'consolidated': []
    }
    
    sample_df = load_sample_info(config)
    
    # First, iterate through all samples to generate per-sample outputs
    for _, row in sample_df.iterrows():
        pool = row['pool']
        sample_id = row['sample_id']  # e.g., "pool1:gex_1"
        sample_name = extract_sample(sample_id)  # e.g., "gex_1"
        sample_type = row['sample_type']
        
        for source, processing in combinations:
            # Read statistics - use sample_id in directories
            if sample_type == 'gex':
                outputs_dict['read_stats'].append(
                    f"{get_results_path(config=config)}/{sample_id}/qc/{sample_id}_all_{source}_{processing}_read_statistics.tsv"
                )
            else:  # guide
                outputs_dict['read_stats'].append(
                    f"{get_results_path(config=config)}/{sample_id}/qc/{sample_id}_guide_{source}_{processing}_read_statistics.tsv"
                )
            
            # Cell calling and QC metrics for GEX samples only
            if sample_type == 'gex':
                # Cell calling outputs - use sample_id for directory structure
                outputs_dict['cell_calling'].append(
                    f"{get_results_path(config=config)}/{report_dir}/data/per_sample/{source}_{processing}/{sample_id}/cell_calling"
                )
                outputs_dict['cell_calling'].append(
                    f"{get_results_path(config=config)}/{report_dir}/plots/cell_calling/{source}_{processing}/{sample_id}"
                )
                outputs_dict['cell_calling'].append(
                    f"{get_results_path(config=config)}/{report_dir}/plots/cell_calling/{source}_{processing}/{sample_id}.complete"
                )
                
                # QC metrics - all three levels (TSVs serve as sentinels for plot completion)
                for level in ['by_sample.tsv', 'by_biological_sample.tsv', 'by_well.tsv']:
                    outputs_dict['qc_metrics'].append(
                        f"{get_results_path(config=config)}/{report_dir}/data/per_sample/{source}_{processing}/{sample_id}/qc_metrics/{level}"
                    )
    
    # Saturation analysis for main/raw only - iterate through samples
    if ('main', 'raw') in combinations:
        for _, row in sample_df.iterrows():
            pool = row['pool']
            sample_id = row['sample_id']  # e.g., "pool1:gex_1"
            sample_name = extract_sample(sample_id)  # e.g., "gex_1"
            sample_type = row['sample_type']
            
            if sample_type == 'gex':
                # GEX saturation - use sample_id for directory structure
                outputs_dict['saturation'].append(
                    f"{get_results_path(config=config)}/{report_dir}/data/per_sample/main_raw/{sample_id}/saturation/umi_saturation.tsv"
                )
                outputs_dict['saturation'].append(
                    f"{get_results_path(config=config)}/{report_dir}/plots/saturation/main_raw/gex/{sample_id}"
                )
            elif sample_type == 'guide':
                # Guide saturation - use sample_id for directory structure
                outputs_dict['saturation'].append(
                    f"{get_results_path(config=config)}/{report_dir}/data/per_sample/main_raw/{sample_id}/saturation/guide_umi_saturation.tsv"
                )
                outputs_dict['saturation'].append(
                    f"{get_results_path(config=config)}/{report_dir}/plots/saturation/main_raw/guide/{sample_id}"
                )
        
        # Pool statistics - uses the filtered sample_df so pool1000/1001 are excluded when needed
        pools = sample_df['pool'].unique()
        outputs_dict['qc_metrics'].extend(
            [f"{get_results_path(config=config)}/{report_dir}/data/per_pool/main_raw/{pool}/pool_statistics.tsv"
             for pool in pools]
        )
        outputs_dict['consolidated'].append(f"{get_results_path(config=config)}/{report_dir}/data/consolidated/main_raw/by_pool.tsv")
    
    # Consolidated outputs for each combination
    for source, processing in combinations:
        outputs_dict['consolidated'].extend([
            f"{get_results_path(config=config)}/{report_dir}/data/consolidated/{source}_{processing}/all_metrics.tsv",
            f"{get_results_path(config=config)}/{report_dir}/data/consolidated/{source}_{processing}/by_biological_sample.tsv",
            f"{get_results_path(config=config)}/{report_dir}/data/consolidated/{source}_{processing}/by_well.tsv",
            f"{get_results_path(config=config)}/{report_dir}/plots/consolidated_general/{source}_{processing}",
            f"{get_results_path(config=config)}/{report_dir}/plots/consolidated_cell_based/{source}_{processing}",
            f"{get_results_path(config=config)}/{report_dir}/plots/consolidated_{source}_{processing}.complete"
        ])
    
    if as_dict:
        return outputs_dict
    else:
        return [item for sublist in outputs_dict.values() for item in sublist]


def get_downstream_outputs(config, combinations=None, as_dict=False, report_dir="downstream_report"):
    """Generate all expected output files for the downstream analysis pipeline
    
    This function is used by downstream.smk to define all outputs:
    - UMAP plots from standard analyses
    - Clustering results (future)
    - Differential expression results (future)
    
    Only generates outputs when combined_sublibraries is configured.
    """
    if combinations is None:
        # Get combinations from config (top level)
        combinations = config['combinations']
        # Convert list of lists to list of tuples
        combinations = [tuple(combo) for combo in combinations]
    
    # Initialize with downstream categories only
    outputs_dict = {
        'preprocessed': [],
        'pseudobulk': [],
        'umap': [],
        'differential_expression': [],
    }
    
    # Add preprocessed files and UMAPs if combined sublibraries are configured
    if 'combined_sublibraries' in config:
        for source, processing in combinations:
            # Preprocessed h5ad files
            outputs_dict['preprocessed'].append(
                f"{get_results_path(config=config)}/{config['combined_sublibraries']['output_dir']}/preprocessed_{source}_{processing}.h5ad"
            )
            # UMAP plots
            outputs_dict['umap'].append(f"{get_results_path(config=config)}/{report_dir}/plots/umap/{source}_{processing}")
            outputs_dict['umap'].append(f"{get_results_path(config=config)}/{report_dir}/plots/umap/{source}_{processing}.complete")
    
    # Add pseudobulk outputs
    for source, processing in combinations:
        outputs_dict['pseudobulk'].append(
            f"{get_results_path(config=config)}/pseudobulk/{source}_{processing}/pseudobulk_tmm_cpm.tsv"
        )
    
    # Add differential expression outputs if configured
    if 'differential_expression' in config and 'contrasts' in config['differential_expression']:
        for source, processing in combinations:
            outputs_dict['differential_expression'].append(
                f"{get_results_path(config=config)}/differential_expression/{source}_{processing}/analysis.complete"
            )
    
    if as_dict:
        return outputs_dict
    else:
        return [item for sublist in outputs_dict.values() for item in sublist]

