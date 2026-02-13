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



def get_results_path_origin(*parts, config=None):
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
import os

def get_results_path(*parts, config=None):
    """Build results directory paths automatically from analysis_name
    If config['results_dir'] is set, use it as base; otherwise use cwd/results_{analysis_name}.
    """
    if not config:
        raise ValueError("config is required for get_results_path")

    analysis_name = config["analysis_name"]
    results_dir = config.get("results_dir", None)

    if results_dir:
        results_base = os.path.join(results_dir, f"results_{analysis_name}")
    else:
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


def get_raw_data_path(pool, config=None):
    """Get raw data directory path for a pool from sample_info Excel
    
    This function looks up the fastq_dir from sample_info.xlsx for any sample in the given pool.
    All samples in a pool should have the same fastq_dir.
    """
    if config:
        sample_df = load_sample_info(config)
        pool_samples = sample_df[sample_df['pool'] == pool]
        if not pool_samples.empty:
            # All samples in a pool should have the same fastq_dir
            fastq_dir = pool_samples.iloc[0]['fastq_dir']
            return fastq_dir
    
    # Fallback if no config or no samples found
    raise ValueError(f"Could not find fastq_dir for pool {pool} in sample_info.xlsx")


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
    """Load sample information from Excel file"""
    global _sample_info_cache
    if _sample_info_cache is not None:
        return _sample_info_cache
    
    sample_info_file = config["sample_info_file"]
    if not os.path.exists(sample_info_file):
        raise FileNotFoundError(f"Sample info file not found: {sample_info_file}")
    
    # Read the first sheet of the Excel file
    df = pd.read_excel(sample_info_file)
    
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

def find_fastq_file(sample_id, read, source=None, processing=None, config=None, scratch_base=None):
    """Find FASTQ file based on source and processing state
    
    Args:
        sample_id: Full sample ID in format 'pool:sample' (e.g., 'pool1:gex_1')
        read: R1 or R2
        source: main, undetermined, or all
        processing: raw, recovered, or merged
        config: Configuration dictionary
        scratch_base: Base scratch directory (optional, will compute if not provided)
    """
    # Extract pool and sample name from sample_id
    pool = extract_pool(sample_id)
    sample_name = extract_sample(sample_id)
    
    # Use provided scratch_base or compute it
    if scratch_base is None:
        scratch_base = get_scratch_path(config=config)
    
    # Define path patterns for most cases - use sample_name for filenames
    path_patterns = {
        ("main", "raw"): None,  # Special handling below
        ("main", "recovered"): os.path.join(scratch_base, "barcode_recovery", "main", f"{sample_id}_recovered_{read}.fastq.gz"),
        ("undetermined", "raw"): os.path.join(scratch_base, "undetermined_fastqs", f"{sample_id}_{read}.fastq.gz"),
        ("undetermined", "recovered"): os.path.join(scratch_base, "barcode_recovery", "undetermined", f"{sample_id}_recovered_{read}.fastq.gz"),
        ("all", "merged"): os.path.join(scratch_base, "merged_fastqs", f"{sample_id}_all_merged_{read}.fastq.gz")
    }
    
    # Special handling for main/raw
    if source == "main" and processing == "raw":
        if sample_name == "Undetermined":
            # Get the fastq_dir from any sample in this pool
            sample_df = load_sample_info(config)
            pool_samples = sample_df[sample_df['pool'] == pool]
            if not pool_samples.empty:
                fastq_dir = pool_samples.iloc[0]['fastq_dir']
                cmd = f"find {fastq_dir} -name '*Undetermined*{read}*.fastq.gz' | grep -v 'I1\\|I2' | head -1"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                if result.stdout.strip():
                    return result.stdout.strip()
        
        sample_df = load_sample_info(config)
        if not sample_df.empty:
            # Look up by full sample_id (which includes pool)
            sample_row = sample_df[sample_df['sample_id'] == sample_id]
            if not sample_row.empty:
                fastq_dir = sample_row.iloc[0]['fastq_dir']
                # Try pattern with lane info first (e.g., _L001_) - use sample_name for filename
                search_pattern = os.path.join(fastq_dir, f"{sample_name}*_S*_L*_{read}_001.fastq.gz")
                matches = glob.glob(search_pattern)
                if matches:
                    return matches[0]
                # Try pattern without lane info (e.g., gex_1_S1_R1_001.fastq.gz)
                search_pattern = os.path.join(fastq_dir, f"{sample_name}*_S*_{read}_001.fastq.gz")
                matches = glob.glob(search_pattern)
                if matches:
                    return matches[0]
                # Try combined naming (works for both GEX and guide if they share naming)
                search_pattern = os.path.join(fastq_dir, f"{sample_name}*_combined_{read}.fastq.gz")
                matches = glob.glob(search_pattern)
                if matches:
                    return matches[0]
                    
                # Try alternative naming patterns for GEX/gRNA files
                # Pattern: gex_R1_sub1.fastq.gz or gRNA_R1_sub1.fastq.gz
                if '_gex' in sample_name:
                    # Extract the sub number from sample_name (e.g., sub1_gex -> sub1)
                    sub_num = sample_name.replace('_gex', '')
                    search_pattern = os.path.join(fastq_dir, f"GEX_{read}_{sub_num}.fastq.gz")
                    matches = glob.glob(search_pattern)
                    if matches:
                        return matches[0]
                elif '_grna' in sample_name or '_guide' in sample_name:
                    # Extract the sub number from sample_name (e.g., sub1_grna -> sub1)
                    sub_num = sample_name.replace('_grna', '').replace('_guide', '')
                    search_pattern = os.path.join(fastq_dir, f"gRNA_{read}_{sub_num}.fastq.gz")
                    matches = glob.glob(search_pattern)
                    if matches:
                        return matches[0]
    
    # Use dictionary lookup for standard patterns
    key = (source, processing)
    if key in path_patterns:
        return path_patterns[key]
    
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

def get_all_outputs(config, combinations=None, as_dict=False):
    """Generate all expected output files including QC metrics
    
    IMPORTANT: Directory Structure for QC Dashboard
    -----------------------------------------------
    The Streamlit dashboard parses metadata from directory structure following these principles:
    
    1. Top-level categories define plot/data types (e.g., cell_calling, per_cell, saturation)
    2. Second level is always {source}_{processing} (e.g., main_raw, all_merged)
    3. Subsequent levels vary by category but follow consistent patterns within each category
    4. All metadata comes from directory names, not filenames (except for backwards compatibility)
    5. Scale (linear/log) should be a directory level, not part of the filename
    
    When adding new outputs:
    - Keep consistent directory depth within each category
    - Use descriptive directory names that will become filter options
    - Consider if your output fits an existing category before creating a new one
    - If creating a new category, update the dashboard's parse_path_metadata() function
    
    The dashboard automatically creates filters from directory structure, so thoughtful
    organization improves user experience.
    """
    if combinations is None:
        # Get combinations from config
        combinations = config['analysis']['combinations']
        # Convert list of lists to list of tuples
        combinations = [tuple(combo) for combo in combinations]
    
    outputs_dict = {
        'read_stats': [],
        'cell_calling': [],
        'qc_metrics': [],
        'saturation': [],
        'consolidated': [],
        'report': [f"{get_results_path(config=config)}/qc_report/DONE.txt"]
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
                    f"{get_results_path(config=config)}/qc_report/data/per_sample/{source}_{processing}/{sample_id}/cell_calling"
                )
                outputs_dict['cell_calling'].append(
                    f"{get_results_path(config=config)}/qc_report/plots/cell_calling/{source}_{processing}/{sample_id}"
                )
                outputs_dict['cell_calling'].append(
                    f"{get_results_path(config=config)}/qc_report/plots/cell_calling/{source}_{processing}/{sample_id}.complete"
                )
                
                # QC metrics - all three levels
                for level in ['by_sample.tsv', 'by_biological_sample.tsv', 'by_well.tsv']:
                    outputs_dict['qc_metrics'].append(
                        f"{get_results_path(config=config)}/qc_report/data/per_sample/{source}_{processing}/{sample_id}/qc_metrics/{level}"
                    )
                
                # Per-cell plots
                outputs_dict['qc_metrics'].append(
                    f"{get_results_path(config=config)}/qc_report/plots/per_cell/{source}_{processing}/{sample_id}"
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
                    f"{get_results_path(config=config)}/qc_report/data/per_sample/main_raw/{sample_id}/saturation/umi_saturation.tsv"
                )
                outputs_dict['saturation'].append(
                    f"{get_results_path(config=config)}/qc_report/plots/saturation/main_raw/gex/{sample_id}"
                )
            elif sample_type == 'guide':
                # Guide saturation - use sample_id for directory structure
                outputs_dict['saturation'].append(
                    f"{get_results_path(config=config)}/qc_report/data/per_sample/main_raw/{sample_id}/saturation/guide_umi_saturation.tsv"
                )
                outputs_dict['saturation'].append(
                    f"{get_results_path(config=config)}/qc_report/plots/saturation/main_raw/guide/{sample_id}"
                )
        
        # Pool statistics - uses the filtered sample_df so pool1000/1001 are excluded when needed
        pools = sample_df['pool'].unique()
        outputs_dict['qc_metrics'].extend(
            [f"{get_results_path(config=config)}/qc_report/data/per_pool/main_raw/{pool}/pool_statistics.tsv"
             for pool in pools]
        )
        outputs_dict['consolidated'].append(f"{get_results_path(config=config)}/qc_report/data/consolidated/main_raw/by_pool.tsv")
    
    # Consolidated outputs for each combination
    for source, processing in combinations:
        outputs_dict['consolidated'].extend([
            f"{get_results_path(config=config)}/qc_report/data/consolidated/{source}_{processing}/all_metrics.tsv",
            f"{get_results_path(config=config)}/qc_report/data/consolidated/{source}_{processing}/by_biological_sample.tsv",
            f"{get_results_path(config=config)}/qc_report/data/consolidated/{source}_{processing}/by_well.tsv",
            f"{get_results_path(config=config)}/qc_report/plots/consolidated_general/{source}_{processing}",
            f"{get_results_path(config=config)}/qc_report/plots/consolidated_cell_based/{source}_{processing}",
            f"{get_results_path(config=config)}/qc_report/plots/consolidated_{source}_{processing}.complete"
        ])
    
    if as_dict:
        return outputs_dict
    else:
        return [item for sublist in outputs_dict.values() for item in sublist]