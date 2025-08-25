#!/usr/bin/env python3
"""
Fast package QC outputs for laptop viewing with Streamlit dashboard
Uses direct tar creation without intermediate file copying
Modified to work with directory structure instead of JSON metadata
"""

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from datetime import datetime
import sys
import re
import time
import tarfile
import io

def parse_args():
    parser = argparse.ArgumentParser(description="Package QC outputs for laptop viewing")
    parser.add_argument('--qc-report-dir', required=True,
                        help='Path to QC report directory (contains plots/ and data/ subdirectories)')
    parser.add_argument('--sample-info', required=True,
                        help='Path to sample_info.tsv')
    parser.add_argument('--output-archive', required=True,
                        help='Output tar.gz file path')
    parser.add_argument('--dashboard-script', 
                        default='scripts/streamlit_qc_dashboard_no_json.py',
                        help='Path to Streamlit dashboard script')
    parser.add_argument('--per-cell-method-filter', default=None,
                        help='Filter per_cell plots to only include this method (e.g., BarcodeRanks_Knee)')
    parser.add_argument('--guide-cutoff-filter', default=None,
                        help='Comma-separated list of guide cutoffs to include (e.g., 1,2)')
    parser.add_argument('--threads', type=int, default=8,
                        help='Number of threads for pigz compression')
    return parser.parse_args()

def should_include_file(file_path, per_cell_method_filter, allowed_cutoffs):
    """Check if a file should be included based on filters"""
    path_str = str(file_path)
    
    # Skip JSON metadata files (no longer needed)
    if path_str.endswith('.json') and 'metadata' in path_str:
        return False, 'json'
    
    # Check per_cell method filter
    if per_cell_method_filter and 'per_cell' in path_str:
        if per_cell_method_filter not in path_str:
            return False, 'method'
    
    # Check guide cutoff filter
    if allowed_cutoffs and 'cutoff' in path_str:
        # Look for cutoff patterns
        cutoff_match = re.search(r'cutoff(\d+)', path_str)
        if cutoff_match:
            cutoff = cutoff_match.group(1)
            if cutoff not in allowed_cutoffs:
                return False, 'cutoff'
    
    return True, None

def create_text_file(arcname, content):
    """Create a TarInfo and fileobj for text content"""
    encoded = content.encode('utf-8')
    info = tarfile.TarInfo(name=arcname)
    info.size = len(encoded)
    info.mtime = time.time()
    fileobj = io.BytesIO(encoded)
    return info, fileobj

def create_archive_directly(args):
    """Create tar.gz archive directly without intermediate copying"""
    
    # Parse filters
    per_cell_method_filter = args.per_cell_method_filter
    allowed_cutoffs = set(args.guide_cutoff_filter.split(',')) if args.guide_cutoff_filter else None
    
    print(f"\nCreating archive: {args.output_archive}")
    if per_cell_method_filter:
        print(f"  Filtering per_cell plots to method: {per_cell_method_filter}")
    if allowed_cutoffs:
        print(f"  Filtering guide plots to cutoffs: {allowed_cutoffs}")
    
    # Scan QC report directory for all plots and data
    qc_report_dir = Path(args.qc_report_dir)
    if not qc_report_dir.exists():
        raise FileNotFoundError(f"QC report directory not found: {qc_report_dir}")
    
    plot_files = []
    data_files = []
    
    # Find all files in plots/ and data/ subdirectories
    plots_dir = qc_report_dir / 'plots'
    data_dir = qc_report_dir / 'data'
    
    if plots_dir.exists():
        for plot_file in plots_dir.rglob('*'):
            if plot_file.is_file() and not plot_file.name.startswith('.'):
                plot_files.append(str(plot_file))
    
    if data_dir.exists():
        for data_file in data_dir.rglob('*'):
            if data_file.is_file() and not data_file.name.startswith('.'):
                data_files.append(str(data_file))
    
    # Count files before filtering
    total_plots = len(plot_files)
    total_data = len(data_files)
    
    # Apply filters
    included_plots = []
    excluded_plots = {'method': 0, 'cutoff': 0, 'json': 0}
    
    for plot_file in plot_files:
        include, reason = should_include_file(plot_file, per_cell_method_filter, allowed_cutoffs)
        if include:
            included_plots.append(plot_file)
        else:
            excluded_plots[reason] += 1
    
    print(f"\nFile counts:")
    print(f"  Plots: {len(included_plots)}/{total_plots} included")
    if excluded_plots['method'] > 0:
        print(f"    - {excluded_plots['method']} excluded by method filter")
    if excluded_plots['cutoff'] > 0:
        print(f"    - {excluded_plots['cutoff']} excluded by cutoff filter")
    if excluded_plots['json'] > 0:
        print(f"    - {excluded_plots['json']} JSON metadata files excluded")
    print(f"  Data files: {total_data}")
    
    # Create output directory
    output_path = Path(args.output_archive)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Extract the base folder name from the archive name
    # e.g., qc_dashboard_downstream_20241115_143022.tar.gz -> qc_dashboard_downstream_20241115_143022
    folder_name = output_path.stem  # Removes .gz
    if folder_name.endswith('.tar'):
        folder_name = folder_name[:-4]  # Remove .tar if present
    
    # Create archive
    start_time = time.time()
    
    try:
        # First create uncompressed tar
        temp_tar = output_path.with_suffix('')  # Remove .gz
        
        # Create tar file
        with tarfile.open(temp_tar, 'w') as tar:
            # Add static files
            print("  Adding dashboard files...")
            
            # Add dashboard script
            tar.add(args.dashboard_script, arcname=f'{folder_name}/streamlit_qc_dashboard.py')
            
            # Add sample info
            tar.add(args.sample_info, arcname=f'{folder_name}/sample_info.tsv')
            
            # Add metric glossary if it exists
            metric_glossary_path = Path(args.dashboard_script).parent / 'metric_glossary.tsv'
            if metric_glossary_path.exists():
                tar.add(metric_glossary_path, arcname=f'{folder_name}/metric_glossary.tsv')
            
            # Create and add run script
            run_script_content = """#!/bin/bash
# QC Dashboard Launcher

echo "Starting QC Dashboard..."
echo "Make sure you have streamlit installed: pip install streamlit pandas numpy"
echo ""

# Run streamlit
streamlit run streamlit_qc_dashboard.py --server.port 8501

echo "Dashboard stopped."
"""
            info, fileobj = create_text_file(f'{folder_name}/run_dashboard.sh', run_script_content)
            info.mode = 0o755  # Make executable
            tar.addfile(tarinfo=info, fileobj=fileobj)
            
            # Create and add README
            readme_content = f"""QC Dashboard Package (Directory Structure Version)
Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

REQUIREMENTS:
- Python 3.7+
- Install dependencies: pip install streamlit pandas numpy openpyxl

TO RUN:
1. Extract this archive
2. cd {folder_name}
3. ./run_dashboard.sh
   (or: streamlit run streamlit_qc_dashboard.py)

The dashboard will open in your browser at http://localhost:8501

FILES INCLUDED:
- streamlit_qc_dashboard.py: Interactive dashboard
- sample_info.tsv: Sample metadata
- metric_glossary.tsv: Metric definitions and descriptions
- plots/: All QC plots organized by category
- data/: Tabular QC data

NOTES:
- This version uses directory structure instead of JSON metadata
- Plot metadata is inferred from file paths and names
- Metric descriptions are defined in metric_glossary.tsv
"""
            info, fileobj = create_text_file(f'{folder_name}/README.txt', readme_content)
            tar.addfile(tarinfo=info, fileobj=fileobj)
            
            # Add plots
            print(f"  Adding {len(included_plots)} plot files...")
            for i, plot_file in enumerate(included_plots, 1):
                if i % 100 == 0:
                    print(f"    Progress: {i}/{len(included_plots)} plots")
                
                # Convert path to be relative to analysis_results
                plot_path = Path(plot_file)
                if 'qc_report' in plot_path.parts:
                    # Find index of 'qc_report' and take everything after it
                    idx = plot_path.parts.index('qc_report')
                    arcname = f'{folder_name}/' + '/'.join(plot_path.parts[idx+1:])
                else:
                    # Fallback to preserving structure after 'plots'
                    if 'plots' in plot_path.parts:
                        idx = plot_path.parts.index('plots')
                        arcname = f'{folder_name}/' + '/'.join(plot_path.parts[idx:])
                    else:
                        arcname = f'{folder_name}/plots/' + plot_path.name
                
                tar.add(plot_file, arcname=arcname)
            
            # Add data files
            print(f"  Adding {total_data} data files...")
            for i, data_file in enumerate(data_files, 1):
                if i % 50 == 0:
                    print(f"    Progress: {i}/{total_data} data files")
                
                # Convert path to be relative
                data_path = Path(data_file)
                if 'qc_report' in data_path.parts:
                    idx = data_path.parts.index('qc_report')
                    arcname = f'{folder_name}/' + '/'.join(data_path.parts[idx+1:])
                else:
                    if 'data' in data_path.parts:
                        idx = data_path.parts.index('data')
                        arcname = f'{folder_name}/' + '/'.join(data_path.parts[idx:])
                    else:
                        arcname = f'{folder_name}/data/' + data_path.name
                
                tar.add(data_file, arcname=arcname)
        
        # Compress the tar file
        print(f"  Compressing archive with pigz ({args.threads} threads)...")
        subprocess.run(['pigz', '-f', '-p', str(args.threads), str(temp_tar)], check=True)
        # pigz -f creates temp_tar.gz which is exactly our output_path
        
    except Exception as e:
        print(f"Error creating archive: {e}")
        # Cleanup
        if temp_tar.exists():
            temp_tar.unlink()
        if temp_tar.with_suffix('.gz').exists():
            temp_tar.with_suffix('.gz').unlink()
        raise
    
    elapsed = time.time() - start_time
    size_mb = output_path.stat().st_size / (1024 * 1024)
    
    print(f"\n✅ Archive created successfully!")
    print(f"   File: {output_path}")
    print(f"   Size: {size_mb:.1f} MB")
    print(f"   Time: {elapsed:.1f} seconds")
    
    return output_path

def main():
    args = parse_args()
    
    print("=== QC Dashboard Packaging (No JSON Version) ===")
    print(f"QC report directory: {args.qc_report_dir}")
    print(f"Output archive: {args.output_archive}")
    
    # Create archive directly
    archive_path = create_archive_directly(args)
    
    print("\n📦 Package ready for distribution!")
    print(f"   Transfer {archive_path} to your laptop and extract it")
    print("   Then run: ./run_dashboard.sh")

if __name__ == "__main__":
    main()