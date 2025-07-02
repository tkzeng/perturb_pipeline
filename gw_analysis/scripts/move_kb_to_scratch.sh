#!/bin/bash
set -euo pipefail

# Script to move kb outputs from analysis_results to scratch
echo "Moving kb outputs to scratch directory..."

# Find all sample directories with kb outputs
for sample_dir in ../analysis_results/*/; do
    if [ -d "$sample_dir" ]; then
        sample_name=$(basename "$sample_dir")
        
        # Skip non-sample directories
        if [[ ! "$sample_name" =~ ^pool[0-9]+.*$ ]]; then
            continue
        fi
        
        # Create target directory in scratch
        target_dir="/scratch/users/tkzeng/gw_analysis/$sample_name"
        mkdir -p "$target_dir"
        
        # Move all kb_* directories for this sample
        for kb_dir in "$sample_dir"kb_*; do
            if [ -d "$kb_dir" ]; then
                kb_name=$(basename "$kb_dir")
                echo "Moving $sample_name/$kb_name to scratch..."
                
                # Move the directory
                mv "$kb_dir" "$target_dir/"
            fi
        done
    fi
done

echo "Done moving kb outputs to scratch!"

# Show disk usage summary
echo -e "\nDisk usage summary:"
echo "Scratch usage:"
du -sh /scratch/users/tkzeng/gw_analysis/ 2>/dev/null || echo "Unable to calculate"
echo "Analysis results remaining:"
du -sh ../analysis_results/ 2>/dev/null || echo "Unable to calculate"