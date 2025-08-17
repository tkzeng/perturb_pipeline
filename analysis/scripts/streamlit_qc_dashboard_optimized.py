#!/usr/bin/env python3
"""
Streamlit QC Dashboard - Optimized for Large Scale

Handles 10,000+ plots efficiently using pre-computed indices and pagination.
"""

import streamlit as st
import pandas as pd
from pathlib import Path
import re
from collections import defaultdict
import numpy as np
from typing import Dict, List, Set, Tuple
from datetime import datetime
import os


# Page configuration
st.set_page_config(
    page_title="QC Dashboard",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS - removed plot container styling since we're using native Streamlit components

# Constants
PLOTS_PER_PAGE = 60  # Divisible by 3 for grid layout
MAX_FILTER_PREVIEW = 100  # Max items to show in filter dropdowns

def snake_to_title(text):
    """Convert snake_case to Title Case"""
    return ' '.join(word.capitalize() for word in text.split('_'))

def parse_path_metadata(file_path, base_path):
    """Extract metadata from file path relative to base path"""
    try:
        rel_path = Path(file_path).relative_to(base_path)
        parts = list(rel_path.parts[:-1])  # Exclude filename
        
        metadata = {
            'path_parts': parts,
            'depth': len(parts),
            'full_path': str(file_path),
            'rel_path': str(rel_path),
            'filename': rel_path.name,
            'stem': rel_path.stem
        }
        
        # Parse based on known directory structure
        if len(parts) > 0:
            category = parts[0]
            metadata['category'] = category
            
            # Parse source_processing from parts[1] if it exists
            if len(parts) > 1:
                source_proc = parts[1]
                if '_' in source_proc:
                    source, processing = source_proc.split('_', 1)
                    metadata['source'] = source
                    metadata['processing'] = processing
                else:
                    metadata['source_processing'] = source_proc
            
            # Parse remaining parts based on category with new structure
            if category == 'consolidated_general':
                # Structure: consolidated_general/{source}_{processing}/{stratify_by}/{metric}/{scale}/plot.png
                if len(parts) > 2:
                    metadata['stratify_by'] = parts[2]
                if len(parts) > 3:
                    metadata['metric'] = parts[3]
                if len(parts) > 4:
                    metadata['scale'] = parts[4]
                        
            elif category == 'consolidated_cell_based':
                # Structure: consolidated_cell_based/{source}_{processing}/{stratify_by}/{metric}/{method}/{scale}/plot.png
                if len(parts) > 2:
                    metadata['stratify_by'] = parts[2]
                if len(parts) > 3:
                    metadata['metric'] = parts[3]
                if len(parts) > 4:
                    metadata['method'] = parts[4]
                if len(parts) > 5:
                    metadata['scale'] = parts[5]
                    
            elif category == 'per_cell':
                # Structure: per_cell/{source}_{processing}/{sample_id}/{metric}/{method}/{stratify_by}/{scale}/plot.png
                # where sample_id is {pool}:{sample}
                if len(parts) > 2:
                    sample_id = parts[2]
                    metadata['sample_id'] = sample_id
                    # Split pool:sample format
                    if ':' in sample_id:
                        pool, sample = sample_id.split(':', 1)
                        metadata['pool'] = pool
                        metadata['sample'] = sample
                if len(parts) > 3:
                    metadata['metric'] = parts[3]
                if len(parts) > 4:
                    metadata['method'] = parts[4]
                if len(parts) > 5:
                    metadata['stratify_by'] = parts[5]
                if len(parts) > 6:
                    metadata['scale'] = parts[6]
                    
            elif category == 'saturation':
                # Structure: saturation/{source}_{processing}/{sample_type}/{sample_id}/{metric}/{scale}/plot.png
                # where sample_id is {pool}:{sample}
                if len(parts) > 1:
                    source_proc = parts[1]
                    if '_' in source_proc:
                        source, processing = source_proc.split('_', 1)
                        metadata['source'] = source
                        metadata['processing'] = processing
                    else:
                        metadata['source_processing'] = source_proc
                if len(parts) > 2:
                    metadata['sample_type'] = parts[2]  # gex or guide
                if len(parts) > 3:
                    sample_id = parts[3]
                    metadata['sample_id'] = sample_id
                    # Split pool:sample format
                    if ':' in sample_id:
                        pool, sample = sample_id.split(':', 1)
                        metadata['pool'] = pool
                        metadata['sample'] = sample
                if len(parts) > 4:
                    metadata['metric'] = parts[4]  # umi_saturation or guide_umi_saturation
                if len(parts) > 5:
                    metadata['scale'] = parts[5]  # Always linear for saturation plots
                    
            elif category == 'cell_calling':
                # Structure: cell_calling/{source}_{processing}/{sample_id}/{metric}/{scale}/plot.png
                # where sample_id is {pool}:{sample}
                if len(parts) > 2:
                    sample_id = parts[2]
                    metadata['sample_id'] = sample_id
                    # Split pool:sample format
                    if ':' in sample_id:
                        pool, sample = sample_id.split(':', 1)
                        metadata['pool'] = pool
                        metadata['sample'] = sample
                if len(parts) > 3:
                    metadata['metric'] = parts[3]
                if len(parts) > 4:
                    metadata['scale'] = parts[4]  # Always linear for cell calling plots
            
            elif category.startswith('umap'):
                # Fixed structure: umap/{source}_{processing}/{umap_subset}/{de_subset}/{metric}/plot.png
                
                metadata['umap_subset'] = parts[2]  # e.g., "all_cells_umap", "low_mito", "high_mito"
                metadata['de_subset'] = parts[3]    # e.g., "all_cells_scores", "low_mito", "high_mito" 
                metadata['metric'] = parts[4]       # e.g., "pct_counts_mt", "guide_moi_score_pos"
            
            # Handle any other categories generically
            else:
                for i in range(2, len(parts)):
                    metadata[f'level_{i}'] = parts[i]
        
        return metadata
        
    except Exception as e:
        # If parsing fails, return basic metadata with error
        return {
            'full_path': str(file_path),
            'filename': Path(file_path).name,
            'stem': Path(file_path).stem,
            'parse_error': str(e),
            'path_parts': [],
            'depth': 0
        }


@st.cache_data
def scan_and_index_plots(base_path="plots"):
    """Scan plots directory and build optimized indices"""
    plots = []
    metadata_indices = {}
    metadata_values = {}
    # Category-specific indices for tab-specific filters
    category_metadata_indices = {}
    category_metadata_values = {}
    base_path = Path(base_path)
    
    if not base_path.exists():
        return plots, metadata_indices, metadata_values, category_metadata_indices, category_metadata_values
    
    # Find all PNG files and build indices
    plot_id = 0
    for png_file in base_path.rglob("*.png"):
        # Skip hidden files and thumbnails
        if any(part.startswith('.') for part in png_file.parts):
            continue
            
        # Extract path metadata - this now contains everything we need
        path_meta = parse_path_metadata(png_file, base_path)
        
        # Get file modification time
        try:
            mtime = os.path.getmtime(png_file)
            mod_time = datetime.fromtimestamp(mtime)
            path_meta['mod_time'] = mod_time
            path_meta['mod_time_str'] = mod_time.strftime("%Y-%m-%d %H:%M:%S")
        except:
            path_meta['mod_time'] = None
            path_meta['mod_time_str'] = "Unknown"
        
        # All metadata is now in the path
        plot_info = {**path_meta, 'id': plot_id}
        
        # Create display name - always have one even if parsing failed
        plot_info['display_name'] = snake_to_title(plot_info.get('stem', Path(png_file).stem))
        
        # If there's no category due to parse error, mark it as 'uncategorized'
        if 'category' not in plot_info and not plot_info.get('parse_error'):
            plot_info['category'] = 'uncategorized'
        elif plot_info.get('parse_error'):
            plot_info['category'] = 'parse_errors'
        
        # Get the category for this plot
        category = plot_info.get('category', 'uncategorized')
        
        # Build global indices - include parse_error plots too
        for key, value in plot_info.items():
            if key not in ['full_path', 'rel_path', 'filename', 'display_name', 
                          'path_parts', 'depth', 'stem', 'id', 'parse_error', 
                          'filename_parse_error'] and value is not None:
                # Initialize structures if needed
                if key not in metadata_indices:
                    metadata_indices[key] = {}
                if value not in metadata_indices[key]:
                    metadata_indices[key][value] = set()
                if key not in metadata_values:
                    metadata_values[key] = set()
                    
                metadata_indices[key][value].add(plot_id)
                metadata_values[key].add(value)
                
                # Also add to category-specific indices
                if category not in category_metadata_indices:
                    category_metadata_indices[category] = {}
                if category not in category_metadata_values:
                    category_metadata_values[category] = {}
                if key not in category_metadata_indices[category]:
                    category_metadata_indices[category][key] = {}
                if value not in category_metadata_indices[category][key]:
                    category_metadata_indices[category][key][value] = set()
                if key not in category_metadata_values[category]:
                    category_metadata_values[category][key] = set()
                    
                category_metadata_indices[category][key][value].add(plot_id)
                category_metadata_values[category][key].add(value)
        
        # Special handling for parse errors - add them to a special category
        if 'parse_error' in plot_info or 'filename_parse_error' in plot_info:
            if 'has_error' not in metadata_indices:
                metadata_indices['has_error'] = {}
            if 'true' not in metadata_indices['has_error']:
                metadata_indices['has_error']['true'] = set()
            if 'has_error' not in metadata_values:
                metadata_values['has_error'] = set()
                
            metadata_indices['has_error']['true'].add(plot_id)
            metadata_values['has_error'].add('true')
        
        plots.append(plot_info)
        plot_id += 1
    
    # Convert sets to lists for serialization
    for key in metadata_indices:
        for value in metadata_indices[key]:
            metadata_indices[key][value] = list(metadata_indices[key][value])
        metadata_values[key] = sorted(list(metadata_values[key]))
    
    # Convert category-specific sets to lists
    for cat in category_metadata_indices:
        for key in category_metadata_indices[cat]:
            for value in category_metadata_indices[cat][key]:
                category_metadata_indices[cat][key][value] = list(category_metadata_indices[cat][key][value])
            category_metadata_values[cat][key] = sorted(list(category_metadata_values[cat][key]))
    
    return plots, metadata_indices, metadata_values, category_metadata_indices, category_metadata_values


def get_filtered_indices(filters: Dict[str, str], metadata_indices: Dict) -> Set[int]:
    """Get plot indices matching all filters using set operations"""
    if not filters:
        return None  # No filtering needed
    
    # Start with all plots that match the first filter
    result_set = None
    
    for key, value in filters.items():
        if key in metadata_indices and value in metadata_indices[key]:
            indices = set(metadata_indices[key][value])
            if result_set is None:
                result_set = indices
            else:
                result_set = result_set.intersection(indices)
    
    return result_set

# Removed get_available_values function - no longer needed for independent filters

# Define internal fields globally to ensure consistency
INTERNAL_FIELDS = {
    'category', 'id', 'display_name', 'full_path', 'rel_path', 
    'path_parts', 'filename', 'stem', 'depth', 'parse_error',
    'mod_time', 'mod_time_str', 'table_name', 'size_bytes', 
    'size_mb', 'table_type', 'sample_id'  # Exclude sample_id since pool and sample are already available
}

def create_tab_filters(metadata_indices: Dict, 
                      metadata_values: Dict, key_prefix: str) -> Tuple[Dict, Set[int]]:
    """Create filters at the top of each tab - shows ALL available filters"""
    filters = {}
    
    # Get metadata keys (excluding internal fields)
    available_keys = [k for k in metadata_values.keys() 
                     if k not in INTERNAL_FIELDS and len(metadata_values[k]) > 1]
    
    if not available_keys:
        st.info("No filters available for this category")
        return filters, None
    
    # Reorder keys to put 'metric' first if it exists
    if 'metric' in available_keys:
        available_keys.remove('metric')
        available_keys.insert(0, 'metric')
    
    # Sort remaining keys alphabetically but put level_X at the end
    regular_keys = [k for k in available_keys if not k.startswith('level_')]
    level_keys = sorted([k for k in available_keys if k.startswith('level_')], 
                       key=lambda x: int(x.split('_')[1]) if len(x.split('_')) > 1 else 0)
    available_keys = regular_keys + level_keys
    
    # Create expandable filter section
    with st.expander(f"Filters ({len(available_keys)} available)", expanded=True):
        # Show filters in rows of 4 columns
        filters_per_row = 4
        num_rows = (len(available_keys) + filters_per_row - 1) // filters_per_row
        
        # Clear all button at the top
        col1, col2, col3, col4 = st.columns([1, 1, 1, 1])
        with col1:
            if st.button("Clear All Filters", key=f"{key_prefix}_clear"):
                st.rerun()
        
        # Create filter rows
        for row in range(num_rows):
            cols = st.columns(filters_per_row)
            for col_idx in range(filters_per_row):
                filter_idx = row * filters_per_row + col_idx
                if filter_idx < len(available_keys):
                    key = available_keys[filter_idx]
                    values = metadata_values[key]
                    
                    if values:
                        with cols[col_idx]:
                            display_key = snake_to_title(key)
                            # Show number of options in help text
                            selected = st.selectbox(
                                display_key,
                                ['All'] + sorted(list(values))[:MAX_FILTER_PREVIEW],
                                key=f"{key_prefix}_{key}",
                                help=f"{len(values)} options available"
                            )
                            
                            if selected != 'All':
                                filters[key] = selected
        
        # Show active filters summary
        if filters:
            st.markdown("**Active filters:** " + ", ".join([f"{snake_to_title(k)}: {v}" for k, v in filters.items()]))
    
    # Apply all filters with AND logic to get final indices
    if not filters:
        return filters, None
    
    # Start with first filter
    result_indices = None
    for key, value in filters.items():
        if key in metadata_indices and value in metadata_indices[key]:
            indices = set(metadata_indices[key][value])
            if result_indices is None:
                result_indices = indices
            else:
                result_indices = result_indices.intersection(indices)
    
    return filters, result_indices

def display_plots_paginated(plots: List, filtered_indices: Set[int], 
                           metadata_values: Dict, key_prefix: str):
    """Display plots with pagination and grouping controls"""
    # Get filtered plots
    if filtered_indices is not None:
        filtered_plots = [plots[i] for i in sorted(filtered_indices)]
    else:
        filtered_plots = plots
    
    total_plots = len(filtered_plots)
    
    if total_plots == 0:
        st.info("No plots match the current filters")
        return
    
    # Pagination controls
    col1, col2, col3 = st.columns([2, 4, 2])
    with col1:
        page = st.number_input("Page", min_value=1, 
                              max_value=max(1, (total_plots - 1) // PLOTS_PER_PAGE + 1),
                              value=1, key=f"{key_prefix}_plots_page")
    
    with col2:
        st.write(f"Showing {min((page-1)*PLOTS_PER_PAGE + 1, total_plots)}-"
                f"{min(page*PLOTS_PER_PAGE, total_plots)} of {total_plots} plots")
    
    # Grouping controls
    st.markdown("---")
    col1, col2, col3 = st.columns([2, 2, 6])
    
    # Extract metadata keys from current page plots
    start_idx = (page - 1) * PLOTS_PER_PAGE
    end_idx = min(start_idx + PLOTS_PER_PAGE, total_plots)
    page_plots = filtered_plots[start_idx:end_idx]
    
    available_group_keys = set()
    for plot in page_plots:
        available_group_keys.update(k for k, v in plot.items() 
                                  if k not in INTERNAL_FIELDS
                                  and v is not None)
    
    # Primary grouping
    with col1:
        group_options = ['None'] + sorted(list(available_group_keys))
        primary_group = st.selectbox("Group by:", group_options, 
                                     key=f"{key_prefix}_plots_primary_group")
    
    # Secondary grouping
    with col2:
        if primary_group != 'None':
            secondary_options = ['None'] + [k for k in sorted(available_group_keys) if k != primary_group]
            secondary_group = st.selectbox("Then by:", secondary_options,
                                         key=f"{key_prefix}_plots_secondary_group")
        else:
            secondary_group = 'None'
    
    # Display plots for current page
    display_grouped_plots(page_plots, primary_group, secondary_group)

def display_grouped_plots(plots: List, primary_group: str, secondary_group: str):
    """Display plots with grouping"""
    if primary_group == 'None':
        # No grouping - display flat grid
        # Create new columns for each row to maintain alignment
        for i in range(0, len(plots), 3):
            cols = st.columns(3)
            for j in range(3):
                if i + j < len(plots):
                    with cols[j]:
                        display_single_plot(plots[i + j])
    
    elif secondary_group == 'None':
        # Single-level grouping
        groups = defaultdict(list)
        for plot in plots:
            group_value = plot.get(primary_group, 'Other')
            groups[group_value].append(plot)
        
        # Display each group
        for group_name in sorted(groups.keys()):
            st.subheader(f"{snake_to_title(primary_group)}: {group_name}")
            group_plots = groups[group_name]
            # Create new columns for each row
            for i in range(0, len(group_plots), 3):
                cols = st.columns(3)
                for j in range(3):
                    if i + j < len(group_plots):
                        with cols[j]:
                            display_single_plot(group_plots[i + j])
            st.markdown("---")
    
    else:
        # Two-level grouping
        primary_groups = defaultdict(lambda: defaultdict(list))
        for plot in plots:
            primary_value = plot.get(primary_group, 'Other')
            secondary_value = plot.get(secondary_group, 'Other')
            primary_groups[primary_value][secondary_value].append(plot)
        
        # Display nested groups
        for primary_name in sorted(primary_groups.keys()):
            st.subheader(f"{snake_to_title(primary_group)}: {primary_name}")
            
            for secondary_name in sorted(primary_groups[primary_name].keys()):
                if len(primary_groups[primary_name]) > 1:
                    st.markdown(f"**{snake_to_title(secondary_group)}: {secondary_name}**")
                
                secondary_plots = primary_groups[primary_name][secondary_name]
                # Create new columns for each row
                for i in range(0, len(secondary_plots), 3):
                    cols = st.columns(3)
                    for j in range(3):
                        if i + j < len(secondary_plots):
                            with cols[j]:
                                display_single_plot(secondary_plots[i + j])
                
            st.markdown("---")

def display_single_plot(plot):
    """Display a single plot with lazy loading"""
    with st.container():
        # Show parse errors prominently
        if 'parse_error' in plot:
            st.error(f"Path parse error: {plot['parse_error']}")
        if 'filename_parse_error' in plot:
            st.warning(f"Filename parse error: {plot['filename_parse_error']}")
        
        # Build metadata display - show all available metadata
        metadata_parts = []
        # First show standard fields including metric
        for key in ['pool', 'sample', 'source', 'processing', 'metric', 
                    'stratify_by', 'sample_type', 'method', 'scale']:
            if key in plot and plot[key]:
                metadata_parts.append(f"{snake_to_title(key)}: {plot[key]}")
        
        # Then show any level_X fields
        for key in sorted(plot.keys()):
            if key.startswith('level_') and plot[key]:
                metadata_parts.append(f"{snake_to_title(key)}: {plot[key]}")
        
        
        if metadata_parts:
            st.caption(" | ".join(metadata_parts))
        
        # Display image with lazy loading
        try:
            st.image(plot['full_path'], use_container_width=True)
        except Exception as e:
            st.error(f"Error loading image: {str(e)}")
        
        # Show file info: path and modification time
        file_info = f"{plot['full_path']}\nModified: {plot.get('mod_time_str', 'Unknown')}"
        st.code(file_info, language=None)
        
        # Add some spacing between plots
        st.markdown("---")


def load_metric_glossary():
    """Load metric glossary from TSV file."""
    glossary_path = Path(__file__).parent / "metric_glossary.tsv"
    
    if glossary_path.exists():
        return pd.read_csv(glossary_path, sep='\t')
    else:
        st.error(f"Metric glossary file not found: {glossary_path}")
        return pd.DataFrame()  # Return empty dataframe



def main():
    
    # Scan and index plots and tables
    with st.spinner("Building indices..."):
        all_plots, plot_metadata_indices, plot_metadata_values, category_plot_indices, category_plot_values = scan_and_index_plots()
    
    if not all_plots:
        st.error("No plots found. Please check that the 'plots' directory exists.")
        return
    
    # Get plot categories
    plot_categories = sorted(plot_metadata_values.get('category', ['unknown']))
    
    # Put parse_errors first if it exists
    if 'parse_errors' in plot_categories:
        plot_categories.remove('parse_errors')
        plot_categories.insert(0, 'parse_errors')
    
    # Create category options - plot categories + data tables
    category_options = []
    category_counts = {}
    
    for cat in plot_categories:
        if cat in plot_metadata_indices.get('category', {}):
            count = len(plot_metadata_indices['category'][cat])
            category_counts[cat] = count
            # Special styling for parse errors
            if cat == 'parse_errors':
                category_options.append(f"Parse Errors ({count})")
            else:
                category_options.append(f"{snake_to_title(cat)} ({count})")
    
    # Add metric glossary option
    category_options.append("Metric Glossary")
    
    # Sidebar navigation
    with st.sidebar:
        st.title("QC Report")
        st.markdown("---")
        
        # Category selector
        selected_category_display = st.radio(
            "Select Category:",
            category_options,
            key="category_selector"
        )
        
        # Extract the actual category name from display name
        if selected_category_display.startswith("Parse Errors"):
            selected_category = "parse_errors"
        elif selected_category_display == "Metric Glossary":
            selected_category = "metric_glossary"
        else:
            # Extract category name by removing count suffix
            selected_category = None
            for cat in plot_categories:
                display_name = f"{snake_to_title(cat)} ({category_counts.get(cat, 0)})"
                if selected_category_display == display_name:
                    selected_category = cat
                    break
        
        st.markdown("---")
        st.caption(f"Total plots: {len(all_plots)}")
    
    # Main content area
    
    # Display selected category
    if selected_category == "metric_glossary":
        # Metric glossary section
        st.subheader("QC Metric Glossary")
        st.info("Reference guide for all quality control metrics used in this pipeline")
        
        # Load glossary
        glossary_df = load_metric_glossary()
        
        if not glossary_df.empty:
            # Add filter by category
            col1, col2, col3 = st.columns([2, 2, 6])
            with col1:
                categories = ['All'] + sorted(glossary_df['Category'].unique().tolist())
                selected_cat = st.selectbox("Filter by category:", categories)
            
            # Filter dataframe if category selected
            if selected_cat != 'All':
                display_df = glossary_df[glossary_df['Category'] == selected_cat]
            else:
                display_df = glossary_df
            
            # Display as interactive dataframe
            st.dataframe(
                display_df,
                use_container_width=True,
                hide_index=True,
                column_config={
                    "Category": st.column_config.TextColumn("Category", width="small"),
                    "Metric": st.column_config.TextColumn("Metric", width="medium"),
                    "Calculation": st.column_config.TextColumn("Calculation", width="large"),
                    "Description": st.column_config.TextColumn("Description", width="large"),
                }
            )
            
            # Add download button
            tsv_data = display_df.to_csv(sep='\t', index=False)
            st.download_button(
                label="Download glossary as TSV",
                data=tsv_data,
                file_name="qc_metric_glossary.tsv",
                mime="text/tab-separated-values"
            )
    
    else:
        # Plot category section
        if selected_category and 'category' in plot_metadata_indices and selected_category in plot_metadata_indices['category']:
            category_indices = set(plot_metadata_indices['category'][selected_category])
            
            # Create filters using category-specific indices
            cat_indices = category_plot_indices.get(selected_category, {})
            cat_values = category_plot_values.get(selected_category, {})
            
            filters, filtered_indices = create_tab_filters(
                cat_indices, cat_values, f"cat_{selected_category}"
            )
            
            # Add separator between filters and plots
            st.markdown("---")
            
            # Apply category filter to results
            if filtered_indices is not None:
                filtered_indices = filtered_indices.intersection(category_indices)
            else:
                filtered_indices = category_indices
            
            # Display with pagination using category-specific values
            display_plots_paginated(all_plots, filtered_indices, 
                                  cat_values, f"display_{selected_category}")
        else:
            st.error("Selected category not found")

if __name__ == '__main__':
    main()