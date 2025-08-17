"""
⚠️  CRITICAL WARNING: NEVER USE adata.raw IN THIS PIPELINE ⚠️

Many scanpy functions (sc.tl.rank_genes_groups, sc.tl.score_genes, etc.) 
automatically use .raw when it exists, causing them to analyze raw counts 
instead of normalized data. This creates silent bugs where:

- DE analysis runs on raw counts (wrong)
- Gene scoring runs on raw counts (wrong) 
- Results don't correlate with expected normalized analysis

SOLUTION: Keep original adata object throughout pipeline and transfer 
results back directly. Never set adata.raw = anything.
"""

import os
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
import gc
# Removed gffutils import - gene annotations done in sublibrary_annotation
# Removed PdfPages import - UMAP plotting moved to separate script

# GPU acceleration support
try:
    import rapids_singlecell as rsc
    GPU_AVAILABLE = True
except ImportError:
    GPU_AVAILABLE = False
    rsc = None
# Add parent directory to path for imports
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.pipeline_utils import (
    log_print, setup_gpu, cleanup_gpu_memory,
    normalize_and_preprocess, ensure_cpu_data, perform_dimensionality_reduction, perform_clustering,
    add_comprehensive_gene_annotations, add_mitochondrial_metrics
)


def setup_directories(outdir):
    """Create output directory if it doesn't exist."""
    os.makedirs(outdir, exist_ok=True)
    return outdir


# Removed add_comprehensive_gene_annotations_original - already done in sublibrary_annotation.py


def load_combined_data(input_file):
    """Load combined h5ad file from combine_sublibraries output."""
    log_print(f"📁 Loading combined data from: {input_file}")
    adata = sc.read_h5ad(input_file)
    log_print(f"📊 Combined data shape: {adata.shape}")
    log_print(f"📋 Available layers: {list(adata.layers.keys())}")
    
    # Check that guide data is present
    if 'guide_counts' not in adata.obsm:
        raise RuntimeError("CRITICAL: guide_counts not found in obsm. Expected combined file from combine_sublibraries.")
    
    log_print(f"📊 Guides in obsm: {adata.obsm['guide_counts'].shape[1]} guides")
    
    return adata

    # X should already contain total counts from combine_sublibraries
    log_print("✅ Using X as total counts from combine_sublibraries")
    
    if adata.layers:
        log_print(f"📋 Available layers: {list(adata.layers.keys())}")
    else:
        log_print("📋 No layers present")
    return adata


def load_and_filter_data(input_file):
    """Load combined data and apply comprehensive annotation."""

    # Load combined file from combine_sublibraries
    adata = load_combined_data(input_file)

    # Verify required annotations are present from sublibrary_annotation
    required_var = ['gene_type', 'gene', 'cell_cycle', 'cell_cycle_phase']
    missing_var = [col for col in required_var if col not in adata.var.columns]
    if missing_var:
        raise RuntimeError(f"CRITICAL: Required gene annotations missing: {missing_var}. These should be added by sublibrary_annotation.py")
    
    log_print("✅ Required gene annotations present from sublibrary_annotation")
    log_print(f"📊 Cell cycle genes available: {adata.var['cell_cycle'].sum()}")

    log_print("✨ Using combined format from combine_sublibraries")
    
    # Verify required annotations are present from sublibrary_annotation
    required_obs = ['biological_sample', 'pct_counts_mt', 'guides_per_cell', 'guide_umi_counts']
    missing_obs = [col for col in required_obs if col not in adata.obs.columns]
    if missing_obs:
        log_print(f"⚠️ WARNING: Missing required obs columns: {missing_obs}")
    else:
        log_print(f"✅ All required obs annotations present: {required_obs}")
    
    # Verify gene annotations are present
    if 'gene' not in adata.var.columns:
        log_print("⚠️ WARNING: Gene symbols missing - some analysis may not work")
    
    # Cell filtering now happens at sublibrary combination level
    # No additional filtering applied here
    gc.collect()

    return adata


# The old guide processing and scoring functions have been replaced by the unified scoring system

def _score_gene_list(adata_main, gene_list, score_name, use_gpu=False):
    """Score a gene list on normalized, log-transformed counts (no scaling needed).
    
    Note: Scaling is NOT needed for gene scoring as recommended in:
    https://satijalab.org/seurat/articles/cell_cycle_vignette
    Gene scoring should be performed on normalized, log-transformed counts.
    """
    if len(gene_list) == 0:
        log_print(f"⚠️ No genes found for {score_name}")
        return

    # Use all genes as universe (no subset needed)
    log_print(f"📊 Scoring {len(gene_list)} genes for {score_name} on all {adata_main.n_vars} genes")

    # Score genes directly on normalized, log-transformed data (no scaling)
    log_print(f"🔍 DEBUG: adata_main.X.max() RIGHT before score_genes call: {adata_main.X.max()}")
    if use_gpu:
        rsc.tl.score_genes(adata_main, gene_list=gene_list, score_name=score_name)
    else:
        sc.tl.score_genes(adata_main, gene_list=gene_list, score_name=score_name)

    log_print(f"✅ {score_name}: scored {len(gene_list)} genes")

# All deprecated scoring functions removed - now using unified scoring system


# ============================================================================
# SIMPLIFIED UNIFIED SCORING SYSTEM
# ============================================================================

def apply_subset_filters(adata, subset_configs, min_cells=0):
    """Apply subset filters to create cell masks.
    
    Unified function used by both DE analysis and UMAP subset generation.
    Each config contains a 'name' and 'filter' expression.
    
    Args:
        adata: AnnData object
        subset_configs: List of subset configs with 'name' and 'filter'
        min_cells: Minimum cells required for a subset (default 0)
    
    Returns dict of subset_name -> cell_mask (boolean array)
    """
    import numpy as np
    import pandas as pd
    
    subset_masks = {}
    
    # Always include all_cells_scores as a subset (represents no DE filtering)
    subset_masks["all_cells_scores"] = np.ones(len(adata), dtype=bool)
    log_print(f"  Subset 'all_cells_scores': {len(adata)} cells (no DE filtering)")
    
    if not subset_configs:
        # No additional subsetting - only all cells
        return subset_masks
    
    for config in subset_configs:
        subset_name = config['name']
        subset_filter = config['filter']
        
        # Apply filter to get cell mask
        eval_context = {
            "adata": adata,
            "np": np,
            "pd": pd
        }
        mask = eval(subset_filter, eval_context)
        
        # Ensure mask is boolean array
        mask = np.asarray(mask, dtype=bool)
        n_cells = mask.sum()
        
        if n_cells >= min_cells:
            subset_masks[subset_name] = mask
            log_print(f"  Subset '{subset_name}': {n_cells} cells")
        else:
            log_print(f"  Subset '{subset_name}': {n_cells} cells (below minimum {min_cells}, skipping)")
    
    return subset_masks


def target_gene_strategy(adata, group_def, cutoff, non_targeting_strings):
    """Strategy for target gene groupings: target_genes vs non_targeting"""
    # Get guide names and counts from obsm
    guide_names = adata.uns["guide_names"]
    guide_counts = adata.obsm["guide_counts"]

    if hasattr(guide_counts, 'toarray'):
        guide_counts = guide_counts.toarray()

    guide_names_clean = [x.replace("+", ".").replace("-", ".") for x in guide_names]
    guide_df = pd.DataFrame(guide_counts, columns=guide_names_clean, index=adata.obs.index)

    # Find non-targeting controls
    nt_cols = [x for x in guide_names_clean if any(nt_str in x for nt_str in non_targeting_strings)]
    NT = (guide_df[nt_cols] >= cutoff).sum(axis=1) >= 1 if nt_cols else pd.Series([False] * len(adata), index=adata.obs.index)

    # Handle single gene or gene combination
    if isinstance(group_def, str):
        gene = group_def
        gene_cols = [x for x in guide_names_clean if gene in x and not any(nt_str in x for nt_str in non_targeting_strings)]
        target_mask = (guide_df[gene_cols] >= cutoff).sum(axis=1) >= 1 if gene_cols else pd.Series([False] * len(adata), index=adata.obs.index)
        group_name = f"{gene}_category"
    else:
        genes = group_def
        target_masks = []
        for gene in genes:
            gene_cols = [x for x in guide_names_clean if gene in x and not any(nt_str in x for nt_str in non_targeting_strings)]
            mask = (guide_df[gene_cols] >= cutoff).sum(axis=1) >= 1 if gene_cols else pd.Series([False] * len(adata), index=adata.obs.index)
            target_masks.append(mask)
        target_mask = np.logical_or.reduce(target_masks) if target_masks else np.array([False] * len(adata.obs))
        group_name = f"{'_'.join(genes)}_category"

    # Create categories
    categories = pd.Series("other", index=adata.obs.index)
    categories.loc[NT] = "non_targeting"
    categories.loc[target_mask] = "target_genes"
    categories = categories.astype('category')

    log_print(f"  Created {group_name}: {categories.value_counts().to_dict()}")

    return {
        "group_name": group_name,
        "categories": categories,
        "test_group": "target_genes",
        "reference_group": "non_targeting",
        "score_prefix": group_name.replace('_category', '')
    }


def filter_group_strategy(adata, group_def):
    """Filter-based grouping strategy for differential expression.
    
    group_def should be a list of dicts with 'name' and 'filter'.
    First group is reference, second group is test.
    """
    if len(group_def) != 2:
        raise ValueError(f"filter_group_strategy requires exactly 2 groups, got {len(group_def)}")
    
    # First group is reference, second is test
    reference_group = group_def[0]['name']
    test_group = group_def[1]['name']
    
    # Create categories series
    categories = pd.Series('other', index=adata.obs.index)
    
    # Apply each filter
    for group in group_def:
        group_name = group['name']
        filter_expr = group['filter']
        
        # Evaluate filter expression
        mask = eval(filter_expr)
        categories.loc[mask] = group_name
        log_print(f"    Group '{group_name}': {mask.sum()} cells")
    
    categories = categories.astype('category')
    
    # Generate group name based on test vs reference
    group_name = f"{test_group}_vs_{reference_group}"
    score_prefix = group_name
    
    log_print(f"  Created {group_name}: {categories.value_counts().to_dict()}")
    
    return {
        "group_name": group_name,
        "categories": categories,
        "test_group": test_group,
        "reference_group": reference_group,
        "score_prefix": score_prefix
    }


def cluster_strategy(adata, group_def):
    """Strategy for cluster groupings: specific cluster vs rest"""
    # Find clustering result closest to target number of clusters
    target_clusters = group_def
    best_col = None
    best_n_clusters = None
    min_diff = float('inf')

    for col in adata.obs.columns:
        if col.startswith('clusters') and '_res' in col:
            # Extract number of clusters from column name (e.g., clusters3_res0.15 -> 3)
            n_clusters = int(col.split('_')[0].replace('clusters', ''))
            diff = abs(n_clusters - target_clusters)
            
            # Update if this is closer to target, or same distance but fewer clusters (prefer simpler)
            if diff < min_diff or (diff == min_diff and n_clusters < best_n_clusters):
                min_diff = diff
                best_col = col
                best_n_clusters = n_clusters
                
            # If exact match found, use it immediately
            if diff == 0:
                break

    if best_col is None:
        raise ValueError(f"No clustering results found in adata.obs")

    categories = adata.obs[best_col].astype('category')

    if best_n_clusters != target_clusters:
        log_print(f"  ⚠️  No exact match for {target_clusters} clusters - using closest: {best_col} with {best_n_clusters} clusters")
    else:
        log_print(f"  Using {best_col}: {categories.value_counts().to_dict()}")

    # Column name already has the format we want
    return {
        "group_name": best_col,
        "categories": categories,
        "test_group": None,  # Will be set per cluster
        "reference_group": "rest",
        "score_prefix": best_col
    }




def _perform_single_de_comparison(adata, adata_subset, group_name, test_group, reference_group, score_prefix, split_label, use_gpu):
    """Helper function to perform a single DE comparison and scoring"""
    # Create comparison description for logging
    comparison_desc = f"{test_group} vs {reference_group}"
    log_print(f"    DE analysis: {comparison_desc}...")

    # Create CPU copy for DE
    adata_subset_cpu = ensure_cpu_data(adata_subset, f"{comparison_desc} DE analysis")

    # Perform DE analysis
    sc.tl.rank_genes_groups(adata_subset_cpu, group_name, groups=[test_group], reference=reference_group, method='t-test', use_raw=False)
    de_genes = sc.get.rank_genes_groups_df(adata_subset_cpu, group=test_group, pval_cutoff=None)  # Get all results
    
    # Store DE results in main adata object
    de_key = f"de_{score_prefix}_{test_group}_vs_{reference_group}_{split_label}" if split_label != "all_cells_scores" else f"de_{score_prefix}_{test_group}_vs_{reference_group}"
    if 'de_results' not in adata.uns:
        adata.uns['de_results'] = {}
    adata.uns['de_results'][de_key] = de_genes

    # DEBUG: Print the actual DE dataframe (significant only)
    de_genes_sig = de_genes[de_genes['pvals_adj'] < 0.05] if 'pvals_adj' in de_genes.columns else de_genes[de_genes['pvals'] < 0.05]
    log_print(f"    🔍 DEBUG: DE results for {comparison_desc} (p<0.05, top 15 and bottom 15 rows):")
    if len(de_genes_sig) > 30:
        log_print(f"\n{de_genes_sig.head(15).to_string()}")
        log_print("\n... [middle rows omitted] ...\n")
        log_print(f"\n{de_genes_sig.tail(15).to_string()}")
    else:
        log_print(f"\n{de_genes_sig.to_string()}")

    # DEBUG: Verify negative gene approach by comparing with swapped groups
    if reference_group != "rest":  # Only for pairwise comparisons, not cluster vs rest
        log_print(f"    🔍 DEBUG: Verifying negative gene approach by swapping groups...")
        
        # Perform swapped DE analysis
        sc.tl.rank_genes_groups(adata_subset_cpu, group_name, groups=[reference_group], reference=test_group, method='t-test', use_raw=False)
        de_genes_swapped = sc.get.rank_genes_groups_df(adata_subset_cpu, group=reference_group, pval_cutoff=None)
        de_genes_swapped_sig = de_genes_swapped[de_genes_swapped['pvals_adj'] < 0.05] if 'pvals_adj' in de_genes_swapped.columns else de_genes_swapped[de_genes_swapped['pvals'] < 0.05]
        
        # Get negative genes from original approach
        neg_genes_original = de_genes_sig[de_genes_sig["logfoldchanges"] <= -0.25].sort_values("logfoldchanges", ascending=True).iloc[:100]
        
        # Get positive genes from swapped approach
        pos_genes_swapped = de_genes_swapped_sig[de_genes_swapped_sig["logfoldchanges"] >= 0.25].sort_values("logfoldchanges", ascending=False).iloc[:100]
        
        # Compare gene lists
        orig_gene_set = set(neg_genes_original["names"].tolist())
        swap_gene_set = set(pos_genes_swapped["names"].tolist())
        
        overlap = len(orig_gene_set & swap_gene_set)
        orig_only = len(orig_gene_set - swap_gene_set)
        swap_only = len(swap_gene_set - orig_gene_set)
        
        log_print(f"    📊 Negative gene approach verification:")
        log_print(f"      Original approach (negative LFC): {len(orig_gene_set)} genes")
        log_print(f"      Swapped approach (positive LFC): {len(swap_gene_set)} genes")
        log_print(f"      Overlap: {overlap} genes")
        log_print(f"      Original only: {orig_only} genes")
        log_print(f"      Swapped only: {swap_only} genes")
        
        if overlap > 0:
            # Show LFC correlation for overlapping genes
            overlap_genes = list(orig_gene_set & swap_gene_set)[:10]  # First 10 overlapping genes
            if overlap_genes:
                log_print(f"    📊 LFC comparison for overlapping genes (first 10):")
                for gene in overlap_genes:
                    orig_lfc = neg_genes_original[neg_genes_original["names"] == gene]["logfoldchanges"].iloc[0]
                    swap_lfc = pos_genes_swapped[pos_genes_swapped["names"] == gene]["logfoldchanges"].iloc[0]
                    log_print(f"      {gene}: original={orig_lfc:.3f}, swapped={swap_lfc:.3f}, sum={orig_lfc + swap_lfc:.6f}")
        
        # Restore original DE results for downstream processing
        sc.tl.rank_genes_groups(adata_subset_cpu, group_name, groups=[test_group], reference=reference_group, method='t-test', use_raw=False)

    # Create score names based on comparison type
    suffix = f"_{split_label}" if split_label != "all_cells_scores" else ""
    if reference_group == "rest":
        # Cluster vs rest naming: include cluster name with 'clust' prefix
        pos_score_name = f"{score_prefix}_clust{test_group}_pos{suffix}"
        neg_score_name = f"{score_prefix}_clust{test_group}_neg{suffix}"
    else:
        # Pairwise comparison naming: use generic "score" 
        pos_score_name = f"{score_prefix}_score_pos{suffix}"
        neg_score_name = f"{score_prefix}_score_neg{suffix}"

    # Score genes (use minimal fold change cutoffs, sort by score for best ranking)
    created_scores = []
    pos_genes = de_genes_sig[de_genes_sig["logfoldchanges"] >= 0.25].sort_values("scores", ascending=False).iloc[:100]
    if len(pos_genes) > 0:
        log_print(f"    📊 Using top {len(pos_genes)} upregulated genes for {pos_score_name} (sorted by score)")
        log_print(f"    🔍 DEBUG: Positive genes for scoring ({pos_score_name}):")
        log_print(f"\n{pos_genes.to_string()}")
        _score_gene_list(adata, pos_genes["names"].tolist(), pos_score_name, use_gpu)
        created_scores.append(pos_score_name)
        log_print(f"    ✅ {pos_score_name}: {len(pos_genes)} genes")

    neg_genes = de_genes_sig[de_genes_sig["logfoldchanges"] <= -0.25].sort_values("scores", ascending=True).iloc[:100]
    if len(neg_genes) > 0:
        log_print(f"    📊 Using top {len(neg_genes)} downregulated genes for {neg_score_name} (sorted by score)")
        log_print(f"    🔍 DEBUG: Negative genes for scoring ({neg_score_name}):")
        log_print(f"\n{neg_genes.to_string()}")
        _score_gene_list(adata, neg_genes["names"].tolist(), neg_score_name, use_gpu)
        created_scores.append(neg_score_name)
        log_print(f"    ✅ {neg_score_name}: {len(neg_genes)} genes")

    return created_scores


def perform_de_and_scoring(adata, group_config, cell_mask, split_label, use_gpu):
    """Clean DE and scoring function - works for all strategies"""
    group_name = group_config["group_name"]
    test_group = group_config["test_group"]
    reference_group = group_config["reference_group"]
    score_prefix = group_config["score_prefix"]

    # Subset to relevant cells using mask (create copy to avoid view warnings)
    adata_subset = adata[cell_mask].copy()

    # Determine comparison pairs
    if reference_group == "rest" and test_group is None:
        # Multi-cluster case: generate list of (cluster, "rest") pairs
        unique_clusters = adata_subset.obs[group_name].unique()
        comparisons = [(cluster, "rest") for cluster in unique_clusters]
    else:
        # Single pairwise case: validate groups exist
        value_counts = adata_subset.obs[group_name].value_counts()
        if test_group not in value_counts.index or reference_group not in value_counts.index:
            log_print(f"    ⚠️  Skipping {group_name}: missing required groups")
            return []
        comparisons = [(test_group, reference_group)]

    # Execute all comparisons using the helper function
    all_created_scores = []
    for test, ref in comparisons:
        scores = _perform_single_de_comparison(
            adata, adata_subset, group_name, test, ref, 
            score_prefix, split_label, use_gpu
        )
        all_created_scores.extend(scores)
        log_print(f"    📋 DE comparison {test} vs {ref} created {len(scores)} scores: {scores}")

    log_print(f"    📋 Total scores created for {score_prefix}: {len(all_created_scores)} scores: {all_created_scores}")
    return all_created_scores


def calculate_unified_scores(adata, scoring_config, use_gpu, cutoff, non_targeting_strings, de_analysis_subsets=None):
    """Simplified unified scoring system with strategy functions"""
    log_print("🚀 Starting simplified unified scoring system...")

    # Strategy function mapping
    strategy_functions = {
        "target_gene_strategy": target_gene_strategy,
        "cluster_strategy": cluster_strategy,
        "filter_group_strategy": filter_group_strategy
    }

    created_categories = []
    created_scores = []
    
    # Track DE subset structure for plotting
    de_subset_structure = {}
    de_subsets_found = set()

    # Process each configuration
    for config_name, config in scoring_config.items():
        strategy = config["strategy"]
        strategy_func = strategy_functions[strategy]

        log_print(f"🎯 Processing {config_name} using {strategy}")

        for group_def in config["groups"]:
            # Create grouping using strategy function
            if strategy == "target_gene_strategy":
                group_config = strategy_func(adata, group_def, cutoff=cutoff, non_targeting_strings=non_targeting_strings)
            else:
                group_config = strategy_func(adata, group_def)

            # Skip if strategy function returns None (e.g., missing cluster count)
            if group_config is None:
                continue

            # Add categories to adata
            adata.obs[group_config["group_name"]] = group_config["categories"]
            created_categories.append(group_config["group_name"])

            # Apply DE splitting - returns cell masks for separate DE analyses
            de_split_masks = apply_subset_filters(
                adata,
                de_analysis_subsets  # Use the parameter passed to main(), not from config
            )

            # Perform DE and scoring for each DE split
            for split_label, cell_mask in de_split_masks.items():
                if cell_mask.sum() == 0:  # No cells in this DE split
                    continue

                scores = perform_de_and_scoring(adata, group_config, cell_mask, split_label, use_gpu)
                log_print(f"    📋 perform_de_and_scoring returned {len(scores)} scores: {scores}")
                if scores:
                    created_scores.extend(scores)
                    
                    # Track DE subset structure (only real subsets, not all_cells_scores)
                    if split_label != "all_cells_scores":
                        de_subsets_found.add(split_label)
                        if split_label not in de_subset_structure:
                            de_subset_structure[split_label] = []
                        de_subset_structure[split_label].extend(scores)

    # Store for automatic UMAP registration
    adata.uns["unified_scoring_categories"] = created_categories
    adata.uns["unified_scoring_scores"] = created_scores
    
    # Store DE subset structure for plotting
    if de_subsets_found:
        adata.uns["de_subsets"] = sorted(list(de_subsets_found))
        adata.uns["de_subset_scores"] = de_subset_structure
        log_print(f"  📊 DE subsets found: {adata.uns['de_subsets']}")
        for subset, scores in de_subset_structure.items():
            log_print(f"    - {subset}: {len(scores)} scores")

    log_print(f"🎆 Simplified unified scoring complete!")
    log_print(f"  📊 Created {len(created_categories)} categories: {created_categories}")
    log_print(f"  📊 Created {len(created_scores)} scores: {created_scores}")

    return adata


def calculate_subset_umaps(adata_hvg, subset_configs, use_gpu, n_neighbors, n_pcs):
    """Calculate UMAP embeddings for cell subsets using existing PCA coordinates.
    
    This function computes separate UMAP embeddings for defined cell subsets,
    allowing focused visualization of specific cell populations. It reuses the
    PCA coordinates from the full dataset but recalculates neighbors and UMAP
    for each subset.
    
    Args:
        adata_hvg: AnnData object with PCA coordinates in obsm['X_pca']
        subset_configs: List of subset configurations with 'name' and 'filter'
        use_gpu: Whether to use GPU acceleration
    """
    import numpy as np
    
    if not subset_configs:
        log_print("ℹ️  No subset configurations provided for UMAP calculation")
        return
    
    log_print(f"📊 Calculating UMAP embeddings for {len(subset_configs)} cell subsets...")
    
    # Import GPU libraries if needed
    if use_gpu and GPU_AVAILABLE:
        import rapids_singlecell as rsc
    
    # Use unified subset filter function with minimum cell requirement
    subset_masks = apply_subset_filters(adata_hvg, subset_configs, min_cells=100)
    
    successful_subsets = []
    
    for subset_name, mask in subset_masks.items():
        if subset_name == "all_cells_scores":
            continue  # Skip all_cells_scores for UMAP subsets (it's for DE, not UMAP)
            
        n_cells = mask.sum()
        
        # Create subset
        subset_adata = adata_hvg[mask].copy()
        log_print(f"  🗺️  Computing UMAP for subset '{subset_name}': {n_cells:,} cells")
        
        # Recalculate neighbors and UMAP using existing PCA
        # The PCA coordinates are already in obsm['X_pca']
        if use_gpu and GPU_AVAILABLE:
            rsc.pp.neighbors(subset_adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
            rsc.tl.umap(subset_adata)
        else:
            sc.pp.neighbors(subset_adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
            sc.tl.umap(subset_adata)
        
        # Store subset UMAP coordinates in main object
        # Use NaN for cells not in subset to maintain array shape
        umap_key = f'X_umap_{subset_name}'
        adata_hvg.obsm[umap_key] = np.full((adata_hvg.n_obs, 2), np.nan)
        adata_hvg.obsm[umap_key][mask] = subset_adata.obsm['X_umap']
        
        successful_subsets.append(subset_name)
        log_print(f"    ✅ Successfully computed UMAP for '{subset_name}'")
    
    if successful_subsets:
        log_print(f"✅ Successfully calculated {len(successful_subsets)} subset UMAPs: {successful_subsets}")
        
        # Store list of subset names in uns for downstream use
        adata_hvg.uns['umap_subsets'] = successful_subsets
    else:
        log_print("⚠️  No subset UMAPs were successfully calculated")


# UMAP plotting function moved to scripts/plot_standard_analyses_umap.py
# This separation allows for:
# - Independent plot regeneration without re-running analysis
# - Dashboard-compatible individual PNG files instead of multi-page PDF
# - Better memory management and parallelization


def prepare_final_dataset(adata, final_output_file):
    """Prepare final dataset for downstream analysis."""

    log_print("🧬 Preparing final dataset")

    # Use the current data (already has total counts in X and all analysis)
    adata_final = adata.copy()

    # X contains raw counts from combine_sublibraries
    log_print("✅ Using raw counts in X from combine_sublibraries")

    # Add metadata describing the data structure
    adata_final.uns["data_description"] = {
        "X_contains": "raw total RNA counts",
        "normalization": "Analysis performed on normalized subsets (target_sum=1e4 + log1p), but X contains raw counts",
        "data_type": "single-cell perturbation experiment",
        "created_by": "standard_analyses.py prepare_final_dataset()"
    }

    # Remove .raw if it exists
    adata_final.raw = None

    # Clean up obs columns (remove guide data now stored in obsm)
    guide_cols = [col for col in adata_final.obs.columns if "_guide" in col]
    # Get target gene columns from the uns metadata
    target_gene_names = adata_final.uns.get("guide_target_names", [])
    target_cols = [col for col in target_gene_names if col in adata_final.obs.columns]
    cols_to_remove = guide_cols + target_cols
    if cols_to_remove:
        log_print(f"🧽 Removing {len(cols_to_remove)} guide columns from obs (now in obsm): {cols_to_remove[:5]}...")
        adata_final.obs = adata_final.obs.drop(columns=cols_to_remove)

    log_print(f"🎆 Created final dataset with shape: {adata_final.shape}")
    log_print(f"📋 Layers: {list(adata_final.layers.keys())}")
    if adata_final.X is not None:
        log_print(f"📊 Total counts range: {adata_final.X.sum(axis=1).min():.0f} - {adata_final.X.sum(axis=1).max():.0f}")


    log_print(f"💾 Writing final dataset to {final_output_file}...")
    adata_final.write_h5ad(final_output_file)
    log_print("✅ Final dataset saved successfully!")

    return adata_final


def main(cutoff, input_file, outdir, final_output_file, 
         use_gpu, target_clusters,
         target_genes=None, de_analysis_subsets=None, 
         non_targeting_strings=None,
         gpu_threshold_gb=2.0, hvg_reference_file=None,
         standard_analyses_config=None):
    """Main processing pipeline for unified data format.

    Parameters:
    -----------
    cutoff : int, default=5
        Threshold for guide assignment
    target_genes : list, optional
        List of target genes to analyze. If None, no target gene analysis is performed
    input_file : str
        Path to input h5ad file with annotated data from combine_sublibraries
    outdir : str
        Output directory for plots and results
    final_output_file : str
        Path for final h5ad output ready for training
    use_gpu : bool, default=False
        Whether to use GPU acceleration for computations
    gpu_threshold_gb : float, default=2.0
        Maximum used GPU memory allowed (in GB) - GPUs with less usage will be selected
    """
    # Setup output paths and GPU if requested
    outdir = setup_directories(outdir)
    gpu_enabled = setup_gpu(use_gpu, gpu_threshold_gb=gpu_threshold_gb)

    log_print(f"🚀 Processing mode: {'GPU' if gpu_enabled else 'CPU'}")
    log_print(f"📁 Loading combined data from: {input_file}")
    log_print(f"📁 Output directory: {outdir}")
    log_print(f"📁 Final output: {final_output_file}")

    try:
        # Load and preprocess data
        adata = load_and_filter_data(input_file)
        log_print(f"📊 Available layers: {list(adata.layers.keys())}")

        # Gene annotations already added in load_and_filter_data

        # Library labels should already be correct from combine_sublibraries
        log_print("📊 Using existing library labels from combine_sublibraries")

        # Normalize and create HVG subset for memory-efficient processing
        adata_hvg, adata_expressed = normalize_and_preprocess(adata, use_gpu=gpu_enabled)
        log_print(f"📊 Working with HVG subset: {adata_hvg.shape}")
        log_print(f"📊 Working with expressed subset: {adata_expressed.shape}")

        # Clean GPU memory after preprocessing
        if gpu_enabled:
            cleanup_gpu_memory()

        # Perform clustering on HVG dataset FIRST (before unified scoring)
        log_print("🎯 Performing clustering on HVG dataset...")
        n_neighbors = standard_analyses_config['n_neighbors']
        n_pcs = standard_analyses_config['n_pcs']
        adata_hvg = perform_dimensionality_reduction(adata_hvg, use_gpu=gpu_enabled, n_neighbors=n_neighbors, n_pcs=n_pcs)

        # Clean GPU memory after dimensionality reduction
        if gpu_enabled:
            cleanup_gpu_memory()

        adata_hvg = perform_clustering(adata_hvg, use_gpu=gpu_enabled, target_clusters=target_clusters)

        # Clean GPU memory after clustering
        if gpu_enabled:
            cleanup_gpu_memory()

        # Calculate subset UMAPs if configured
        # Use the same analysis_subsets for UMAP
        if de_analysis_subsets:
            calculate_subset_umaps(adata_hvg, de_analysis_subsets, use_gpu=gpu_enabled, n_neighbors=n_neighbors, n_pcs=n_pcs)

        # Use expressed dataset for DE and scoring (already available from normalize_and_preprocess)

        # Transfer clustering results from HVG to expressed dataset
        log_print("🔄 Transferring clustering results to expressed dataset...")
        cluster_cols = [col for col in adata_hvg.obs.columns if col.startswith('clusters') and '_res' in col]
        for col in cluster_cols:
            adata_expressed.obs[col] = adata_hvg.obs[col]
        log_print(f"📊 Transferred {len(cluster_cols)} clustering results: {cluster_cols}")

        # Delete layers from expressed dataset (not needed for DE/scoring)
        if adata_expressed.layers:
            for layer in list(adata_expressed.layers.keys()):
                del adata_expressed.layers[layer]
            log_print("🧹 Deleted layers from expressed genes dataset")

        # Add cell cycle scores using standard scoring function
        log_print("🧬 Adding cell cycle scores...")
        
        # Use cell cycle gene annotations from sublibrary_annotation
        # Get S and G2M genes from the annotations
        s_genes_mask = (adata_expressed.var['cell_cycle'] == True) & (adata_expressed.var['cell_cycle_phase'] == 'S')
        g2m_genes_mask = (adata_expressed.var['cell_cycle'] == True) & (adata_expressed.var['cell_cycle_phase'] == 'G2M')
        
        s_genes_ensembl = adata_expressed.var_names[s_genes_mask].tolist()
        g2m_genes_ensembl = adata_expressed.var_names[g2m_genes_mask].tolist()
        
        log_print(f"📊 Found {len(s_genes_ensembl)} S phase genes and {len(g2m_genes_ensembl)} G2M genes from annotations")
        
        # Score using standard scoring function (same as perturbation genes)
        _score_gene_list(adata_expressed, s_genes_ensembl, 'S_score', use_gpu=gpu_enabled)
        _score_gene_list(adata_expressed, g2m_genes_ensembl, 'G2M_score', use_gpu=gpu_enabled)
        
        log_print(f"📊 Cell cycle scores added: S_score ({len(s_genes_ensembl)} genes), G2M_score ({len(g2m_genes_ensembl)} genes)")

        # Clean GPU memory after cell cycle scoring
        if gpu_enabled:
            cleanup_gpu_memory()

        # Verify guide counts are already calculated from sublibrary_annotation
        if 'guides_per_cell' not in adata_expressed.obs:
            raise RuntimeError("CRITICAL: guides_per_cell not found. Should be calculated in sublibrary_annotation.py")
        
        # Use existing guide counts (rename for consistency)
        adata_expressed.obs['n_guides_total'] = adata_expressed.obs['guides_per_cell']
        log_print(f"📊 Using existing guide counts: mean = {adata_expressed.obs['n_guides_total'].mean():.2f}")

        # Configure unified scoring system
        scoring_config = {}
        
        # Only add target_genes config if target genes are specified
        if target_genes:
            # Create list of individual genes and optionally gene pairs
            groups = list(target_genes)  # Individual genes
            # If there are at least 2 genes, add them as a pair for combined analysis
            if len(target_genes) >= 2:
                groups.append(target_genes[:2])  # Add first two genes as a pair
            
            scoring_config["target_genes"] = {
                "groups": groups,
                "strategy": "target_gene_strategy",
                "de_analysis_subsets": de_analysis_subsets if de_analysis_subsets else []
            }
        
        # Add configured scoring groups from config
        scoring_groups = standard_analyses_config.get('scoring_groups', {})
        for group_name, group_list in scoring_groups.items():
            scoring_config[group_name] = {
                "groups": [group_list],  # Wrap in list for compatibility with existing code
                "strategy": "filter_group_strategy"
            }
            # First is reference, second is test
            log_print(f"✅ Added scoring group '{group_name}': {group_list[1]['name']} vs {group_list[0]['name']}")
        
        # Use target_clusters to determine range of clusters to test
        # Test from 2 up to target_clusters + 1
        scoring_config["leiden_clusters"] = {
            "groups": list(range(2, target_clusters + 2)),
            "strategy": "cluster_strategy"
        }

        # DEBUG: Save data before DE/scoring for debugging
        debug_file = f"{outdir}/debug_before_de_scoring.h5ad"
        log_print(f"💾 DEBUG: Saving data before DE/scoring to {debug_file}")
        #adata_expressed.write_h5ad(debug_file)
        log_print("✅ Debug checkpoint saved - use this file with debug_unified_scoring.py")

        # Check data state before scoring
        log_print(f"🔍 DEBUG: adata_expressed.X.max() before scoring: {adata_expressed.X.max()}")
        log_print(f"🔍 DEBUG: adata_expressed.X data type: {type(adata_expressed.X)}")
        
        # Run unified scoring system
        # Use empty list if non_targeting_strings is None
        nt_strings = non_targeting_strings if non_targeting_strings else []
        adata_expressed = calculate_unified_scores(adata_expressed, scoring_config, use_gpu=gpu_enabled, cutoff=cutoff, non_targeting_strings=nt_strings, de_analysis_subsets=de_analysis_subsets)

        # Clean GPU memory after unified scoring
        if gpu_enabled:
            cleanup_gpu_memory()

        # Transfer analysis results directly from expressed dataset to full dataset
        log_print("🔄 Transferring analysis results from expressed dataset to full dataset...")

        # Transfer ALL obs columns from expressed dataset (only if they don't exist)
        obs_transferred = 0
        obs_skipped = 0
        for col in adata_expressed.obs.columns:
            if col not in adata.obs.columns:
                adata.obs[col] = adata_expressed.obs[col]
                log_print(f"  Transferred obs column: {col}")
                obs_transferred += 1
            else:
                log_print(f"  Skipped obs column: {col} (already exists in full dataset)")
                obs_skipped += 1

        # Transfer ALL obs columns from HVG dataset (clustering, etc.)
        for col in adata_hvg.obs.columns:
            if col not in adata.obs.columns:
                adata.obs[col] = adata_hvg.obs[col]
                log_print(f"  Transferred obs column from HVG: {col}")
                obs_transferred += 1
            else:
                log_print(f"  Skipped obs column from HVG: {col} (already exists in full dataset)")
                obs_skipped += 1

        # Transfer ALL obsm keys from HVG (PCA, UMAP, etc.)
        obsm_transferred = 0
        obsm_skipped = 0
        for key in adata_hvg.obsm.keys():
            if key not in adata.obsm:
                adata.obsm[key] = adata_hvg.obsm[key]
                log_print(f"  Transferred obsm key: {key}")
                obsm_transferred += 1
            else:
                log_print(f"  Skipped obsm key: {key} (already exists in full dataset)")
                obsm_skipped += 1

        # Transfer ALL uns keys from expressed dataset (DE results, scoring metadata, etc.)
        uns_transferred = 0
        uns_skipped = 0
        for key in adata_expressed.uns.keys():
            if key not in adata.uns:
                adata.uns[key] = adata_expressed.uns[key]
                log_print(f"  Transferred uns key: {key}")
                uns_transferred += 1
            else:
                log_print(f"  Skipped uns key: {key} (already exists in full dataset)")
                uns_skipped += 1

        # Transfer ALL uns keys from HVG dataset
        for key in adata_hvg.uns.keys():
            if key not in adata.uns:
                adata.uns[key] = adata_hvg.uns[key]
                log_print(f"  Transferred uns key from HVG: {key}")
                uns_transferred += 1
            else:
                log_print(f"  Skipped uns key from HVG: {key} (already exists in full dataset)")
                uns_skipped += 1

        # Transfer ALL obs columns to HVG for UMAP plotting
        log_print("🔄 Transferring all analysis results to HVG dataset for UMAP plotting...")
        for col in adata.obs.columns:
            if col not in adata_hvg.obs.columns:
                adata_hvg.obs[col] = adata.obs[col]
        
        # Transfer ALL uns keys to HVG for UMAP plotting
        for key in adata.uns.keys():
            if key not in adata_hvg.uns:
                adata_hvg.uns[key] = adata.uns[key]

        # Delete layers from HVG dataset (not needed for UMAP plotting)
        if adata_hvg.layers:
            for layer in list(adata_hvg.layers.keys()):
                del adata_hvg.layers[layer]
            log_print("🧹 Deleted layers from HVG dataset")

        # UMAP plotting moved to separate rule for dashboard integration
        log_print("🎨 UMAP plotting will be done by plot_standard_analyses_umap rule for dashboard integration")

        log_print(f"📊 Transfer summary - obs: {obs_transferred} transferred, {obs_skipped} skipped")
        log_print(f"📊 Transfer summary - obsm: {obsm_transferred} transferred, {obsm_skipped} skipped") 
        log_print(f"📊 Transfer summary - uns: {uns_transferred} transferred, {uns_skipped} skipped")

        log_print(f"✅ All analysis results transferred to full dataset: {adata.shape}")
        adata_full = adata

        # Prepare final dataset using full gene set (this will delete the raw count layers)
        adata_final = prepare_final_dataset(adata_full, final_output_file)

        log_print(f"✅ Processing complete! Final output saved to: {final_output_file}")
        return adata_final

    finally:
        # Clean up GPU memory if it was used
        if gpu_enabled:
            cleanup_gpu_memory()


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Preprocess perturbation data for training')
    parser.add_argument('--config-file', required=True,
                        help='Configuration YAML file')
    parser.add_argument('--input-file', required=True,
                        help='Input combined h5ad file from combine_sublibraries')
    parser.add_argument('--outdir', required=True,
                        help='Output directory for plots and results')
    parser.add_argument('--final-output-file', required=True,
                        help='Path for final h5ad output ready for training')
    args = parser.parse_args()
    
    # Load configuration file
    import yaml
    with open(args.config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    # Extract relevant config sections
    standard_analyses_config = config.get('standard_analyses', {})
    # Use unified analysis_subsets for both DE and UMAP
    analysis_subsets = standard_analyses_config.get('analysis_subsets', [])

    # Note: Snakemake handles dependency tracking, so we always regenerate when called
    if os.path.exists(args.final_output_file):
        log_print(f"🔄 Final output file {args.final_output_file} exists - regenerating as requested by Snakemake...")

    # Run the pipeline - get all parameters from config
    adata_final = main(
        cutoff=standard_analyses_config.get('guide_assignment_cutoff', 5),
        input_file=args.input_file,
        outdir=args.outdir,
        final_output_file=args.final_output_file,
        use_gpu=standard_analyses_config.get('use_gpu', False),
        target_clusters=standard_analyses_config.get('target_clusters', 3),
        target_genes=standard_analyses_config.get('target_genes', []),
        de_analysis_subsets=analysis_subsets,
        non_targeting_strings=standard_analyses_config.get('non_targeting_strings', []),
        gpu_threshold_gb=2.0,
        hvg_reference_file=None,
        standard_analyses_config=standard_analyses_config
    )
