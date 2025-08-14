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
from matplotlib.backends.backend_pdf import PdfPages

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

def apply_stratification(adata, base_samples, stratify_by, stratify_values):
    """Apply stratification strategy to sample selection."""
    if stratify_by is None:
        return {"all": base_samples}

    if stratify_by not in adata.obs.columns:
        raise ValueError(f"Stratification column '{stratify_by}' not found in adata.obs")

    stratified_samples = {}
    for value in stratify_values:
        mask = adata.obs[stratify_by].astype(str).str.contains(str(value), na=False)
        filtered_samples = list(set(adata.obs.loc[mask, "sample"]))
        filtered_samples = [s for s in filtered_samples if s in base_samples]

        if filtered_samples:
            stratified_samples[value] = filtered_samples
            log_print(f"  Stratification '{value}': {len(filtered_samples)} samples")

    return stratified_samples


def target_gene_strategy(adata, group_def, cutoff=5, non_targeting_strings=["non-targeting", "non.targeting"]):
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


def guide_strategy(adata, group_def):
    """Strategy for guide MOI groupings: high vs low guide count"""
    # Use existing guides_per_cell from sublibrary_annotation
    # Use existing guides_per_cell from sublibrary_annotation
    guide_col = 'n_guides_total' if 'n_guides_total' in adata.obs else 'guides_per_cell'
    categories = (adata.obs[guide_col] > 1).map({True: 'high', False: 'low'}).astype('category')

    log_print(f"  Created guide_category: {categories.value_counts().to_dict()}")

    return {
        "group_name": "guide_category",
        "categories": categories,
        "test_group": "high",
        "reference_group": "low",
        "score_prefix": "guide"
    }


def cluster_strategy(adata, group_def):
    """Strategy for leiden cluster groupings: specific cluster vs rest"""
    # Find leiden result with target number of clusters
    target_clusters = group_def
    best_resolution = None

    for col in adata.obs.columns:
        if col.startswith('leiden_'):
            n_clusters = adata.obs[col].nunique()
            if n_clusters == target_clusters:
                best_resolution = col.replace('leiden_', '')
                break

    if best_resolution is None:
        log_print(f"  ⚠️  No leiden clustering found with {target_clusters} clusters - skipping")
        return None  # Signal to skip this grouping

    leiden_col = f"leiden_{best_resolution}"
    categories = adata.obs[leiden_col].astype('category')

    log_print(f"  Using {leiden_col}: {categories.value_counts().to_dict()}")

    return {
        "group_name": f"clusters_{target_clusters}",
        "categories": categories,
        "test_group": None,  # Will be set per cluster
        "reference_group": "rest",
        "score_prefix": f"cluster_{target_clusters}"
    }


def _perform_single_de_comparison(adata, adata_subset, group_name, test_group, reference_group, score_prefix, strat_label, use_gpu):
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
    de_key = f"de_{score_prefix}_{test_group}_vs_{reference_group}_{strat_label}" if strat_label != "all" else f"de_{score_prefix}_{test_group}_vs_{reference_group}"
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
    suffix = f"_{strat_label}" if strat_label != "all" else ""
    if reference_group == "rest":
        # Cluster vs rest naming: include cluster name
        pos_score_name = f"{score_prefix}_{test_group}_pos{suffix}"
        neg_score_name = f"{score_prefix}_{test_group}_neg{suffix}"
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


def perform_de_and_scoring(adata, group_config, strat_samples, strat_label, use_gpu):
    """Clean DE and scoring function - works for all strategies"""
    group_name = group_config["group_name"]
    test_group = group_config["test_group"]
    reference_group = group_config["reference_group"]
    score_prefix = group_config["score_prefix"]

    # Subset to relevant samples (create copy to avoid view warnings)
    adata_subset = adata[adata.obs["sample"].isin(strat_samples)].copy()

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
            score_prefix, strat_label, use_gpu
        )
        all_created_scores.extend(scores)
        log_print(f"    📋 DE comparison {test} vs {ref} created {len(scores)} scores: {scores}")

    log_print(f"    📋 Total scores created for {score_prefix}: {len(all_created_scores)} scores: {all_created_scores}")
    return all_created_scores


def calculate_unified_scores(adata, scoring_config, use_gpu=False, cutoff=5, non_targeting_strings=["non-targeting", "non.targeting"]):
    """Simplified unified scoring system with strategy functions"""
    log_print("🚀 Starting simplified unified scoring system...")

    # Strategy function mapping
    strategy_functions = {
        "target_gene_strategy": target_gene_strategy,
        "guide_strategy": guide_strategy,
        "cluster_strategy": cluster_strategy
    }

    created_categories = []
    created_scores = []
    available_samples = list(adata.obs["sample"].unique())

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

            # Apply stratification
            stratified_samples = apply_stratification(
                adata, available_samples,
                config.get("stratify_by"),
                config.get("stratify_values", [])
            )

            # Perform DE and scoring for each stratification
            for strat_label, samples in stratified_samples.items():
                if not samples:
                    continue

                scores = perform_de_and_scoring(adata, group_config, samples, strat_label, use_gpu)
                log_print(f"    📋 perform_de_and_scoring returned {len(scores)} scores: {scores}")
                if scores:
                    created_scores.extend(scores)

    # Store for automatic UMAP registration
    adata.uns["unified_scoring_categories"] = created_categories
    adata.uns["unified_scoring_scores"] = created_scores

    log_print(f"🎆 Simplified unified scoring complete!")
    log_print(f"  📊 Created {len(created_categories)} categories: {created_categories}")
    log_print(f"  📊 Created {len(created_scores)} scores: {created_scores}")

    return adata


def save_umap_plots(adata, outdir):
    """Save UMAP plots for visualization using clean unified system."""
    umap_colors_candidate = [
        # Basic metrics
        "pct_counts_mt", "pct_counts_ribo", "rep", "library", "lib_for_batch",
        "guide_umi_counts", "total_counts", "guides_per_cell", "n_guides_total",

        # Gene expression scores
        "G2M_score", "S_score"
    ]

    # AUTOMATIC REGISTRATION: Add all categories and scores created by unified scoring system
    if "unified_scoring_categories" in adata.uns:
        unified_categories = adata.uns["unified_scoring_categories"]
        umap_colors_candidate.extend(unified_categories)
        log_print(f"📊 Auto-registered {len(unified_categories)} unified categories: {unified_categories}")

    if "unified_scoring_scores" in adata.uns:
        unified_scores = adata.uns["unified_scoring_scores"]
        umap_colors_candidate.extend(unified_scores)
        log_print(f"📊 Auto-registered {len(unified_scores)} unified scores: {unified_scores}")

    # Add leiden clustering results
    leiden_cols = [col for col in adata.obs.columns if col.startswith('leiden_')]
    umap_colors_candidate.extend(leiden_cols)

    # Filter to only existing columns
    umap_colors = [col for col in umap_colors_candidate if col in adata.obs.columns]

    missing_cols = [col for col in umap_colors_candidate if col not in adata.obs.columns]

    if missing_cols:
        log_print(f"⚠️  WARNING: Missing columns for UMAP plotting: {missing_cols}")
    log_print(f"🗺️ Plotting UMAP with available columns: {umap_colors}")

    # Generate multi-page PDF with each UMAP on its own page
    output_file = f"{outdir}/umap.final_analysis.pdf"

    # Set scanpy settings for rasterized output (not vector)
    sc.settings.set_figure_params(dpi_save=150, facecolor='white')

    with PdfPages(output_file) as pdf:
        for i, color in enumerate(umap_colors):
            log_print(f"📊 Generating UMAP {i+1}/{len(umap_colors)}: {color}")

            # Set figure size before scanpy creates the plot
            plt.figure(figsize=(12, 8))
            
            # Check if this is a target gene category column
            is_target_gene_column = (color in adata.obs.columns and 
                                   adata.obs[color].dtype.name == 'category' and
                                   'target_genes' in adata.obs[color].cat.categories)
            
            if is_target_gene_column:
                # For target gene columns, create a binary version for cleaner visualization
                # Create temporary binary column: target_genes vs everything else
                binary_col = f"{color}_binary"
                adata.obs[binary_col] = adata.obs[color].map(
                    lambda x: 'target_genes' if x == 'target_genes' else 'background'
                ).astype('category')
                
                # Plot using scanpy with custom colors
                sc.pl.umap(adata, color=binary_col, 
                          palette={'background': 'lightgray', 'target_genes': 'red'},
                          size=40, alpha=0.8, show=False, frameon=False)
                
                # Clean up temporary column
                adata.obs.drop(columns=[binary_col], inplace=True)
                
                # Add title
                plt.title(f'{color} (target_genes highlighted)')
                    
                target_count = (adata.obs[color] == 'target_genes').sum()
                log_print(f"  Target_genes only plot: {target_count} cells highlighted on gray background")
            else:
                # Normal plotting for other columns
                sc.pl.umap(adata, color=color, vmin="p5", vmax="p95",
                          show=False, frameon=False)

            # Get the current figure that scanpy created
            fig = plt.gcf()

            # Ensure points are rasterized (not vector) for smaller file size
            for ax in fig.get_axes():
                for child in ax.get_children():
                    if hasattr(child, 'set_rasterized'):
                        child.set_rasterized(True)

            # Save page to PDF
            pdf.savefig(fig, bbox_inches='tight', dpi=150)
            plt.close(fig)  # Close figure to free memory

    log_print(f"📄 Multi-page UMAP saved to: {output_file}")
    log_print(f"📊 Generated {len(umap_colors)} UMAP pages")



def prepare_final_dataset(adata, final_output_file):
    """Prepare final dataset using unified data with combined layers."""

    # Check if we have the combined layers (created earlier in pipeline)
    required_combined_layers = ["total_counts", "nascent_counts"]
    missing_combined_layers = [layer for layer in required_combined_layers if layer not in adata.layers]

    if missing_combined_layers:
        available_layers = list(adata.layers.keys())
        raise RuntimeError(f"CRITICAL: Required combined layers missing: {missing_combined_layers}. "
                          f"Available layers: {available_layers}. "
                          f"Combined layers should have been created earlier in load_and_combine_data().")

    if len(missing_combined_layers) == 0:
        log_print("🧬 Using combined layers created earlier in pipeline")

        # Use the current data (already has total counts in X and all analysis)
        adata_final = adata.copy()

        # Combined layers already exist from load_and_combine_data()
        # X contains raw counts, layers contain the same raw counts
        log_print("✅ Combined layers already present: total_counts and nascent_counts")

        # Add metadata describing the data structure
        adata_final.uns["data_description"] = {
            "X_contains": "raw total RNA counts (mature + ambiguous + nascent fractions)",
            "total_counts_layer": "raw total RNA counts (mature + ambiguous + nascent fractions) - same as X",
            "nascent_counts_layer": "raw nascent RNA counts (from nascent libraries)",
            "normalization": "Analysis performed on normalized subsets (target_sum=1e4 + log1p), but X contains raw counts",
            "data_type": "single-cell perturbation experiment with total and nascent RNA",
            "created_by": "preprocess_perturb.py prepare_final_dataset()"
        }

        # Remove .raw since we now have raw counts in layers
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
        log_print(f"📊 Total counts range: {adata_final.layers['total_counts'].sum(axis=1).min():.0f} - {adata_final.layers['total_counts'].sum(axis=1).max():.0f}")
        log_print(f"📊 Nascent counts range: {adata_final.layers['nascent_counts'].sum(axis=1).min():.0f} - {adata_final.layers['nascent_counts'].sum(axis=1).max():.0f}")

    else:
        # This should never happen due to the check above, but keeping for completeness
        raise RuntimeError(f"Required combined layers not found. Available layers: {list(adata.layers.keys())}. Expected: {required_combined_layers}")


    log_print(f"💾 Writing final dataset to {final_output_file}...")
    adata_final.write_h5ad(final_output_file)
    log_print("✅ Final dataset saved successfully!")

    return adata_final


def main(cutoff=5, target_genes=["TP53", "CDKN1A", "FDPS", "RABGGTA", "MRPL34", "QARS"],
         input_file=None,
         outdir="../analysis_results",
         final_output_file="./adata.ready_for_train.h5ad",
         use_gpu=False, gpu_threshold_gb=2.0, hvg_reference_file=None, target_clusters=3,
         stratification_column="condition", stratification_values=["72hr_wt", "168hr_wt"],
         non_targeting_strings=["non-targeting", "non.targeting"]):
    """Main processing pipeline for unified data format.

    Parameters:
    -----------
    cutoff : int, default=5
        Threshold for guide assignment
    target_genes : list, default=["TP53", "CDKN1A"]
        List of target genes to analyze
    input_file : str
        Path to input h5ad file with annotated unified data (includes nascent layers)
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
        adata_hvg = perform_dimensionality_reduction(adata_hvg, use_gpu=gpu_enabled)

        # Clean GPU memory after dimensionality reduction
        if gpu_enabled:
            cleanup_gpu_memory()

        adata_hvg = perform_clustering(adata_hvg, use_gpu=gpu_enabled, target_clusters=target_clusters)

        # Clean GPU memory after clustering
        if gpu_enabled:
            cleanup_gpu_memory()

        # Use expressed dataset for DE and scoring (already available from normalize_and_preprocess)

        # Transfer leiden clustering results from HVG to expressed dataset
        log_print("🔄 Transferring leiden clustering results to expressed dataset...")
        leiden_cols = [col for col in adata_hvg.obs.columns if col.startswith('leiden_')]
        for col in leiden_cols:
            adata_expressed.obs[col] = adata_hvg.obs[col]
        log_print(f"📊 Transferred {len(leiden_cols)} leiden clustering results: {leiden_cols}")

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
        scoring_config = {
            "target_genes": {
                "groups": target_genes + [["TP53", "CDKN1A"]],  # Individual genes + TP53/CDKN1A pair only
                "strategy": "target_gene_strategy",
                "stratify_by": stratification_column,
                "stratify_values": stratification_values
            },
            "guide_moi": {
                "groups": ["guide_category"],
                "strategy": "guide_strategy"
            },
            "leiden_clusters": {
                "groups": [2, 3, 4],
                "strategy": "cluster_strategy"
            }
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
        adata_expressed = calculate_unified_scores(adata_expressed, scoring_config, use_gpu=gpu_enabled, cutoff=cutoff, non_targeting_strings=non_targeting_strings)

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

        # Skip UMAP plotting in main pipeline - will be done separately
        log_print("🎨 Skipping UMAP plots in main pipeline - use plot_umap rule or --plot-only for visualization")

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
    parser.add_argument('--cutoff', type=int, default=5,
                        help='Threshold for guide assignment (default: 5)')
    parser.add_argument('--target-genes', nargs='+', default=["TP53", "CDKN1A", "FDPS", "RABGGTA", "MRPL34", "QARS"],
                        help='Target genes to analyze (default: TP53 CDKN1A FDPS RABGGTA MRPL34 QARS)')
    parser.add_argument('--input-file', required=True,
                        help='Input combined h5ad file from combine_sublibraries')
    parser.add_argument('--outdir', default="../analysis_results",
                        help='Output directory for plots and results')
    parser.add_argument('--final-output-file', default="./adata.ready_for_train.h5ad",
                        help='Path for final h5ad output ready for training')
    parser.add_argument('--force', action='store_true',
                        help='Force regeneration of output files even if they exist')
    parser.add_argument('--use-gpu', action='store_true',
                        help='Use GPU acceleration for computations (requires RAPIDS)')
    parser.add_argument('--gpu-threshold-gb', type=float, default=2.0,
                        help='Maximum used GPU memory allowed in GB (default: 2.0)')
    parser.add_argument('--hvg-reference-file', default=None,
                        help='Reference h5ad file to compute HVGs from (e.g., adata_all.annot.ultima.filt.h5ad)')
    parser.add_argument('--target-clusters', type=int, default=3,
                        help='Target number of leiden clusters to identify (default: 3)')
    parser.add_argument('--stratification-column', default="condition",
                        help='Column for stratification analysis (default: condition)')
    parser.add_argument('--stratification-values', nargs='+', default=["72hr_wt", "168hr_wt"],
                        help='Values for stratification analysis (default: 72hr_wt 168hr_wt)')
    parser.add_argument('--non-targeting-strings', nargs='+', default=["non-targeting", "non.targeting"],
                        help='Strings to identify non-targeting guides (default: non-targeting non.targeting)')
    args = parser.parse_args()

    # Note: Snakemake handles dependency tracking, so we always regenerate when called
    if os.path.exists(args.final_output_file):
        log_print(f"🔄 Final output file {args.final_output_file} exists - regenerating as requested by Snakemake...")

    # Run the pipeline
    adata_final = main(
        cutoff=args.cutoff,
        target_genes=args.target_genes,
        input_file=args.input_file,
        outdir=args.outdir,
        final_output_file=args.final_output_file,
        use_gpu=args.use_gpu,
        gpu_threshold_gb=args.gpu_threshold_gb,
        hvg_reference_file=args.hvg_reference_file,
        target_clusters=args.target_clusters,
        stratification_column=args.stratification_column,
        stratification_values=args.stratification_values,
        non_targeting_strings=args.non_targeting_strings
    )
