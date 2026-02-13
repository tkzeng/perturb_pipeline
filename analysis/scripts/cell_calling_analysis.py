#!/usr/bin/env python3
"""
Cell calling analysis using multiple approaches:
1. Expected cell count (top N barcodes)
2. Simple UMI threshold
3. EmptyDrops (via R script)

This script performs cell calling analysis and saves results.
Plotting is handled separately by cell_calling_plots.py.
"""

import argparse
import subprocess
import sys
import yaml
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import scanpy as sc
import os

def load_h5ad_matrix(h5ad_file):
    """Load count matrix from kallisto h5ad output."""
    adata = sc.read_h5ad(h5ad_file)
    
    # Delete guide data from obsm to save memory (not needed for cell calling)
    if 'guide_counts' in adata.obsm:
        print(f"Removing guide_counts matrix from memory ({adata.obsm['guide_counts'].shape})")
        del adata.obsm['guide_counts']
    if 'guide_assignment' in adata.obsm:
        print(f"Removing guide_assignment matrix from memory")
        del adata.obsm['guide_assignment']
    
    # Sum layers to guarantee X is total counts
    adata.X = adata.layers['mature'] + adata.layers['nascent'] + adata.layers['ambiguous']
    
    # Calculate basic QC metrics if not already present
    if 'total_counts' not in adata.obs.columns:
        adata.var['mt'] = adata.var_names.str.match('(?i)^mt-') # ignore capital letter for both human and mouse
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    return adata

def run_dropletutils_r(kb_dir, sample_id, output_dir, emptydrops_lower=100, ncores=1, expected_cells=0, 
                      run_emptydrops=False, run_barcoderanks=False, fdr_cutoffs=None,
                      run_by_group=False, plot_dir = False):
    """Run R script for DropletUtils analysis (EmptyDrops and/or BarcodeRanks).
    
    NOTE: EmptyDrops is DEPRECATED. Use BarcodeRanks methods instead.
    EmptyDrops code is preserved for backwards compatibility but is not actively maintained.
    """
    # Get the path to the R script relative to this Python script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path = os.path.join(script_dir, "dropletutils_r.R")
    
    cmd = [
        "Rscript", "--vanilla", r_script_path,  # Use --vanilla to ensure only conda R packages are used
        "--kb_dir", str(kb_dir),
        "--sample_id", sample_id,
        "--output_dir", str(output_dir),
        "--emptydrops_lower", str(emptydrops_lower),
        "--ncores", str(ncores),
        "--expected_cells", str(expected_cells),
        "--plot_dir", str(plot_dir)
    ]
    
    if run_emptydrops:
        cmd.append("--run_emptydrops")
        if fdr_cutoffs:
            # Pass FDR cutoffs as comma-separated string
            cmd.extend(["--fdr_cutoffs", ",".join(map(str, fdr_cutoffs))])
    if run_barcoderanks:
        cmd.append("--run_barcoderanks")
    if run_by_group:
        cmd.append("--run_by_group")

    result = subprocess.run(cmd, check=True)
    print(f"DropletUtils R script completed successfully")

def expected_cell_method(adata, expected_cells):
    """Cell calling using expected cell count - take top N barcodes."""
    if 'total_counts' in adata.obs.columns:
        total_counts = adata.obs['total_counts'].values
    else:
        total_counts = np.array(adata.X.sum(axis=1)).flatten()
    
    sorted_indices = np.argsort(total_counts)[::-1]
    
    # Take top N barcodes
    n_cells = min(expected_cells, len(sorted_indices))
    cell_indices = sorted_indices[:n_cells]
    
    # Create boolean mask
    is_cell = np.zeros(adata.n_obs, dtype=bool)
    is_cell[cell_indices] = True
    
    return is_cell, total_counts[sorted_indices[n_cells-1]] if n_cells > 0 else 0

def threshold_method(adata, min_umi_threshold):
    """Cell calling using simple UMI threshold."""
    if 'total_counts' in adata.obs.columns:
        total_counts = adata.obs['total_counts'].values
    else:
        total_counts = np.array(adata.X.sum(axis=1)).flatten()
    
    is_cell = total_counts >= min_umi_threshold
    return is_cell, min_umi_threshold

def load_dropletutils_results(output_dir, sample_id, adata):
    """Load DropletUtils results from R script output (EmptyDrops + barcodeRanks)."""
    # Load DropletUtils results
    results_file = Path(output_dir) / "dropletutils_results.tsv"
    df = pd.read_csv(results_file, sep='\t')
    
    # Map back to adata
    barcode_to_idx = {bc: i for i, bc in enumerate(adata.obs_names)}
    
    # Load results for all methods
    results = {}
    
    # 1. Load EmptyDrops results from summary file to get actual FDR thresholds tested
    summary_file = Path(output_dir) / "cell_calling_summary.tsv"
    if summary_file.exists():
        summary_df = pd.read_csv(summary_file, sep='\t')
        
        # Process EmptyDrops results for each FDR threshold
        emptydrops_rows = summary_df[summary_df['method'].str.startswith('EmptyDrops_FDR')]
        for _, row in emptydrops_rows.iterrows():
            method_name = row['method']
            fdr_threshold = row['fdr_threshold']
            
            is_cell = np.zeros(adata.n_obs, dtype=bool)
            
            # Get cells below FDR threshold
            cells_df = df[df['fdr'] <= fdr_threshold]
            for barcode in cells_df['barcode']:
                if barcode in barcode_to_idx:
                    is_cell[barcode_to_idx[barcode]] = True
            
            results[method_name] = (is_cell, fdr_threshold)
    
    # 2. BarcodeRanks knee and inflection points
    # The R script saves a summary file with these values
    summary_file = Path(output_dir) / "cell_calling_summary.tsv"
    if summary_file.exists():
        summary_df = pd.read_csv(summary_file, sep='\t')
        
        # Get knee and inflection thresholds
        knee_row = summary_df[summary_df['method'] == 'BarcodeRanks_Knee']
        inflection_row = summary_df[summary_df['method'] == 'BarcodeRanks_Inflection']
        
        if not knee_row.empty:
            knee_threshold = knee_row.iloc[0]['threshold_used']
            # Get cells above knee threshold
            if 'total_counts' in adata.obs.columns:
                total_counts = adata.obs['total_counts'].values
            else:
                total_counts = np.array(adata.X.sum(axis=1)).flatten()
            
            is_cell_knee = total_counts >= knee_threshold
            results['BarcodeRanks_Knee'] = (is_cell_knee, knee_threshold)
        
        if not inflection_row.empty:
            inflection_threshold = inflection_row.iloc[0]['threshold_used']
            # Get cells above inflection threshold
            if 'total_counts' in adata.obs.columns:
                total_counts = adata.obs['total_counts'].values
            else:
                total_counts = np.array(adata.X.sum(axis=1)).flatten()
            
            is_cell_inflection = total_counts >= inflection_threshold
            results['BarcodeRanks_Inflection'] = (is_cell_inflection, inflection_threshold)
    
    return results

def load_dropletutils_results_by_group(output_dir, sample_id, adata):
    """
    Load per-group BarcodeRanks results produced by R script.

    Expected files:
      - dropletutils_results_by_group.tsv
      - cell_calling_summary_by_group.tsv (optional but nice for thresholds)

    Returns:
      methods_by_group: defaultdict(dict)
    """
    output_dir = Path(output_dir)
    by_group_file = output_dir / "dropletutils_results_by_group.tsv"
    if not by_group_file.exists():
        return defaultdict(dict)  # nothing

    df = pd.read_csv(by_group_file, sep="\t")

    # map barcode to global index in adata
    barcode_to_idx = {bc: i for i, bc in enumerate(adata.obs_names)}

    methods_by_group = defaultdict(dict)

    # Use per-row threshold columns (knee_threshold / inflection_threshold) if present
    # (they are repeated per row in that group; we'll take the first non-null)
    for g, gdf in df.groupby("group", sort=False):
        # global masks
        mask_knee = np.zeros(adata.n_obs, dtype=bool)
        mask_inf  = np.zeros(adata.n_obs, dtype=bool)

        # thresholds
        knee_thr = None
        inf_thr  = None

        if "knee_threshold" in gdf.columns:
            vals = gdf["knee_threshold"].dropna().unique()
            knee_thr = float(vals[0]) if len(vals) else np.nan
        else:
            knee_thr = np.nan

        if "inflection_threshold" in gdf.columns:
            vals = gdf["inflection_threshold"].dropna().unique()
            inf_thr = float(vals[0]) if len(vals) else np.nan
        else:
            inf_thr = np.nan

        # fill masks by barcode
        # guard: columns names in TSV
        has_knee_col = "is_cell_knee" in gdf.columns
        has_inf_col  = "is_cell_inflection" in gdf.columns
        if not (has_knee_col and has_inf_col):
            # If missing, we can't reconstruct masks reliably
            # Return empty for this group
            continue

        for _, row in gdf.iterrows():
            bc = row["barcode"]
            idx = barcode_to_idx.get(bc, None)
            if idx is None:
                continue
            if bool(row["is_cell_knee"]):
                mask_knee[idx] = True
            if bool(row["is_cell_inflection"]):
                mask_inf[idx] = True

        methods_by_group[g]["BarcodeRanks_Knee"] = (mask_knee, knee_thr)
        methods_by_group[g]["BarcodeRanks_Inflection"] = (mask_inf, inf_thr)

    return methods_by_group

def calculate_umis_in_cells_pct(adata, methods_results):
    """Calculate percentage of total UMIs that fall in called cells."""
    if 'total_counts' in adata.obs.columns:
        total_counts = adata.obs['total_counts'].values
    else:
        total_counts = np.array(adata.X.sum(axis=1)).flatten()
    
    total_umis = np.sum(total_counts)
    
    percentages = {}
    for method, (is_cell, threshold) in methods_results.items():
        umis_in_cells = np.sum(total_counts[is_cell])
        percentages[method] = (umis_in_cells / total_umis * 100) if total_umis > 0 else 0.0
    
    return percentages

def calculate_umis_in_cells_pct_by_group(adata, methods_by_group, group_col):
    if "total_counts" in adata.obs.columns:
        total_counts = adata.obs["total_counts"].to_numpy()
    else:
        total_counts = np.array(adata.X.sum(axis=1)).flatten()

    out = defaultdict(dict)

    for g, mdict in methods_by_group.items():
        g_mask = (adata.obs[group_col].values == g)         # global mask
        total_umis_g = float(total_counts[g_mask].sum())

        for method, (is_cell_global, _thr) in mdict.items():
            is_cell_global = np.asarray(is_cell_global, dtype=bool)  # global length
            umis_in_cells_g = float(total_counts[g_mask & is_cell_global].sum())
            out[g][method] = (umis_in_cells_g / total_umis_g * 100.0) if total_umis_g > 0 else 0.0

    return out
    
def get_cell_calling_params(config, sample_id):
    """Get cell calling parameters for a specific sample."""
    
    # Load sample-specific parameters from sample_info file (REQUIRED)
    sample_info_file = config['sample_info_file']
    if not os.path.exists(sample_info_file):
        raise FileNotFoundError(f"Sample info file not found: {sample_info_file}")
    
    sample_df = pd.read_excel(sample_info_file)
    
    sample_row = sample_df[sample_df['sample_id'] == sample_id]
    if sample_row.empty:
        raise ValueError(f"Sample {sample_id} not found in {sample_info_file}")
    
    row = sample_row.iloc[0]
    
    # REQUIRED: expected_cells must be specified
    if 'expected_cells' not in sample_df.columns or pd.isna(row['expected_cells']):
        raise ValueError(f"Expected cells must be specified for sample {sample_id} in {sample_info_file}")
    
    params = {'expected_cells': int(row['expected_cells'])}
    
    # Optional parameters with config defaults
    defaults = config['cell_calling']['defaults']
    if pd.notna(row.get('min_umi_threshold')):
        params['min_umi_threshold'] = int(row['min_umi_threshold'])
    else:
        params['min_umi_threshold'] = defaults.get('min_umi_threshold', 100)
    
    params['emptydrops_lower'] = defaults.get('emptydrops_lower', 100)

    # add cell calling by self-defined group info
    groupby_col = defaults.get('cell_calling_groupby_col')
    if groupby_col is not None:
        params['cell_calling_groupby_col'] = groupby_col
        params['expected_cells_per_biosample'] = defaults.get('expected_cells_per_biosample', 100)
        params['min_umi_threshold_per_biosample'] = defaults.get('min_umi_threshold_per_biosample', 100)
    else:
        params['cell_calling_groupby_col'] = None
    return params

def main():
    parser = argparse.ArgumentParser(description="Cell calling analysis")
    parser.add_argument("--h5ad_file", required=True, help="Kallisto h5ad output file")
    parser.add_argument("--kb_dir", required=True, help="Kallisto bus output directory")
    parser.add_argument("--sample-id", required=True, help="Full sample ID (format: pool:sample)")
    parser.add_argument("--config", required=True, help="Config YAML file")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--plot_dir", required=True, help="Output directory")
    parser.add_argument("--ncores", type=int, required=True, help="Number of cores for parallel processing")
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)
    
    # Get sample-specific parameters
    params = get_cell_calling_params(config, args.sample_id)
    
    # Load h5ad matrix
    print(f"📥 Loading h5ad file: {args.h5ad_file}")
    adata = load_h5ad_matrix(args.h5ad_file)
    
    print(f"Loaded matrix: {adata.n_obs} barcodes x {adata.n_vars} genes")
    
    # Apply different cell calling methods
    methods_results = {}
    methods_results_by_biosample = defaultdict(dict)  # group -> {method: (is_cell, threshold)}

    # Get list of methods to run from config
    methods_to_run = config['cell_calling']['methods_to_run']
    print(f"Running cell calling methods: {methods_to_run}")
    
    # check if adata actually contains target groupby col
    group_col = params['cell_calling_groupby_col']
    if group_col is not None and group_col not in adata.obs.columns:
        print(f"[WARN] groupby col '{group_col}' not in adata.obs; fallback to per-sample calling.")
        group_col = None
    # write barcode with biological group information
    if group_col is not None:
        # filtering barcode with NA in group_col
        group_values = adata.obs[group_col].copy()
        valid_mask = group_values.notna() & (group_values.astype(str) != "nan") & (group_values.astype(str) != "")
        groups = sorted(group_values[valid_mask].unique().tolist())
        print(f"cell calling by group enabled: {group_col}, n_groups={len(groups)}")
        # write barcode with corresponding sample for cell calling
        metadata = adata.obs.loc[valid_mask, [group_col]].copy()
        metadata.insert(0, "barcode", metadata.index)  # ensure barcode column exists
        metadata = metadata.rename(columns={group_col: "cell_calling_by_group"})
        metadata = metadata[["barcode", "cell_calling_by_group"]]
        counts_dir = Path(args.kb_dir) / "counts_filtered"
        out_csv = counts_dir / "barcode_with_group.csv"
        metadata.to_csv(out_csv, index=False)


    # Method 1: Expected cell count
    if 'Expected_Cells' in methods_to_run:
        print(f"Running expected cell method (n={params['expected_cells']})...")
        methods_results['Expected_Cells'] = expected_cell_method(adata, params['expected_cells'])
        
        if group_col is not None:
            for g in groups:
                g_mask = (adata.obs[group_col].values == g)   # global bool mask, len = adata.n_obs
                adata_g = adata[g_mask]                       # local view

                local_mask, thr = expected_cell_method(adata_g, params['expected_cells_per_biosample'])
                global_mask = np.zeros(adata.n_obs, dtype=bool)
                global_mask[g_mask] = local_mask

                methods_results_by_biosample[g]["Expected_Cells"] = (global_mask, thr)

    # Method 2: Threshold method
    if 'UMI_Threshold' in methods_to_run:
        print(f"Running threshold method (min_umi={params['min_umi_threshold']})...")
        methods_results['UMI_Threshold'] = threshold_method(adata, params['min_umi_threshold'])
    
        if group_col is not None:
            for g in groups:
                g_mask = (adata.obs[group_col].values == g)   # global bool mask, len = adata.n_obs
                adata_g = adata[g_mask]                       # local view

                local_mask, thr = threshold_method(adata_g, params["min_umi_threshold_per_biosample"])
                global_mask = np.zeros(adata.n_obs, dtype=bool)
                global_mask[g_mask] = local_mask
                
                methods_results_by_biosample[g]['UMI_Threshold'] = (global_mask, thr)

    # flag to detect if setting run by group
    run_by_group = (group_col is not None)

    # Method 3-7: DropletUtils methods (via R)
    # Parse EmptyDrops methods and their FDR cutoffs
    emptydrops_fdr_cutoffs = []
    for method in methods_to_run:
        if method.startswith('EmptyDrops_FDR_'):
            # Parse FDR value from method name (e.g., EmptyDrops_FDR_0.001 -> 0.001)
            fdr_str = method.replace('EmptyDrops_FDR_', '')
            fdr_value = float(fdr_str)
            emptydrops_fdr_cutoffs.append(fdr_value)
    
    barcoderanks_methods = ['BarcodeRanks_Knee', 'BarcodeRanks_Inflection']
    
    run_emptydrops = len(emptydrops_fdr_cutoffs) > 0
    run_barcoderanks = any(method in methods_to_run for method in barcoderanks_methods)
    
    if run_emptydrops or run_barcoderanks:
        print("Running DropletUtils analysis on pre-filtered matrix...")
        if run_emptydrops:
            print(f"  - EmptyDrops enabled with FDR cutoffs: {emptydrops_fdr_cutoffs}")
            print("  WARNING: EmptyDrops is DEPRECATED. Use BarcodeRanks methods instead.")
        if run_barcoderanks:
            print("  - BarcodeRanks enabled")
            
        run_dropletutils_r(args.kb_dir, args.sample_id, output_dir, 
                          params['emptydrops_lower'], 
                          args.ncores,
                          params['expected_cells'],
                          run_emptydrops=run_emptydrops,
                          run_barcoderanks=run_barcoderanks,
                          fdr_cutoffs=emptydrops_fdr_cutoffs if run_emptydrops else None,
                          run_by_group=run_by_group,
                          plot_dir=args.plot_dir
                          )
        
        dropletutils_results = load_dropletutils_results(output_dir, args.sample_id, adata)
        # Only include methods that were requested
        for method in methods_to_run:
            if method in dropletutils_results:
                methods_results[method] = dropletutils_results[method]
    
        if run_by_group:
            dropletutils_results_by_group = load_dropletutils_results_by_group(output_dir, args.sample_id, adata)
            for g, mdict in dropletutils_results_by_group.items():
                methods_results_by_biosample[g].update(mdict)

    # Calculate percentage of UMIs in cells
    umi_percentages = calculate_umis_in_cells_pct(adata, methods_results)

    if run_by_group and group_col is not None:
        umi_percentages_by_biosample = calculate_umis_in_cells_pct_by_group(
            adata,
            methods_results_by_biosample,   # group -> {method: (global_mask, thr)}
            group_col
        )

    # Save results
    results_data = []
    for method, (is_cell, threshold) in methods_results.items():
        n_cells = np.sum(is_cell)
        umis_in_cells_pct = umi_percentages[method]
        
        results_data.append({
            'sample_id': args.sample_id,
            'method': method,
            'n_cells_called': n_cells,
            'threshold_used': threshold,
            'umis_in_cells_pct': umis_in_cells_pct,
            'expected_cells_param': params['expected_cells'],
            'min_umi_threshold_param': params['min_umi_threshold']
        })
    
    # Save to TSV
    results_df = pd.DataFrame(results_data)
    results_df.to_csv(output_dir / 'results.tsv', 
                      sep='\t', index=False)
    
    if group_col is not None:
        rows = []
        for g, mdict in methods_results_by_biosample.items():
            for method, (mask, thr) in mdict.items():
                rows.append({
                    "sample_id": args.sample_id,
                    "group": g,
                    "method": method,
                    "n_cells_called": int(np.sum(mask)),
                    "threshold_used": thr,
                    "umis_in_cells_pct": umi_percentages_by_biosample.get(g, {}).get(method, np.nan),
                })

        pd.DataFrame(rows).to_csv(output_dir / "results_by_group.tsv", sep="\t", index=False)
        
    # Save cell barcodes for each method
    for method, (is_cell, threshold) in methods_results.items():
        cell_barcodes = adata.obs_names[is_cell]
        with open(output_dir / f'{args.sample_id}_{method}_cell_barcodes.txt', 'w') as f:
            for barcode in cell_barcodes:
                f.write(f"{barcode}\n")
    
    # Save cell barcodes for each method (by group)
    if group_col is not None and len(methods_results_by_biosample) > 0:
        result_per_group_dir = Path(output_dir) / "cell_calling_result_per_group"
        result_per_group_dir.mkdir(parents=True, exist_ok=True)

        for g, mdict in methods_results_by_biosample.items():
            for method, (is_cell_global, thr) in mdict.items():
                is_cell_global = np.asarray(is_cell_global, dtype=bool)  # global length
                cell_barcodes = adata.obs_names[is_cell_global]

                out_txt = result_per_group_dir / f"{args.sample_id}_{g}_{method}_cell_barcodes.txt"
                with open(out_txt, "w") as f:
                    f.write("\n".join(cell_barcodes) + ("\n" if len(cell_barcodes) else ""))

    print(f"Cell calling analysis completed. Results saved to {output_dir}")
    
    # Print summary
    print("\nSummary:")
    for method, (is_cell, threshold) in methods_results.items():
        n_cells = np.sum(is_cell)
        umis_pct = umi_percentages[method]
        print(f"  {method}: {n_cells} cells (threshold: {threshold:.0f}, "
              f"UMIs in cells: {umis_pct:.1f}%)")
    if group_col is not None and len(methods_results_by_biosample) > 0:
        
        print(f"\nSummary (by group sum: {group_col}):")

        # all methods that appear in any group
        methods_in_groups = sorted({m for md in methods_results_by_biosample.values() for m in md.keys()})

        for method in methods_in_groups:
            total_sum = 0
            n_groups_used = 0

            for g, mdict in methods_results_by_biosample.items():
                if method not in mdict:
                    continue
                mask_g, _thr = mdict[method]
                mask_g = np.asarray(mask_g, dtype=bool)  # global length
                total_sum += int(mask_g.sum())
                n_groups_used += 1

            print(f"  {method}: total_called(sum over groups)={total_sum}, n_groups_used={n_groups_used}")
        
    print("\nWriting cell calling result to adata")

    for method, (is_cell, thr) in methods_results.items():
        adata.obs[f"is_cell_{method}"] = np.asarray(is_cell, dtype=bool)
        adata.uns[f"threshold_{method}"] = float(thr) if thr is not None else np.nan

    if group_col is not None and len(methods_results_by_biosample) > 0:
        for g, mdict in methods_results_by_biosample.items():
            g_mask = (adata.obs[group_col].values == g)
            for method, (is_cell_global, thr) in mdict.items():
                col = f"is_cell_{method}_by_group"
                if col not in adata.obs:
                    adata.obs[col] = False
                adata.obs.loc[g_mask, col] = np.asarray(is_cell_global, dtype=bool)[g_mask]

                key = f"threshold_{method}_by_group"
                if key not in adata.uns:
                    adata.uns[key] = {}
                adata.uns[key][str(g)] = float(thr) if thr is not None else np.nan
    
    adata_w_guide = sc.read_h5ad(args.h5ad_file)
    if adata_w_guide.n_obs != adata.n_obs:
        raise RuntimeError(f"n_obs mismatch: original={adata_w_guide.n_obs}, working={adata.n_obs}")
    if not np.all(adata_w_guide.obs_names.values == adata.obs_names.values):
        raise RuntimeError("obs_names order mismatch between original and working adata. Refuse to write results.")

    adata.obsm['guide_counts'] = adata_w_guide.obsm['guide_counts']
    adata.obsm['guide_assignment'] = adata_w_guide.obsm['guide_assignment']
    adata_out_path = (output_dir / f"{args.sample_id}.cell_calling.h5ad").resolve()
    adata.write_h5ad(adata_out_path)
    
    print(f"\n✅ Writing cell-calling AnnData with guide info to: {adata_out_path}")

    # Clean up memory
    del adata
    import gc
    gc.collect()

if __name__ == "__main__":
    main()