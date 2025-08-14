#!/usr/bin/env python3
"""
Generate QC cell filtering decisions for multiple QC methods.

This script calculates various QC filtering decisions including:
- MAD-based methods (median + 2*MAD, median + 3*MAD)
- Percentile-based methods (95th, 99th percentile)
- 2D GMM method using mitochondrial % and gene-UMI residuals

The output includes both boolean filtering decisions and continuous 
posterior probabilities for visualization.
"""

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy.stats import median_abs_deviation
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression


def calculate_2d_gmm_statistics(mito_values, umi_counts, gene_counts):
    """Calculate 2D GMM statistics using mitochondrial % and gene-UMI residuals.
    
    Args:
        mito_values: Array of mitochondrial percentages
        umi_counts: Array of UMI counts 
        gene_counts: Array of gene counts
        
    Returns:
        dict: Dictionary with 2D GMM statistics and cutoffs
    """
    # Remove NaN values
    mask = ~(np.isnan(mito_values) | np.isnan(umi_counts) | np.isnan(gene_counts))
    mito_clean = mito_values[mask]
    umi_clean = umi_counts[mask]
    gene_clean = gene_counts[mask]
    
    if len(mito_clean) < 100:  # Need minimum cells for meaningful analysis
        return {'converged': False, 'error': 'Insufficient data points'}
    
    # Log transform (using log1p to handle zeros)
    log_umis = np.log1p(umi_clean)
    log_genes = np.log1p(gene_clean)
    
    # Fit linear regression: log(genes) ~ log(UMIs)
    lr = LinearRegression()
    lr.fit(log_umis.reshape(-1, 1), log_genes)
    predicted_log_genes = lr.predict(log_umis.reshape(-1, 1))
    
    # Calculate residuals
    residuals = log_genes - predicted_log_genes
    
    # Prepare features for 2D GMM
    features = np.column_stack([mito_clean, residuals])
    
    # Standardize features
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features)
    
    # Fit 2D GMM
    gmm = GaussianMixture(n_components=2, random_state=42, covariance_type='full')
    gmm.fit(features_scaled)
    
    if not gmm.converged_:
        return {'converged': False, 'error': 'GMM did not converge'}
    
    # Get predictions
    labels = gmm.predict(features_scaled)
    posteriors = gmm.predict_proba(features_scaled)
    
    # Identify which component is "compromised" (high mito, low residual)
    # Component with higher mean mito (in scaled space) is likely compromised
    means_original = scaler.inverse_transform(gmm.means_)
    comp0_mito = means_original[0, 0]
    comp1_mito = means_original[1, 0]
    
    if comp0_mito > comp1_mito:
        compromised_idx = 0
        healthy_idx = 1
    else:
        compromised_idx = 1
        healthy_idx = 0
    
    # Calculate statistics
    n_compromised = np.sum(labels == compromised_idx)
    n_healthy = np.sum(labels == healthy_idx)
    
    # Get posterior probability thresholds
    prob_compromised = posteriors[:, compromised_idx]
    
    # Find cells at different posterior thresholds
    thresholds = {}
    for thresh in [0.5, 0.75, 0.9]:
        n_above = np.sum(prob_compromised > thresh)
        thresholds[f'posterior_{int(thresh*100)}'] = {
            'n_cells_compromised': int(n_above),
            'pct_cells_compromised': float(100 * n_above / len(prob_compromised))
        }
    
    results = {
        'converged': True,
        'regression': {
            'slope': float(lr.coef_[0]),
            'intercept': float(lr.intercept_),
            'r_squared': float(lr.score(log_umis.reshape(-1, 1), log_genes))
        },
        'scaler': {
            'mito_mean': float(scaler.mean_[0]),
            'mito_std': float(scaler.scale_[0]),
            'residual_mean': float(scaler.mean_[1]),
            'residual_std': float(scaler.scale_[1])
        },
        'components': {
            'compromised': {
                'idx': int(compromised_idx),
                'mito_mean': float(means_original[compromised_idx, 0]),
                'residual_mean': float(means_original[compromised_idx, 1]),
                'weight': float(gmm.weights_[compromised_idx]),
                'n_cells': int(n_compromised)
            },
            'healthy': {
                'idx': int(healthy_idx),
                'mito_mean': float(means_original[healthy_idx, 0]),
                'residual_mean': float(means_original[healthy_idx, 1]),
                'weight': float(gmm.weights_[healthy_idx]),
                'n_cells': int(n_healthy)
            }
        },
        'thresholds': thresholds,
        'bic': float(gmm.bic(features_scaled)),
        'aic': float(gmm.aic(features_scaled)),
        'posteriors': posteriors,  # Return full posterior array
        'compromised_idx': compromised_idx
    }
    
    return results


def save_qc_cell_lists(adata, output_path, config=None):
    """Save cell filtering decisions for multiple QC methods.
    
    Args:
        adata: AnnData object with QC metrics
        output_path: Path to save TSV file with cell lists
        config: Configuration dictionary for fallback options
    """
    print("\nCalculating QC filtering decisions for all methods...")
    
    # Initialize DataFrame with barcodes
    cell_lists = pd.DataFrame({'barcode': adata.obs.index})
    
    # Get mitochondrial values
    mito_values = adata.obs['pct_counts_mt'].values
    
    # Calculate basic statistics for simple methods
    median = np.median(mito_values)
    mad = median_abs_deviation(mito_values)
    
    # MAD-based methods
    cell_lists['median_plus_2mad'] = mito_values < (median + 2 * mad)
    cell_lists['median_plus_3mad'] = mito_values < (median + 3 * mad)
    
    # Percentile-based methods
    cell_lists['percentile_95'] = mito_values < np.percentile(mito_values, 95)
    cell_lists['percentile_99'] = mito_values < np.percentile(mito_values, 99)
    
    # 2D GMM method (DEFAULT)
    if 'total_counts' in adata.obs.columns and 'n_genes' in adata.obs.columns:
        # Calculate 2D GMM
        gmm_2d_stats = calculate_2d_gmm_statistics(
            mito_values,
            adata.obs['total_counts'].values,
            adata.obs['n_genes'].values
        )
        
        if gmm_2d_stats['converged']:
            # Get the posterior probabilities from the returned stats
            posteriors = gmm_2d_stats['posteriors']
            compromised_idx = gmm_2d_stats['compromised_idx']
            
            # Get probability of being compromised
            prob_compromised = posteriors[:, compromised_idx]
            
            # Set filtering decisions (True = keep cell, False = filter out)
            cell_lists['gmm_2d_posterior_75'] = prob_compromised < 0.75
            cell_lists['gmm_2d_posterior_50'] = prob_compromised < 0.50
            cell_lists['gmm_2d_posterior_90'] = prob_compromised < 0.90
            
            # IMPORTANT: Save the continuous posterior probability for plotting
            cell_lists['gmm_2d_posterior_prob'] = prob_compromised
            
        else:
            # GMM failed - raise error if GMM is the requested method
            error_msg = gmm_2d_stats.get('error', 'Unknown error')
            
            # Check if GMM is the requested QC method in config
            requested_qc_method = config.get('combined_sublibraries', {}).get('qc_method', 'gmm_2d_posterior_75')
            
            if 'gmm' in requested_qc_method.lower():
                # GMM is requested but failed - raise error
                raise RuntimeError(
                    f"2D GMM failed to converge: {error_msg}\n"
                    f"The requested QC method '{requested_qc_method}' cannot be used.\n"
                    f"Please change config['combined_sublibraries']['qc_method'] to one of:\n"
                    f"  - median_plus_2mad\n"
                    f"  - median_plus_3mad\n"
                    f"  - percentile_95\n"
                    f"  - percentile_99\n"
                )
            else:
                # GMM failed but user is using a different method anyway
                print(f"  Note: 2D GMM failed to converge ({error_msg}), but this is OK")
                print(f"  since you're using '{requested_qc_method}' for QC filtering")
            
    else:
        raise ValueError("Missing required data (total_counts or n_genes) for 2D GMM calculation")
    
    # Save to TSV
    cell_lists.to_csv(output_path, sep='\t', index=False)
    
    # Print summary
    print("\nQC filtering summary:")
    for col in cell_lists.columns:
        if col not in ['barcode', 'gmm_2d_posterior_prob']:  # Skip continuous column
            n_keep = cell_lists[col].sum()
            n_total = len(cell_lists)
            pct_keep = 100 * n_keep / n_total
            default_marker = " (DEFAULT)" if col == 'gmm_2d_posterior_75' else ""
            print(f"  {col}{default_marker}: {n_keep:,}/{n_total:,} cells ({pct_keep:.1f}%)")
    
    print(f"\nSaved cell lists to: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Generate QC cell filtering decisions')
    parser.add_argument('--h5ad', required=True, help='Annotated h5ad file')
    parser.add_argument('--output', required=True, help='Output TSV file for cell lists')
    parser.add_argument('--config', help='Config YAML file')
    
    args = parser.parse_args()
    
    # Load config if provided
    config = None
    if args.config:
        with open(args.config) as f:
            config = yaml.safe_load(f)
    
    # Load h5ad file
    print(f"Loading h5ad file: {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    print(f"Loaded {adata.n_obs} barcodes x {adata.n_vars} genes")
    
    # Ensure we have the necessary QC metrics
    required_columns = ['pct_counts_mt', 'total_counts', 'n_genes_by_counts']
    missing_columns = [col for col in required_columns if col not in adata.obs.columns]
    
    # Handle n_genes vs n_genes_by_counts
    if 'n_genes_by_counts' in adata.obs.columns and 'n_genes' not in adata.obs.columns:
        adata.obs['n_genes'] = adata.obs['n_genes_by_counts']
    
    # Re-check required columns
    required_columns = ['pct_counts_mt', 'total_counts', 'n_genes']
    missing_columns = [col for col in required_columns if col not in adata.obs.columns]
    
    if missing_columns:
        raise ValueError(f"Required QC metrics missing from h5ad file: {missing_columns}")
    
    # Generate and save cell lists
    save_qc_cell_lists(adata, args.output, config)


if __name__ == "__main__":
    main()