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
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
import os


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
    # Order: [residuals, mito] to match plotting convention
    features = np.column_stack([residuals, mito_clean])
    
    # Standardize features
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features)
    
    # Fit 2D GMM
    gmm = GaussianMixture(n_components=2, random_state=42, covariance_type='full')
    gmm.fit(features_scaled)
    
    # Get predictions (even if not converged)
    labels = gmm.predict(features_scaled)
    posteriors = gmm.predict_proba(features_scaled)
    
    # Identify which component is "compromised" (high mito, low residual)
    # Component with higher mean mito is likely compromised
    # Note: features are now [residuals, mito] so mito is index 1
    means_original = scaler.inverse_transform(gmm.means_)
    comp0_mito = means_original[0, 1]  # mito is second feature
    comp1_mito = means_original[1, 1]  # mito is second feature
    
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
        'converged': bool(gmm.converged_),
        'regression': {
            'slope': float(lr.coef_[0]),
            'intercept': float(lr.intercept_),
            'r_squared': float(lr.score(log_umis.reshape(-1, 1), log_genes))
        },
        'scaler': {
            'residual_mean': float(scaler.mean_[0]),  # residual is first feature now
            'residual_std': float(scaler.scale_[0]),
            'mito_mean': float(scaler.mean_[1]),  # mito is second feature now
            'mito_std': float(scaler.scale_[1])
        },
        'components': {
            'compromised': {
                'idx': int(compromised_idx),
                'residual_mean': float(means_original[compromised_idx, 0]),  # residual is first
                'mito_mean': float(means_original[compromised_idx, 1]),  # mito is second
                'weight': float(gmm.weights_[compromised_idx]),
                'n_cells': int(n_compromised)
            },
            'healthy': {
                'idx': int(healthy_idx),
                'residual_mean': float(means_original[healthy_idx, 0]),  # residual is first
                'mito_mean': float(means_original[healthy_idx, 1]),  # mito is second
                'weight': float(gmm.weights_[healthy_idx]),
                'n_cells': int(n_healthy)
            }
        },
        'thresholds': thresholds,
        'bic': float(gmm.bic(features_scaled)),
        'aic': float(gmm.aic(features_scaled)),
        'posteriors': posteriors,  # Return full posterior array
        'compromised_idx': compromised_idx,
        'gmm_model': gmm,  # Return the GMM model for plotting
        'scaler': scaler  # Return the scaler for plotting
    }
    
    return results


def save_qc_cell_lists(adata, output_path, config=None, save_plots=True, plot_output=None):
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
    
    # 2D GMM method (DEFAULT)
    if 'total_counts' in adata.obs.columns and 'n_genes' in adata.obs.columns:
        # Calculate 2D GMM
        gmm_2d_stats = calculate_2d_gmm_statistics(
            mito_values,
            adata.obs['total_counts'].values,
            adata.obs['n_genes'].values
        )
        
        # Always get the posterior probabilities
        posteriors = gmm_2d_stats['posteriors']
        compromised_idx = gmm_2d_stats['compromised_idx']
        
        # Get probability of being compromised
        prob_compromised = posteriors[:, compromised_idx]
        
        # Save the continuous posterior probability
        cell_lists['gmm_posterior_prob_compromised'] = prob_compromised
        
        # Create debug plots if requested
        if save_plots and plot_output:
                # Use the provided plot output path
                plot_dir = os.path.dirname(plot_output)
                os.makedirs(plot_dir, exist_ok=True)
                plot_path = plot_output
                
                print(f"\nGenerating combined debug plot...")
                
                # Get the regression residuals (recalculate for plotting)
                log_umis = np.log1p(adata.obs['total_counts'].values)
                log_genes = np.log1p(adata.obs['n_genes'].values)
                
                # Fit regression again for residuals
                lr = LinearRegression()
                lr.fit(log_umis.reshape(-1, 1), log_genes)
                predicted_log_genes = lr.predict(log_umis.reshape(-1, 1))
                residuals = log_genes - predicted_log_genes
                
                # Create combined figure with 1x3 layout (no 3D plot)
                fig, axes = plt.subplots(1, 3, figsize=(18, 6))
                
                # Plot 1: UMIs vs Genes BEFORE regression (left)
                ax1 = axes[0]
                scatter1 = ax1.scatter(
                    log_umis,
                    log_genes,
                    c=prob_compromised,
                    cmap='RdYlBu_r',
                    s=3,
                    alpha=0.6,
                    edgecolors='none'
                )
                ax1.plot(log_umis, predicted_log_genes, 'g-', alpha=0.5, linewidth=2, label='Regression line')
                ax1.set_xlabel('Log(Total UMIs)', fontsize=10)
                ax1.set_ylabel('Log(Number of genes)', fontsize=10)
                ax1.set_title('Before Regression\n(Gene-UMI Relationship)', fontsize=11)
                ax1.legend(fontsize=8)
                plt.colorbar(scatter1, ax=ax1, label='P(Compromised)')
                
                # Plot 2: UMIs vs Residuals AFTER regression (middle)
                ax2 = axes[1]
                scatter2 = ax2.scatter(
                    log_umis,
                    residuals,
                    c=prob_compromised,
                    cmap='RdYlBu_r',
                    s=3,
                    alpha=0.6,
                    edgecolors='none'
                )
                ax2.axhline(y=0, color='g', linestyle='-', alpha=0.5, linewidth=1)
                ax2.set_xlabel('Log(Total UMIs)', fontsize=10)
                ax2.set_ylabel('Gene-UMI Residuals', fontsize=10)
                ax2.set_title('After Regression\n(Residuals vs UMIs)', fontsize=11)
                plt.colorbar(scatter2, ax=ax2, label='P(Compromised)')
                
                # Plot 3: The actual 2D GMM features with contours - Mitochondrial % vs Residuals (right)
                ax3 = axes[2]
                scatter3 = ax3.scatter(
                    residuals,
                    mito_values,
                    c=prob_compromised,
                    cmap='RdYlBu_r',
                    s=3,
                    alpha=0.6,
                    edgecolors='none'
                )
                
                # Add GMM contour ellipses
                gmm_model = gmm_2d_stats['gmm_model']
                scaler = gmm_2d_stats['scaler']
                
                # Draw confidence ellipses for each component
                def draw_ellipse(ax, mean, cov, color, label, n_std=2.0):
                    """Draw confidence ellipse for a 2D Gaussian."""
                    # Get eigenvalues and eigenvectors
                    eigenvalues, eigenvectors = np.linalg.eigh(cov)
                    order = eigenvalues.argsort()[::-1]
                    eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
                    
                    # Calculate angle and dimensions
                    angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))
                    width, height = 2 * n_std * np.sqrt(eigenvalues)
                    
                    # Create and add ellipse
                    ellipse = Ellipse(mean, width, height, angle=angle,
                                     facecolor='none', edgecolor=color, 
                                     linewidth=2, linestyle='-', alpha=0.8,
                                     label=label)
                    ax.add_patch(ellipse)
                    return ellipse
                
                # Get component parameters in original space
                for comp_idx, (comp_name, comp_color) in enumerate([('healthy', 'blue'), ('compromised', 'red')]):
                    comp_key = gmm_2d_stats['components'][comp_name]['idx']
                    
                    # Transform mean and covariance back to original space
                    mean_scaled = gmm_model.means_[comp_key]
                    cov_scaled = gmm_model.covariances_[comp_key]
                    
                    # Inverse transform to get original space parameters
                    mean_original = scaler.inverse_transform(mean_scaled.reshape(1, -1))[0]
                    
                    # For covariance, we need to account for the scaling
                    scale = scaler.scale_
                    cov_original = cov_scaled * np.outer(scale, scale)
                    
                    # Draw 1-sigma and 2-sigma ellipses
                    draw_ellipse(ax3, mean_original, cov_original, comp_color, 
                               f'{comp_name.capitalize()} (2σ)', n_std=2.0)
                    draw_ellipse(ax3, mean_original, cov_original, comp_color, 
                               None, n_std=1.0)
                    
                    # Plot component means
                    ax3.plot(mean_original[0], mean_original[1], 'o', 
                            color=comp_color, markersize=10, markeredgecolor='black',
                            markeredgewidth=1, label=f'{comp_name.capitalize()} mean')
                
                ax3.set_xlabel('Gene-UMI Residuals', fontsize=10)
                ax3.set_ylabel('% Mitochondrial', fontsize=10)
                ax3.set_title('2D GMM Features with Contours\n(Used for clustering)', fontsize=11, fontweight='bold')
                ax3.legend(fontsize=8, loc='best')
                plt.colorbar(scatter3, ax=ax3, label='P(Compromised)')
                
                # Add text annotation with statistics
                n_filtered_75 = np.sum(prob_compromised > 0.75)
                pct_filtered_75 = 100 * n_filtered_75 / len(prob_compromised)
                n_filtered_50 = np.sum(prob_compromised > 0.50)
                pct_filtered_50 = 100 * n_filtered_50 / len(prob_compromised)
                
                # Add statistics as text box on the figure
                textstr = f'Filtering Statistics:\n'
                textstr += f'P>0.50: {n_filtered_50:,} cells ({pct_filtered_50:.1f}%)\n'
                textstr += f'P>0.75: {n_filtered_75:,} cells ({pct_filtered_75:.1f}%) [DEFAULT]\n'
                textstr += f'\nCluster Parameters:\n'
                textstr += f'Compromised: {gmm_2d_stats["components"]["compromised"]["mito_mean"]:.1f}% mito, {gmm_2d_stats["components"]["compromised"]["residual_mean"]:.2f} resid\n'
                textstr += f'Healthy: {gmm_2d_stats["components"]["healthy"]["mito_mean"]:.1f}% mito, {gmm_2d_stats["components"]["healthy"]["residual_mean"]:.2f} resid\n'
                textstr += f'\nGMM Weights: Comp={gmm_2d_stats["components"]["compromised"]["weight"]:.2f}, Healthy={gmm_2d_stats["components"]["healthy"]["weight"]:.2f}'
                
                # Place text box in the right plot
                ax3.text(0.02, 0.98, textstr, transform=ax3.transAxes,
                        fontsize=8, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))
                
                # Overall title
                sample_id = adata.obs.get("sample_id", ["unknown"])[0]
                convergence_text = "Converged" if gmm_2d_stats['converged'] else "Not Converged (Warning)"
                plt.suptitle(f'GMM Quality Control Debug - Sample: {sample_id} - {convergence_text}\n'
                            f'Total cells: {len(prob_compromised):,} | '
                            f'Color: Red = Compromised, Blue = Healthy', 
                            fontsize=12, fontweight='bold')
                
                plt.tight_layout()
                plt.savefig(plot_path, dpi=150, bbox_inches='tight')
                plt.close()
                
                print(f"  Combined debug plot saved to: {plot_path}")
            
    else:
        raise ValueError("Missing required data (total_counts or n_genes) for 2D GMM calculation")
    
    # Save to TSV
    cell_lists.to_csv(output_path, sep='\t', index=False)
    
    # Print summary
    print("\nQC summary:")
    prob_comp = cell_lists['gmm_posterior_prob_compromised']
    convergence_status = "converged" if gmm_2d_stats['converged'] else "did not converge (using best fit)"
    print(f"  GMM {convergence_status}")
    print(f"  Cells with P(compromised) > 0.50: {(prob_comp > 0.50).sum():,}/{len(prob_comp):,} ({100*(prob_comp > 0.50).mean():.1f}%)")
    print(f"  Cells with P(compromised) > 0.75: {(prob_comp > 0.75).sum():,}/{len(prob_comp):,} ({100*(prob_comp > 0.75).mean():.1f}%)")
    print(f"  Mean P(compromised): {prob_comp.mean():.3f}")
    
    print(f"\nSaved cell lists to: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Generate QC cell filtering decisions')
    parser.add_argument('--h5ad', required=True, help='Annotated h5ad file')
    parser.add_argument('--cell-barcodes', required=True, help='File with called cell barcodes')
    parser.add_argument('--output', required=True, help='Output TSV file for cell lists')
    parser.add_argument('--plot-output', help='Output path for GMM debug plot')
    parser.add_argument('--config', help='Config YAML file')
    
    args = parser.parse_args()
    
    # Load config if provided
    config = None
    if args.config:
        with open(args.config) as f:
            config = yaml.safe_load(f)
    
    # Load called cell barcodes
    print(f"Loading called cell barcodes: {args.cell_barcodes}")
    with open(args.cell_barcodes) as f:
        called_cells = set(line.strip() for line in f)
    print(f"Loaded {len(called_cells):,} called cell barcodes")
    
    # Load h5ad file
    print(f"Loading h5ad file: {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    print(f"Loaded {adata.n_obs:,} total barcodes x {adata.n_vars:,} genes")
    
    # Filter to called cells only
    adata = adata[adata.obs_names.isin(called_cells)]
    print(f"Filtered to {adata.n_obs:,} called cells")
    
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
    save_qc_cell_lists(adata, args.output, config, save_plots=True, plot_output=args.plot_output)


if __name__ == "__main__":
    main()