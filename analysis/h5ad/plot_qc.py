#!/usr/bin/env python3

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

# Load data
adata = sc.read_h5ad('sample1_gex_filtered.h5ad')

# Filter to top 1000 cells by total UMI count
adata.obs['total_counts'] = adata.X.sum(axis=1).A1
adata = adata[adata.obs['total_counts'].nlargest(1000).index].copy()

# Create scatter plot
fig, ax = plt.subplots(figsize=(8, 6))
scatter = ax.scatter(adata.obs['n_genes'], adata.obs['total_counts'], 
                    c=adata.obs['pct_counts_mt'], cmap='viridis', alpha=0.7)
ax.set_xlabel('Number of Genes')
ax.set_ylabel('Total UMI Counts')
ax.set_title('Top 1000 Cells: UMIs vs Genes (colored by Mito %)')
plt.colorbar(scatter, label='Mitochondrial %')
plt.tight_layout()
plt.savefig('qc_scatter.png', dpi=150, bbox_inches='tight')
plt.show()

print(f"Plotted {adata.n_obs} cells")
print(f"UMI range: {adata.obs['total_counts'].min():.0f} - {adata.obs['total_counts'].max():.0f}")
print(f"Gene range: {adata.obs['n_genes'].min()} - {adata.obs['n_genes'].max()}")
print(f"Mito % range: {adata.obs['pct_counts_mt'].min():.1f}% - {adata.obs['pct_counts_mt'].max():.1f}%")