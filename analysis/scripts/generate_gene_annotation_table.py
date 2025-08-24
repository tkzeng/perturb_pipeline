#!/usr/bin/env python3
"""
Generate Gene Annotation Table
===============================

This script pre-processes gene annotations from various sources into a single 
TSV table for fast lookups. This replaces slow gffutils database queries
with vectorized pandas operations.

Usage:
    python generate_gene_annotation_table.py \
        --gene-database /path/to/gencode.db \
        --ribosomal-genes /path/to/ribosomal_genes.tsv \
        --cell-cycle-genes /path/to/cell_cycle_genes.txt \
        --output /path/to/gene_annotations.tsv
"""

import argparse
import pandas as pd
import gffutils
import sys
from pathlib import Path

def log_print(*args):
    """Simple logging function."""
    message = " ".join(str(arg) for arg in args)
    print(message, file=sys.stderr, flush=True)


def load_ribosomal_genes(ribosomal_genes_path):
    """Load ribosomal gene HGNC IDs from TSV file."""
    log_print("📋 Loading ribosomal genes...")
    ribosomal_hgnc_ids = set()
    
    with open(ribosomal_genes_path) as f:
        header = f.readline().strip().split('\t')
        hgnc_col = header.index("HGNC ID")
        
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) > hgnc_col:
                ribosomal_hgnc_ids.add(cols[hgnc_col])
    
    log_print(f"  Loaded {len(ribosomal_hgnc_ids)} ribosomal HGNC IDs")
    return ribosomal_hgnc_ids


def load_cell_cycle_genes(cell_cycle_genes_path):
    """Load cell cycle genes and phase information."""
    log_print("📋 Loading cell cycle genes...")
    
    with open(cell_cycle_genes_path) as f:
        cell_cycle_genes_list = [x.strip() for x in f]
    
    cell_cycle_genes = set(cell_cycle_genes_list)
    s_phase_genes = set(cell_cycle_genes_list[:43])  # First 43 are S phase
    g2m_phase_genes = set(cell_cycle_genes_list[43:])  # Rest are G2M
    
    log_print(f"  Loaded {len(cell_cycle_genes)} cell cycle genes")
    log_print(f"    S phase: {len(s_phase_genes)} genes")
    log_print(f"    G2M phase: {len(g2m_phase_genes)} genes")
    
    return cell_cycle_genes, s_phase_genes, g2m_phase_genes


def extract_gene_annotations_from_db(gene_database_path, ribosomal_hgnc_ids, 
                                    cell_cycle_genes, s_phase_genes):
    """Extract all gene annotations from gffutils database into a DataFrame."""
    log_print("🔍 Extracting gene annotations from database...")
    
    # Open database
    db = gffutils.FeatureDB(gene_database_path)
    
    # Lists to store annotation data
    gene_data = []
    
    # Counter for progress
    count = 0
    
    # Iterate through all genes in database
    for gene in db.features_of_type('gene'):
        gene_id = gene.id
        
        # Extract attributes
        gene_type = gene.attributes.get('gene_type', [''])[0]
        hgnc_id = gene.attributes.get('hgnc_id', [''])[0]
        gene_name = gene.attributes.get('gene_name', [''])[0]
        chromosome = gene.seqid
        
        # Check if ribosomal
        is_ribosomal = hgnc_id in ribosomal_hgnc_ids
        
        # Check cell cycle status
        is_cell_cycle = gene_name in cell_cycle_genes
        cell_cycle_phase = ''
        if is_cell_cycle:
            if gene_name in s_phase_genes:
                cell_cycle_phase = 'S'
            else:
                cell_cycle_phase = 'G2M'
        
        # Store all annotations
        gene_data.append({
            'gene_id': gene_id,
            'gene_type': gene_type,
            'hgnc_id': hgnc_id,
            'gene_name': gene_name,
            'chromosome': chromosome,
            'ribosomal': is_ribosomal,
            'cell_cycle': is_cell_cycle,
            'cell_cycle_phase': cell_cycle_phase
        })
        
        count += 1
        if count % 5000 == 0:
            log_print(f"  Processed {count} genes...")
    
    log_print(f"✅ Extracted annotations for {count} genes")
    
    # Create DataFrame
    df = pd.DataFrame(gene_data)
    
    # Set gene_id as index for fast lookups
    df.set_index('gene_id', inplace=True)
    
    return df


def main():
    parser = argparse.ArgumentParser(description='Generate gene annotation table')
    parser.add_argument('--gene-database', required=True,
                       help='Path to gffutils gene database')
    parser.add_argument('--ribosomal-genes', required=True,
                       help='Path to ribosomal genes TSV file')
    parser.add_argument('--cell-cycle-genes', required=True,
                       help='Path to cell cycle genes text file')
    parser.add_argument('--output', required=True,
                       help='Output path for gene annotation table (TSV format)')
    
    args = parser.parse_args()
    
    log_print("=" * 80)
    log_print("🚀 GENERATING GENE ANNOTATION TABLE")
    log_print("=" * 80)
    
    # Load reference data
    ribosomal_hgnc_ids = load_ribosomal_genes(args.ribosomal_genes)
    cell_cycle_genes, s_phase_genes, g2m_phase_genes = load_cell_cycle_genes(args.cell_cycle_genes)
    
    # Extract all gene annotations
    gene_df = extract_gene_annotations_from_db(
        args.gene_database,
        ribosomal_hgnc_ids,
        cell_cycle_genes,
        s_phase_genes
    )
    
    # Show summary statistics
    log_print("\n📊 Annotation summary:")
    log_print(f"  Total genes: {len(gene_df)}")
    log_print(f"  Ribosomal genes: {gene_df['ribosomal'].sum()}")
    log_print(f"  Cell cycle genes: {gene_df['cell_cycle'].sum()}")
    log_print(f"    S phase: {(gene_df['cell_cycle_phase'] == 'S').sum()}")
    log_print(f"    G2M phase: {(gene_df['cell_cycle_phase'] == 'G2M').sum()}")
    
    # Gene type distribution
    log_print("\n📊 Gene type distribution (top 10):")
    gene_type_counts = gene_df['gene_type'].value_counts().head(10)
    for gene_type, count in gene_type_counts.items():
        log_print(f"    {gene_type}: {count}")
    
    # Save as TSV for both human readability and fast loading
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Change extension to .tsv if not already
    if not str(output_path).endswith('.tsv'):
        output_path = output_path.with_suffix('.tsv')
    
    log_print(f"\n💾 Saving to {output_path}...")
    gene_df.to_csv(output_path, sep='\t', index=True)
    
    log_print(f"\n✅ Gene annotation table created successfully!")
    log_print(f"   File size: {output_path.stat().st_size / 1024 / 1024:.1f} MB")
    
    # Show example of how to use it
    log_print("\n📖 Usage example:")
    log_print("   import pandas as pd")
    log_print(f"   gene_annotations = pd.read_csv('{output_path}', sep='\\t', index_col='gene_id')")
    log_print("   # Use with merge for vectorized annotation:")
    log_print("   adata.var = adata.var.merge(gene_annotations, left_index=True, right_index=True, how='left')")


if __name__ == "__main__":
    main()