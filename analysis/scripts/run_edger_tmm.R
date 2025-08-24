#!/usr/bin/env Rscript

# TMM normalization using edgeR
# Usage: Rscript run_edger_tmm.R <input_counts.tsv> <output_prefix>

library(edgeR)
library(readr)
library(tibble)
library(dplyr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  stop("Usage: Rscript run_edger_tmm.R <input_counts.tsv> <output_prefix> [reference_norm_factors.tsv]")
}

input_file <- args[1]
output_prefix <- args[2]
reference_norm_file <- if (length(args) == 3) args[3] else NULL

cat("Loading count matrix from:", input_file, "\n")
if (!is.null(reference_norm_file)) {
  cat("Using reference normalization factors from:", reference_norm_file, "\n")
}

# Read count matrix (genes × samples)
counts_df <- read_tsv(input_file, show_col_types = FALSE)
gene_names <- counts_df$gene
gene_ids <- counts_df$gene_id
counts_matrix <- as.matrix(counts_df %>% select(-gene, -gene_id))
rownames(counts_matrix) <- gene_ids  # Use gene_ids as rownames for processing

cat("Count matrix shape:", nrow(counts_matrix), "genes ×", ncol(counts_matrix), "samples\n")

# Create DGEList
dge <- DGEList(counts = counts_matrix)

if (!is.null(reference_norm_file)) {
  # Use reference normalization factors and gene filtering
  cat("Loading reference normalization factors...\n")
  ref_norm_df <- read_tsv(reference_norm_file, show_col_types = FALSE)
  ref_gene_filter_file <- gsub("_norm_factors.tsv", "_gene_filter.tsv", reference_norm_file)
  ref_gene_filter_df <- read_tsv(ref_gene_filter_file, show_col_types = FALSE)
  
  # Apply reference gene filtering
  keep <- ref_gene_filter_df$kept_by_edger[match(gene_ids, ref_gene_filter_df$gene_id)]
  keep[is.na(keep)] <- FALSE  # Genes not in reference are not kept
  dge <- dge[keep, ]
  
  # Apply reference normalization factors
  # Match sample order between reference and current data
  sample_order <- match(colnames(dge$counts), ref_norm_df$sample_id)
  if (any(is.na(sample_order))) {
    stop("Sample mismatch between reference and current data")
  }
  
  dge$samples$lib.size <- ref_norm_df$lib_size[sample_order]
  dge$samples$norm.factors <- ref_norm_df$norm_factors[sample_order]
  
  cat("Applied reference gene filtering:", sum(keep), "genes kept (from", length(keep), ")\n")
  cat("Applied reference normalization factors\n")
} else {
  # Standard processing: filter genes and calculate normalization factors
  keep <- filterByExpr(dge)
  dge <- dge[keep, ]
  
  cat("Filtered to", sum(keep), "genes (from", length(keep), ")\n")
  
  # Calculate TMM normalization factors
  dge <- calcNormFactors(dge, method = "TMM")
  cat("Calculated TMM normalization factors\n")
}

# Get TMM-normalized CPM
tmm_cpm <- cpm(dge, normalized.lib.sizes = TRUE)

# Get simple CPM (no normalization factors)
simple_cpm <- cpm(dge, normalized.lib.sizes = FALSE)

# Get raw counts for filtered genes
raw_counts <- dge$counts

# Save results
cat("Saving results with prefix:", output_prefix, "\n")

# Create mapping of filtered genes to their names
filtered_gene_ids <- rownames(tmm_cpm)
filtered_gene_names <- gene_names[match(filtered_gene_ids, gene_ids)]

# Save TMM-normalized CPM
tmm_df <- as_tibble(tmm_cpm, rownames = "gene_id")
tmm_df <- tmm_df %>% 
  mutate(gene = filtered_gene_names, .before = 1)
write_tsv(tmm_df, paste0(output_prefix, "_tmm_cpm.tsv"))

# Save simple CPM
simple_cpm_df <- as_tibble(simple_cpm, rownames = "gene_id")
simple_cpm_df <- simple_cpm_df %>% 
  mutate(gene = filtered_gene_names, .before = 1)
write_tsv(simple_cpm_df, paste0(output_prefix, "_simple_cpm.tsv"))

# Save raw counts (filtered)
raw_counts_df <- as_tibble(raw_counts, rownames = "gene_id")
raw_counts_df <- raw_counts_df %>% 
  mutate(gene = filtered_gene_names, .before = 1)
write_tsv(raw_counts_df, paste0(output_prefix, "_raw_counts.tsv"))

# Save gene filtering info
keep_df <- tibble(
  gene_id = gene_ids,
  gene = gene_names,
  kept_by_edger = keep
)
write_tsv(keep_df, paste0(output_prefix, "_gene_filter.tsv"))

# Save normalization factors and library sizes
norm_info_df <- tibble(
  sample_id = colnames(dge$counts),
  lib_size = dge$samples$lib.size,
  norm_factors = dge$samples$norm.factors,
  effective_lib_size = dge$samples$lib.size * dge$samples$norm.factors
)
write_tsv(norm_info_df, paste0(output_prefix, "_norm_factors.tsv"))

cat("edgeR TMM normalization completed successfully!\n")
cat("Output files:\n")
cat("  -", paste0(output_prefix, "_tmm_cpm.tsv"), "\n")
cat("  -", paste0(output_prefix, "_simple_cpm.tsv"), "\n") 
cat("  -", paste0(output_prefix, "_raw_counts.tsv"), "\n")
cat("  -", paste0(output_prefix, "_gene_filter.tsv"), "\n")
cat("  -", paste0(output_prefix, "_norm_factors.tsv"), "\n")