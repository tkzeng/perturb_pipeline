#!/usr/bin/env Rscript

# Differential expression analysis with edgeR including visualization
# Takes pseudobulk outputs and performs DE analysis with plots

library(edgeR)
library(readr)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(argparse)

# Parse command line arguments
parser <- ArgumentParser(description='Differential expression analysis with edgeR including visualization')
parser$add_argument('--counts', required=TRUE, help='Path to raw counts TSV file')
parser$add_argument('--norm-factors', required=TRUE, help='Path to normalization factors TSV file')
parser$add_argument('--metadata', required=TRUE, help='Path to group metadata TSV file')
parser$add_argument('--output-dir', required=TRUE, help='Output directory for results')
parser$add_argument('--contrasts', required=TRUE, help='Comma-separated list of contrasts (e.g., treatment_vs_control,high_vs_low)')
parser$add_argument('--fdr-threshold', type='double', required=TRUE, help='FDR threshold for significance')
parser$add_argument('--logfc-threshold', type='double', required=TRUE, help='Log fold change threshold for volcano plots')
parser$add_argument('--group-column', required=TRUE, help='Column name in metadata that defines groups')
parser$add_argument('--design-formula', required=TRUE, help='Design formula (e.g., ~group)')
parser$add_argument('--plot-dir', required=TRUE, help='Base directory for dashboard plots')
parser$add_argument('--source', required=TRUE, help='Data source identifier')
parser$add_argument('--processing', required=TRUE, help='Processing type identifier')

args <- parser$parse_args()

# No need to validate - argparse handles required arguments

cat("=== edgeR Differential Expression Analysis ===\n")
cat("Counts file:", args$counts, "\n")
cat("Norm factors file:", args$norm_factors, "\n")
cat("Metadata file:", args$metadata, "\n")
cat("Output directory:", args$output_dir, "\n")
cat("Contrasts:", args$contrasts, "\n")
cat("FDR threshold:", args$fdr_threshold, "\n")
cat("LogFC threshold:", args$logfc_threshold, "\n")
cat("Group column:", args$group_column, "\n")
cat("Design formula:", args$design_formula, "\n")

# Create output directories
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(args$output_dir, "data"), showWarnings = FALSE)

# Load data
cat("\nLoading data...\n")
counts_df <- read_tsv(args$counts, show_col_types = FALSE)
norm_factors_df <- read_tsv(args$norm_factors, show_col_types = FALSE)
metadata_df <- read_tsv(args$metadata, show_col_types = FALSE)

# Prepare count matrix
gene_names <- counts_df$gene
gene_ids <- counts_df$gene_id
counts_matrix <- as.matrix(counts_df %>% select(-gene, -gene_id))
rownames(counts_matrix) <- gene_ids

# Align sample order
sample_order <- match(colnames(counts_matrix), norm_factors_df$sample_id)
if (any(is.na(sample_order))) {
  stop("Sample mismatch between counts and normalization factors")
}

# Create DGEList with existing normalization factors
cat("Creating DGEList with TMM normalization factors...\n")
dge <- DGEList(counts = counts_matrix)
dge$samples$lib.size <- norm_factors_df$lib_size[sample_order]
dge$samples$norm.factors <- norm_factors_df$norm_factors[sample_order]

# Add metadata to DGEList
metadata_aligned <- metadata_df[match(colnames(counts_matrix), metadata_df$group_id), ]
if (any(is.na(metadata_aligned$group_id))) {
  stop("Metadata missing for some samples")
}

# Check that specified group column exists
if (!args$group_column %in% colnames(metadata_aligned)) {
  stop(paste("Group column", args$group_column, "not found in metadata. Available columns:",
             paste(colnames(metadata_aligned), collapse=", ")))
}

dge$samples$group <- factor(metadata_aligned[[args$group_column]])

# Design matrix
cat("\nCreating design matrix...\n")
if (args$design_formula == "~group") {
  design <- model.matrix(~0 + group, data = dge$samples)
  colnames(design) <- levels(dge$samples$group)
} else {
  # Parse custom formula if provided
  formula_obj <- as.formula(args$design_formula)
  design <- model.matrix(formula_obj, data = metadata_aligned)
}

cat("Design matrix dimensions:", nrow(design), "x", ncol(design), "\n")
cat("Groups:", paste(colnames(design), collapse=", "), "\n")

# Estimate dispersions
cat("\nEstimating dispersions...\n")
dge <- estimateDisp(dge, design, robust = TRUE)

# Fit GLM
cat("Fitting GLM...\n")
fit <- glmQLFit(dge, design, robust = TRUE)

# Parse contrasts
contrasts_list <- strsplit(args$contrasts, ",")[[1]]

for (contrast_str in contrasts_list) {
  cat("\n=== Processing contrast:", contrast_str, "===\n")
  
  # Parse contrast (format: groupA_vs_groupB)
  contrast_parts <- strsplit(contrast_str, "_vs_")[[1]]
  if (length(contrast_parts) != 2) {
    cat("Warning: Invalid contrast format:", contrast_str, "(expected: groupA_vs_groupB). Skipping.\n")
    next
  }
  
  group1 <- contrast_parts[1]
  group2 <- contrast_parts[2]
  
  # Check if groups exist in design
  if (!group1 %in% colnames(design) || !group2 %in% colnames(design)) {
    cat("Warning: Groups not found in design matrix. Skipping contrast:", contrast_str, "\n")
    cat("Available groups:", paste(colnames(design), collapse=", "), "\n")
    next
  }
  
  # Create contrast
  contrast_vector <- makeContrasts(contrasts = paste0(group1, "-", group2), levels = design)
  
  # Perform differential expression
  qlf <- glmQLFTest(fit, contrast = contrast_vector)
  
  # Get results
  de_results <- topTags(qlf, n = Inf, adjust.method = "BH")$table
  
  # Add gene names
  de_results$gene <- gene_names[match(rownames(de_results), gene_ids)]
  de_results$gene_id <- rownames(de_results)
  
  # Reorder columns
  de_results <- de_results %>%
    select(gene_id, gene, logFC, logCPM, F, PValue, FDR)
  
  # Save full results
  output_file <- file.path(args$output_dir, "data", paste0("de_results_", contrast_str, ".tsv"))
  write_tsv(de_results, output_file)
  cat("Saved results to:", output_file, "\n")
  
  # Count significant genes for summary
  de_significant <- de_results %>%
    filter(FDR < args$fdr_threshold)
  cat("Significant genes (FDR <", args$fdr_threshold, "):", nrow(de_significant), "\n")
  
  
  # 1. Volcano Plot - using dashboard directory structure
  cat("Generating volcano plot...\n")
  # Structure: {plot_dir}/{metric_name}/{source}_{processing}/{contrast}/plot.png
  volcano_dir <- file.path(args$plot_dir, "de_volcano", paste0(args$source, "_", args$processing), 
                          contrast_str)
  dir.create(volcano_dir, showWarnings = FALSE, recursive = TRUE)
  
  de_results$significant <- de_results$FDR < args$fdr_threshold & abs(de_results$logFC) > args$logfc_threshold
  de_results$neg_log10_fdr <- -log10(de_results$FDR + 1e-300)  # Add small value to avoid inf
  
  # Label top genes - 20 from each direction
  top_up <- de_results %>%
    filter(significant & logFC > 0) %>%
    arrange(FDR) %>%
    head(20)
  
  top_down <- de_results %>%
    filter(significant & logFC < 0) %>%
    arrange(FDR) %>%
    head(20)
  
  top_genes <- rbind(top_up, top_down)
  
  p_volcano <- ggplot(de_results, aes(x = logFC, y = neg_log10_fdr)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red")) +
    geom_vline(xintercept = c(-args$logfc_threshold, args$logfc_threshold), 
               linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(args$fdr_threshold), 
               linetype = "dashed", alpha = 0.5) +
    geom_text_repel(data = top_genes, 
                    aes(label = gene),
                    size = 3,
                    max.overlaps = 20) +
    labs(x = "Log2 Fold Change",
         y = "-log10(FDR)",
         title = paste("Volcano Plot:", contrast_str),
         subtitle = paste(nrow(de_significant), "significant genes")) +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(file.path(volcano_dir, "plot.png"), p_volcano, width = 8, height = 6, dpi = 150)
  
  # 2. MA Plot - using dashboard directory structure
  cat("Generating MA plot...\n")
  ma_dir <- file.path(args$plot_dir, "de_ma_plot", paste0(args$source, "_", args$processing), 
                     contrast_str)
  dir.create(ma_dir, showWarnings = FALSE, recursive = TRUE)
  
  p_ma <- ggplot(de_results, aes(x = logCPM, y = logFC)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red")) +
    geom_hline(yintercept = 0, linetype = "solid", alpha = 0.5) +
    geom_hline(yintercept = c(-args$logfc_threshold, args$logfc_threshold), 
               linetype = "dashed", alpha = 0.5) +
    geom_text_repel(data = top_genes, 
                    aes(label = gene),
                    size = 3,
                    max.overlaps = 20) +
    labs(x = "Average log2 CPM",
         y = "Log2 Fold Change",
         title = paste("MA Plot:", contrast_str),
         subtitle = paste(nrow(de_significant), "significant genes")) +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(file.path(ma_dir, "plot.png"), p_ma, width = 8, height = 6, dpi = 150)
  
}

cat("\n=== Differential expression analysis complete ===\n")
cat("Results saved to:", args$output_dir, "\n")