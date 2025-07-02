#!/usr/bin/env Rscript
#
# EmptyDrops cell calling using R's DropletUtils
# Works with kallisto count matrix outputs
#

suppressPackageStartupMessages({
  library(DropletUtils)
  library(Matrix)
  library(methods)
  library(argparse)
  library(BiocParallel)
})

# Parse command line arguments
parser <- ArgumentParser(description='EmptyDrops cell calling using R')
parser$add_argument('--kb_dir', required=TRUE, help='Directory containing kallisto-bustools output')
parser$add_argument('--sample_id', required=TRUE, help='Sample ID')
parser$add_argument('--output_dir', required=TRUE, help='Output directory')
parser$add_argument('--emptydrops_lower', type='integer', default=100, help='emptyDrops lower threshold')
parser$add_argument('--test_ambient', action='store_true', default=FALSE, help='Test ambient genes')
parser$add_argument('--ncores', type='integer', default=4, help='Number of cores for parallel processing')

args <- parser$parse_args()

# Set up parallel processing
if (args$ncores > 1) {
  BPPARAM <- MulticoreParam(workers = args$ncores)
} else {
  BPPARAM <- SerialParam()
}

cat("Using", args$ncores, "cores for parallel processing\n")

# Create output directory
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading count matrix from:", args$kb_dir, "\n")

# Load kallisto count matrix directly
# Note: Using counts_filtered which now contains unfiltered data with annotations
counts_dir <- file.path(args$kb_dir, "counts_filtered")
matrix_file <- file.path(counts_dir, "cells_x_genes.total.mtx")
barcodes_file <- file.path(counts_dir, "cells_x_genes.barcodes.txt")
genes_file <- file.path(counts_dir, "cells_x_genes.genes.txt")

# Read matrix and metadata
count_matrix <- readMM(matrix_file)
barcodes <- readLines(barcodes_file)
genes <- read.table(genes_file, header = FALSE, stringsAsFactors = FALSE)

cat("Matrix dimensions before naming:", nrow(count_matrix), "x", ncol(count_matrix), "\n")
cat("Number of barcodes:", length(barcodes), "\n")
cat("Number of genes:", nrow(genes), "\n")

# Matrix is already in genes x cells format (transposed in Python)
cat("Matrix class:", class(count_matrix), "\n")

# Set row and column names (genes=rows, cells=columns)  
rownames(count_matrix) <- genes[,1]  # Gene IDs
colnames(count_matrix) <- barcodes
rm(genes)  # Free memory
gc()

cat("Loaded matrix:", nrow(count_matrix), "genes x", ncol(count_matrix), "barcodes\n")

# Calculate total UMI counts per barcode
total_counts <- colSums(count_matrix)
cat("UMI count range:", min(total_counts), "to", max(total_counts), "\n")

# Run barcodeRanks
cat("Running barcodeRanks analysis...\n")
br_out <- barcodeRanks(count_matrix)

# Extract knee and inflection points
knee_threshold <- metadata(br_out)$knee
inflection_threshold <- metadata(br_out)$inflection

cat("BarcodeRanks knee threshold:", knee_threshold, "\n")
cat("BarcodeRanks inflection threshold:", inflection_threshold, "\n")

# Call cells using knee and inflection
is_cell_knee <- total_counts >= knee_threshold
is_cell_inflection <- total_counts >= inflection_threshold
n_cells_knee <- sum(is_cell_knee)
n_cells_inflection <- sum(is_cell_inflection)

cat("BarcodeRanks knee:", n_cells_knee, "cells\n")
cat("BarcodeRanks inflection:", n_cells_inflection, "cells\n")

# Run emptyDrops
cat("\nRunning emptyDrops analysis...\n")
cat("  Lower threshold:", args$emptydrops_lower, "\n")
cat("  Testing multiple FDR thresholds: 0.001, 0.01, 0.05\n")

set.seed(100)
e_out <- emptyDrops(count_matrix, lower = args$emptydrops_lower, test.ambient = args$test_ambient, BPPARAM = BPPARAM)

# Clean up large objects after EmptyDrops
rm(count_matrix)
gc()

# Test multiple FDR thresholds
fdr_thresholds <- c(0.001, 0.01, 0.05)
total_umis <- sum(total_counts)

# Initialize results storage
fdr_results <- list()

for (fdr_thresh in fdr_thresholds) {
  # Call cells based on FDR
  is_cell <- e_out$FDR <= fdr_thresh & !is.na(e_out$FDR)
  n_cells <- sum(is_cell, na.rm = TRUE)
  
  # Get threshold (minimum total count among called cells)
  threshold <- NA
  if (n_cells > 0) {
    threshold <- min(total_counts[is_cell], na.rm = TRUE)
  }
  
  # Calculate percentage of UMIs in cells
  umis_in_cells_pct <- sum(total_counts[is_cell], na.rm = TRUE) / total_umis * 100
  
  # Store results
  fdr_results[[as.character(fdr_thresh)]] <- list(
    is_cell = is_cell,
    n_cells = n_cells,
    threshold = threshold,
    umis_in_cells_pct = umis_in_cells_pct
  )
  
  cat("FDR", fdr_thresh, ":", n_cells, "cells (threshold:", threshold, ", UMIs in cells:", round(umis_in_cells_pct, 1), "%)\n")
}


# Save emptyDrops results with FDR columns
# barcodes already loaded above

emptydrops_results <- data.frame(
  barcode = barcodes,
  total_counts = total_counts,
  rank = rank(-total_counts, ties.method = "first"),  # Add barcode rank
  log_prob = e_out$LogProb,
  p_value = e_out$PValue,
  fdr = e_out$FDR,
  is_cell_fdr_0.001 = ifelse(is.na(fdr_results[["0.001"]]$is_cell), FALSE, fdr_results[["0.001"]]$is_cell),
  is_cell_fdr_0.01 = ifelse(is.na(fdr_results[["0.01"]]$is_cell), FALSE, fdr_results[["0.01"]]$is_cell),
  is_cell_fdr_0.05 = ifelse(is.na(fdr_results[["0.05"]]$is_cell), FALSE, fdr_results[["0.05"]]$is_cell),
  is_cell_knee = is_cell_knee,
  is_cell_inflection = is_cell_inflection,
  limited = e_out$Limited
)

write.table(emptydrops_results,
           file = file.path(args$output_dir, "emptydrops_results.tsv"),
           sep = "\t", row.names = FALSE, quote = FALSE)

# Save cell barcodes for each FDR threshold
for (fdr_thresh in fdr_thresholds) {
  fdr_str <- gsub("\\.", "", as.character(fdr_thresh))  # Remove decimal point for filename
  called_cells <- barcodes[fdr_results[[as.character(fdr_thresh)]]$is_cell & !is.na(fdr_results[[as.character(fdr_thresh)]]$is_cell)]
  writeLines(called_cells,
            file.path(args$output_dir, paste0(args$sample_id, "_EmptyDrops_FDR", fdr_str, "_cell_barcodes.txt")))
}

# Save barcodeRanks cell barcodes
knee_cells <- barcodes[is_cell_knee]
inflection_cells <- barcodes[is_cell_inflection]
writeLines(knee_cells,
          file.path(args$output_dir, paste0(args$sample_id, "_BarcodeRanks_Knee_cell_barcodes.txt")))
writeLines(inflection_cells,
          file.path(args$output_dir, paste0(args$sample_id, "_BarcodeRanks_Inflection_cell_barcodes.txt")))

# Save summary statistics for all methods
summary_data <- list()

# Add EmptyDrops FDR results
for (i in seq_along(fdr_thresholds)) {
  fdr_thresh <- fdr_thresholds[i]
  result <- fdr_results[[as.character(fdr_thresh)]]
  
  summary_data[[length(summary_data) + 1]] <- data.frame(
    sample_id = args$sample_id,
    method = paste0("EmptyDrops_FDR", gsub("\\.", "", as.character(fdr_thresh))),
    fdr_threshold = fdr_thresh,
    n_cells_called = result$n_cells,
    threshold_used = result$threshold,
    umis_in_cells_pct = result$umis_in_cells_pct,
    emptydrops_lower_param = args$emptydrops_lower,
    n_barcodes_tested = sum(!is.na(e_out$PValue)),
    n_barcodes_limited = sum(e_out$Limited, na.rm = TRUE)
  )
}

# Add barcodeRanks results
umis_in_cells_pct_knee <- sum(total_counts[is_cell_knee], na.rm = TRUE) / total_umis * 100
umis_in_cells_pct_inflection <- sum(total_counts[is_cell_inflection], na.rm = TRUE) / total_umis * 100

summary_data[[length(summary_data) + 1]] <- data.frame(
  sample_id = args$sample_id,
  method = "BarcodeRanks_Knee",
  fdr_threshold = NA,
  n_cells_called = n_cells_knee,
  threshold_used = knee_threshold,
  umis_in_cells_pct = umis_in_cells_pct_knee,
  emptydrops_lower_param = args$emptydrops_lower,
  n_barcodes_tested = length(total_counts),
  n_barcodes_limited = 0
)

summary_data[[length(summary_data) + 1]] <- data.frame(
  sample_id = args$sample_id,
  method = "BarcodeRanks_Inflection",
  fdr_threshold = NA,
  n_cells_called = n_cells_inflection,
  threshold_used = inflection_threshold,
  umis_in_cells_pct = umis_in_cells_pct_inflection,
  emptydrops_lower_param = args$emptydrops_lower,
  n_barcodes_tested = length(total_counts),
  n_barcodes_limited = 0
)

summary_stats <- do.call(rbind, summary_data)
write.table(summary_stats,
           file = file.path(args$output_dir, "emptydrops_summary.tsv"),
           sep = "\t", row.names = FALSE, quote = FALSE)

cat("DropletUtils analysis completed!\n")
cat("Results saved to:", args$output_dir, "\n")
cat("Cell barcodes saved for:\n")
cat("  - EmptyDrops FDR thresholds: 0.001, 0.01, 0.05\n")
cat("  - BarcodeRanks knee and inflection points\n")