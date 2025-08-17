#!/usr/bin/env Rscript
#
# Minimal BarcodeRanks implementation
# Uses standalone barcodeRanks function without full DropletUtils dependency
#

suppressPackageStartupMessages({
  library(Matrix)
  library(methods)
  library(argparse)
})

# Source the standalone barcodeRanks implementation
script_dir <- dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])
script_dir <- sub("--file=", "", script_dir)
source(file.path(script_dir, "barcodeRanks_standalone.R"))

# Parse command line arguments
parser <- ArgumentParser(description='Minimal BarcodeRanks for knee/inflection detection')
parser$add_argument('--kb_dir', required=TRUE, help='Directory containing kallisto-bustools output')
parser$add_argument('--sample_id', required=TRUE, help='Sample ID')
parser$add_argument('--output_dir', required=TRUE, help='Output directory')
parser$add_argument('--expected_cells', type='integer', default=0, help='Expected number of cells (for adaptive barcodeRanks)')
parser$add_argument('--lower', type='integer', default=100, help='Lower threshold for total UMI counts')

args <- parser$parse_args()

# Create output directory
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading count matrix from:", args$kb_dir, "\n")

# Load kallisto count matrix directly
counts_dir <- file.path(args$kb_dir, "counts_filtered")
matrix_file <- file.path(counts_dir, "cells_x_genes.total.mtx")
barcodes_file <- file.path(counts_dir, "cells_x_genes.barcodes.txt")
genes_file <- file.path(counts_dir, "cells_x_genes.genes.txt")

# Check if files exist
if (!file.exists(matrix_file) || !file.exists(barcodes_file) || !file.exists(genes_file)) {
  stop("Count matrix files not found in ", counts_dir)
}

# Read the matrix (rows are genes, columns are barcodes)
count_matrix <- readMM(matrix_file)
barcodes <- readLines(barcodes_file)
genes <- readLines(genes_file)

# Set dimension names
colnames(count_matrix) <- barcodes
rownames(count_matrix) <- genes

cat("Matrix dimensions: ", nrow(count_matrix), " genes x ", ncol(count_matrix), " barcodes\n")

# Calculate total counts per barcode
total_counts <- Matrix::colSums(count_matrix)

cat("\nRunning barcodeRanks analysis...\n")

# Adaptive exclusion based on expected cells
if (args$expected_cells > 0) {
  # Exclude 10% of expected cells (minimum 50)
  exclude_n <- max(50, ceiling(args$expected_cells * 0.1))
  cat("  Using adaptive exclusion: excluding top", exclude_n, "barcodes (10% of", args$expected_cells, "expected cells)\n")
  br_out <- barcodeRanks(count_matrix, exclude.from=exclude_n, lower=args$lower)
} else {
  # expected_cells is 0 - run without exclude.from
  cat("  No expected_cells provided, using default barcodeRanks parameters\n")
  br_out <- barcodeRanks(count_matrix, lower=args$lower)
}

# Extract knee and inflection points from attributes (not S4 metadata)
knee_threshold <- attr(br_out, "knee")
inflection_threshold <- attr(br_out, "inflection")

cat("BarcodeRanks knee threshold:", knee_threshold, "\n")
cat("BarcodeRanks inflection threshold:", inflection_threshold, "\n")

# Call cells using knee and inflection
is_cell_knee <- total_counts >= knee_threshold
is_cell_inflection <- total_counts >= inflection_threshold
n_cells_knee <- sum(is_cell_knee)
n_cells_inflection <- sum(is_cell_inflection)

cat("BarcodeRanks knee:", n_cells_knee, "cells\n")
cat("BarcodeRanks inflection:", n_cells_inflection, "cells\n")

# Save results
results <- data.frame(
  barcode = barcodes,
  total_counts = total_counts,
  rank = br_out$rank,
  knee_threshold = knee_threshold,
  inflection_threshold = inflection_threshold,
  is_cell_knee = is_cell_knee,
  is_cell_inflection = is_cell_inflection,
  stringsAsFactors = FALSE
)

output_file <- file.path(args$output_dir, paste0(args$sample_id, "_barcoderanks.tsv"))
write.table(results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nResults saved to:", output_file, "\n")

# Save cell barcodes
knee_cells_file <- file.path(args$output_dir, paste0(args$sample_id, "_cells_knee.txt"))
writeLines(barcodes[is_cell_knee], knee_cells_file)
cat("Knee cells saved to:", knee_cells_file, "\n")

inflection_cells_file <- file.path(args$output_dir, paste0(args$sample_id, "_cells_inflection.txt"))
writeLines(barcodes[is_cell_inflection], inflection_cells_file)
cat("Inflection cells saved to:", inflection_cells_file, "\n")

cat("\nDone!\n")