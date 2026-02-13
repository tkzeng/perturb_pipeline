#!/usr/bin/env Rscript
#
# DropletUtils cell calling (EmptyDrops and BarcodeRanks) 
# Works with kallisto count matrix outputs
#
# NOTE: EmptyDrops is DEPRECATED. Use BarcodeRanks methods instead.
# EmptyDrops code is preserved for backwards compatibility but is not actively maintained.
#

suppressPackageStartupMessages({
  library(DropletUtils)
  library(Matrix)
  library(methods)
  library(argparse)
  library(BiocParallel)
})

# Parse command line arguments
parser <- ArgumentParser(description='DropletUtils cell calling (EmptyDrops and BarcodeRanks)')
parser$add_argument('--kb_dir', required=TRUE, help='Directory containing kallisto-bustools output')
parser$add_argument('--sample_id', required=TRUE, help='Sample ID')
parser$add_argument('--output_dir', required=TRUE, help='Output directory')
parser$add_argument('--emptydrops_lower', type='integer', default=100, help='emptyDrops lower threshold')
parser$add_argument('--test_ambient', action='store_true', default=FALSE, help='Test ambient genes')
parser$add_argument('--ncores', type='integer', default=4, help='Number of cores for parallel processing')
parser$add_argument('--expected_cells', type='integer', default=0, help='Expected number of cells (for adaptive barcodeRanks)')
parser$add_argument('--run_emptydrops', action='store_true', default=FALSE, help='Run EmptyDrops analysis')
parser$add_argument('--run_barcoderanks', action='store_true', default=FALSE, help='Run BarcodeRanks analysis')
parser$add_argument('--fdr_cutoffs', type='character', help='Comma-separated FDR cutoffs for EmptyDrops')
parser$add_argument('--run_by_group', action='store_true', default=FALSE, help='If set, run BarcodeRanks per group using counts_filtered/barcode_with_group.csv')
parser$add_argument('--plot_dir', required=TRUE, help='Plot output directory')

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
plot_dir <- file.path(args$plot_dir, "dropletutills_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("Plot storage path:", plot_dir,  "\n")
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

# Initialize placeholders for results
is_cell_knee <- rep(FALSE, ncol(count_matrix))
is_cell_inflection <- rep(FALSE, ncol(count_matrix))
knee_threshold <- NA
inflection_threshold <- NA
fitted_values <- rep(NA, ncol(count_matrix))  # Store fitted spline values
br_rank <- rep(NA, ncol(count_matrix))  # Store BarcodeRanks rank values

# Run barcodeRanks if requested
if (args$run_barcoderanks) {
  cat("Running barcodeRanks analysis...\n")
  
  # Adaptive exclusion based on expected cells
  if (args$expected_cells > 0) {
    # Exclude 10% of expected cells (minimum 50)
    exclude_n <- max(50, ceiling(args$expected_cells * 0.1))
    cat("  Using adaptive exclusion: excluding top", exclude_n, "barcodes (10% of", args$expected_cells, "expected cells)\n")
    cat("  Using df=5 for spline fitting\n")
    br_out <- barcodeRanks(count_matrix, exclude.from=exclude_n, df=5)
  } else {
    # expected_cells is 0 - run without exclude.from
    cat("  No expected_cells provided, using default barcodeRanks parameters\n")
    cat("  Using df=5 for spline fitting\n")
    br_out <- barcodeRanks(count_matrix, df=5)
  }
  
  # Extract knee and inflection points
  knee_threshold <- metadata(br_out)$knee
  inflection_threshold <- metadata(br_out)$inflection
  
  # Extract fitted values and ranks for each barcode
  # Match barcodes to br_out results
  barcode_idx <- match(barcodes, rownames(br_out))
  fitted_values[!is.na(barcode_idx)] <- br_out$fitted[barcode_idx[!is.na(barcode_idx)]]
  br_rank[!is.na(barcode_idx)] <- br_out$rank[barcode_idx[!is.na(barcode_idx)]]
  
  cat("BarcodeRanks knee threshold:", knee_threshold, "\n")
  cat("BarcodeRanks inflection threshold:", inflection_threshold, "\n")
  
  # Call cells using knee and inflection
  is_cell_knee <- total_counts >= knee_threshold
  is_cell_inflection <- total_counts >= inflection_threshold
  n_cells_knee <- sum(is_cell_knee)
  n_cells_inflection <- sum(is_cell_inflection)
  
  cat("BarcodeRanks knee:", n_cells_knee, "cells\n")
  cat("BarcodeRanks inflection:", n_cells_inflection, "cells\n")

  # ---- Save global BarcodeRanks plot ----
  plot_file_all <- file.path(plot_dir, "all_samples_BarcodeRanks.png")

  png(filename = plot_file_all, width = 1800, height = 1300, res = 220)

  plot(br_out$rank, br_out$total,
      log = "xy",
      xlab = "Rank",
      ylab = "Total UMI count",
      main = paste0("BarcodeRanks (ALL): ", args$sample_id,
                    " | n=", ncol(count_matrix)),
      pch = 16, cex = 0.4)

  o <- order(br_out$rank)
  lines(br_out$rank[o], br_out$fitted[o], col="red",lwd=2)

  abline(h=knee_threshold, col="dodgerblue", lty=2, lwd=2)
  abline(h=inflection_threshold, col="forestgreen", lty=2, lwd=2)

  legend("bottomleft",
        lty = 2, lwd = 2,
        col = c("dodgerblue", "forestgreen"),
        legend = c(
          sprintf("knee: UMI>=%.0f | n_cells=%d", knee_threshold, n_cells_knee),
          sprintf("inflection: UMI>=%.0f | n_cells=%d", inflection_threshold, n_cells_inflection)
        ),
        bty = "n")

  dev.off()
} else {
  cat("Skipping barcodeRanks analysis (not requested)\n")
}

# -----------------------------
# Calling cells by Per-group BarcodeRanks (optional, runs AFTER the global barcodeRanks)
# -----------------------------

# loading barcode_with_group_list for cell calling

if (args$run_by_group) {
  barcode_with_group_file <- file.path(counts_dir, "barcode_with_group.csv")
  meta <- read.csv(barcode_with_group_file, stringsAsFactors = FALSE)
  group_col <- "cell_calling_by_group"

  if (!("barcode" %in% colnames(meta))) stop("barcode_with_group.csv missing column: barcode")
  if (!(group_col %in% colnames(meta))) stop(paste0("barcode_with_group.csv missing column: ", group_col))

  idx <- match(colnames(count_matrix), meta$barcode)
  if (any(is.na(idx))) {
    missing <- colnames(count_matrix)[is.na(idx)]
    stop(paste0("barcode_with_group.csv missing ", length(missing), " barcodes. Example: ", missing[1]))
  }
  meta <- meta[idx, , drop=FALSE]
  if (!all(meta$barcode == colnames(count_matrix))) stop("Barcode order not aligned!")

  group_ids <- meta[[group_col]]
  groups <- unique(group_ids[!is.na(group_ids)])

  cat("Cell calling split by group:", length(groups), "groups\n")
}

if (args$run_barcoderanks && args$run_by_group) {
  cat("Running BarcodeRanks by group...\n")

  group_summary_list <- list()
  group_results_rows <- list()  # per-barcode results for plotting

  for (g in groups) {
    g_safe <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(g))
    keep <- (group_ids == g)
    if (sum(keep) == 0) next

    current_counts <- count_matrix[, keep, drop=FALSE]
    current_total  <- colSums(current_counts)
    n_unique_total <- length(unique(as.numeric(current_total)))
    n_barcodes <- length(current_total)
    
    cat("Group:", g,"n_barcodes:", n_barcodes,"n_unique_total:", n_unique_total, "\n")
    
    # if too few barcodes, record zeros and continue
    if (n_barcodes < 3 || n_unique_total < 10) {
      cat("  [WARN] Skip barcodeRanks for group", g, "because n_barcodes < 3 or unique total < 10 \n")

      group_summary_list[[length(group_summary_list) + 1]] <- data.frame(
        sample_id = args$sample_id,
        group = g,
        method = "BarcodeRanks_Knee",
        n_cells_called = 0,
        threshold_used = NA,
        umis_in_cells_pct = 0,
        n_barcodes_tested = n_barcodes,
        n_unique_total = n_unique_total
      )
      group_summary_list[[length(group_summary_list) + 1]] <- data.frame(
        sample_id = args$sample_id,
        group = g,
        method = "BarcodeRanks_Inflection",
        n_cells_called = 0,
        threshold_used = NA,
        umis_in_cells_pct = 0,
        n_barcodes_tested = n_barcodes,
        n_unique_total = n_unique_total
      )

      rm(current_counts, current_total); gc()
      next
    }

    # run barcodeRanks with error capture
    tryCatch({
      # IMPORTANT: for small groups, exclude.from should be small or 0
      # Here keep simplest: no exclude.from in group mode (your current choice).
      br_g <- barcodeRanks(current_counts, df=5)

      knee_g <- metadata(br_g)$knee
      infl_g <- metadata(br_g)$inflection

      is_knee_g <- current_total >= knee_g
      is_infl_g <- current_total >= infl_g

      total_umis_g <- sum(current_total)

      # summary rows
      group_summary_list[[length(group_summary_list) + 1]] <- data.frame(
        sample_id = args$sample_id,
        group = g,
        method = "BarcodeRanks_Knee",
        n_cells_called = sum(is_knee_g),
        threshold_used = knee_g,
        umis_in_cells_pct = ifelse(total_umis_g > 0, sum(current_total[is_knee_g]) / total_umis_g * 100, 0),
        n_barcodes_tested = n_barcodes,
        n_unique_total = n_unique_total
      )
      group_summary_list[[length(group_summary_list) + 1]] <- data.frame(
        sample_id = args$sample_id,
        group = g,
        method = "BarcodeRanks_Inflection",
        n_cells_called = sum(is_infl_g),
        threshold_used = infl_g,
        umis_in_cells_pct = ifelse(total_umis_g > 0, sum(current_total[is_infl_g]) / total_umis_g * 100, 0),
        n_barcodes_tested = n_barcodes,
        n_unique_total = n_unique_total
      )

      # per-barcode table for plotting knee curves
      br_df <- as.data.frame(br_g)
      # align to current barcodes order
      idx <- match(colnames(current_counts), rownames(br_df))
      br_rank_vec <- rep(NA_real_, n_barcodes)
      br_fit_vec  <- rep(NA_real_, n_barcodes)
      if (any(!is.na(idx))) {
        br_rank_vec[!is.na(idx)] <- br_df$rank[idx[!is.na(idx)]]
        br_fit_vec[!is.na(idx)]  <- br_df$fitted[idx[!is.na(idx)]]
      }

      group_results_rows[[length(group_results_rows) + 1]] <- data.frame(
        sample_id = args$sample_id,
        group = g,
        barcode = colnames(current_counts),
        total_counts = as.numeric(current_total),
        rank = rank(-as.numeric(current_total), ties.method = "first"),
        br_rank = br_rank_vec,
        br_fitted = br_fit_vec,
        knee_threshold = knee_g,
        inflection_threshold = infl_g,
        is_cell_knee = as.logical(is_knee_g),
        is_cell_inflection = as.logical(is_infl_g),
        stringsAsFactors = FALSE
      )
      # ---- Save BarcodeRanks plot per group ----
      plot_file <- file.path(plot_dir, paste0(g_safe, "_BarcodeRanks.png"))

      png(filename = plot_file, width = 1600, height = 1200, res = 200)

      # dropletutils doc style plot: rank vs total (log-log)
      plot(br_g$rank, br_g$total,
          log = "xy",
          xlab = "Rank",
          ylab = "Total UMI count",
          main = paste0("BarcodeRanks: ", args$sample_id, " | group: ", g),
          pch = 16, cex = 0.4)

      o <- order(br_g$rank)
      lines(br_g$rank[o], br_g$fitted[o], col = "red", lwd = 2)

      abline(h = knee_g, col = "dodgerblue", lty = 2, lwd = 2)
      abline(h = infl_g, col = "forestgreen", lty = 2, lwd = 2)

      legend("bottomleft",
            lty = 2,
            lwd = 2,
            col = c("dodgerblue", "forestgreen"),
            legend = c(
              sprintf("knee: UMI>=%.0f | n_cells=%d", knee_g, sum(is_knee_g)),
              sprintf("inflection: UMI>=%.0f | n_cells=%d", infl_g, sum(is_infl_g))
            ),
            bty = "n")

      dev.off()

      rm(br_g, br_df)

    }, error = function(e) {
      cat(paste0("  [ERROR] BarcodeRanks FAILED for group ", g, ". Error: ", e$message, "\n"))
      cat("  -> Record n_cells_called=0 and continue.\n")

      group_summary_list[[length(group_summary_list) + 1]] <- data.frame(
        sample_id = args$sample_id,
        group = g,
        method = "BarcodeRanks_Knee",
        n_cells_called = 0,
        threshold_used = NA,
        umis_in_cells_pct = 0,
        n_barcodes_tested = n_barcodes,
        n_unique_total = n_unique_total
      )
      group_summary_list[[length(group_summary_list) + 1]] <- data.frame(
        sample_id = args$sample_id,
        group = g,
        method = "BarcodeRanks_Inflection",
        n_cells_called = 0,
        threshold_used = NA,
        umis_in_cells_pct = 0,
        n_barcodes_tested = n_barcodes,
        n_unique_total = n_unique_total
      )
    })

    rm(current_counts, current_total)
    gc()
  }

  # write summary_by_group
  if (length(group_summary_list) > 0) {
    group_summary <- do.call(rbind, group_summary_list)
    write.table(group_summary,
      file = file.path(args$output_dir, "cell_calling_summary_by_group.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE)
  }

  # write per-barcode results for plotting
  if (length(group_results_rows) > 0) {
    group_results <- do.call(rbind, group_results_rows)
    write.table(group_results,
      file = file.path(args$output_dir, "dropletutils_results_by_group.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE)
  }
}


# Initialize EmptyDrops results with NA
e_out <- list(
  LogProb = rep(NA, ncol(count_matrix)),
  PValue = rep(NA, ncol(count_matrix)),
  FDR = rep(NA, ncol(count_matrix)),
  Limited = rep(FALSE, ncol(count_matrix))
)

# Initialize fdr_thresholds as empty
fdr_thresholds <- numeric(0)

# Run emptyDrops if requested
if (args$run_emptydrops) {
  # Parse FDR cutoffs - REQUIRED when running emptyDrops
  if (is.null(args$fdr_cutoffs)) {
    stop("Error: --fdr_cutoffs must be provided when --run_emptydrops is set")
  }
  fdr_thresholds <- as.numeric(strsplit(args$fdr_cutoffs, ",")[[1]])
  cat("\nWARNING: EmptyDrops is DEPRECATED. Use BarcodeRanks methods instead.\n")
  cat("Running emptyDrops analysis...\n")
  cat("  Lower threshold:", args$emptydrops_lower, "\n")
  cat("  Testing FDR thresholds:", paste(fdr_thresholds, collapse=", "), "\n")
  
  set.seed(100)
  # Increase iterations for better p-value resolution when expecting fewer cells
  # Default is 10000, which gives minimum p-value of 1e-4
  # With 100000, minimum p-value is 1e-5, allowing better discrimination
  e_out <- emptyDrops(count_matrix, 
                      lower = args$emptydrops_lower, 
                      test.ambient = args$test_ambient,
                      niters = 100000,  # Increased from default 10000 for better precision
                      BPPARAM = BPPARAM)
} else {
  cat("Skipping emptyDrops analysis (not requested)\n")
}

# Clean up large objects after EmptyDrops
rm(count_matrix)
gc()

# Test multiple FDR thresholds
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
  br_rank = br_rank,  # BarcodeRanks rank (averaged for ties)
  br_fitted = fitted_values,  # BarcodeRanks fitted spline values
  log_prob = e_out$LogProb,
  p_value = e_out$PValue,
  fdr = e_out$FDR,
  is_cell_knee = is_cell_knee,
  is_cell_inflection = is_cell_inflection,
  limited = e_out$Limited
)

# Add FDR columns dynamically based on what was actually tested
for (fdr_thresh in fdr_thresholds) {
  col_name <- paste0("is_cell_fdr_", gsub("\\.", "_", as.character(fdr_thresh)))
  if (!is.null(fdr_results[[as.character(fdr_thresh)]])) {
    emptydrops_results[[col_name]] <- ifelse(is.na(fdr_results[[as.character(fdr_thresh)]]$is_cell), 
                                               FALSE, 
                                               fdr_results[[as.character(fdr_thresh)]]$is_cell)
  } else {
    emptydrops_results[[col_name]] <- FALSE
  }
}

write.table(emptydrops_results,
           file = file.path(args$output_dir, "dropletutils_results.tsv"),
           sep = "\t", row.names = FALSE, quote = FALSE)

# Save cell barcodes for each FDR threshold (only if EmptyDrops was run)
if (args$run_emptydrops) {
  for (fdr_thresh in fdr_thresholds) {
    fdr_str <- gsub("\\.", "", as.character(fdr_thresh))  # Remove decimal point for filename
    called_cells <- barcodes[fdr_results[[as.character(fdr_thresh)]]$is_cell & !is.na(fdr_results[[as.character(fdr_thresh)]]$is_cell)]
    writeLines(called_cells,
              file.path(args$output_dir, paste0(args$sample_id, "_EmptyDrops_FDR_", fdr_thresh, "_cell_barcodes.txt")))
  }
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

# Add EmptyDrops FDR results (only if EmptyDrops was run)
if (args$run_emptydrops) {
  for (i in seq_along(fdr_thresholds)) {
    fdr_thresh <- fdr_thresholds[i]
    result <- fdr_results[[as.character(fdr_thresh)]]
    
    summary_data[[length(summary_data) + 1]] <- data.frame(
      sample_id = args$sample_id,
      method = paste0("EmptyDrops_FDR_", fdr_thresh),
      fdr_threshold = fdr_thresh,
      n_cells_called = result$n_cells,
      threshold_used = result$threshold,
      umis_in_cells_pct = result$umis_in_cells_pct,
      emptydrops_lower_param = args$emptydrops_lower,
      n_barcodes_tested = sum(!is.na(e_out$PValue)),
      n_barcodes_limited = sum(e_out$Limited, na.rm = TRUE)
    )
  }
}

# Add barcodeRanks results (only if BarcodeRanks was run)
if (args$run_barcoderanks) {
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
}

# Write summary statistics if any methods were run
if (length(summary_data) > 0) {
  summary_stats <- do.call(rbind, summary_data)
  write.table(summary_stats,
             file = file.path(args$output_dir, "cell_calling_summary.tsv"),
             sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  # Create empty summary file with header
  empty_summary <- data.frame(
    sample_id = character(),
    method = character(),
    fdr_threshold = numeric(),
    n_cells_called = integer(),
    threshold_used = numeric(),
    umis_in_cells_pct = numeric(),
    emptydrops_lower_param = integer(),
    n_barcodes_tested = integer(),
    n_barcodes_limited = integer()
  )
  write.table(empty_summary,
             file = file.path(args$output_dir, "cell_calling_summary.tsv"),
             sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("DropletUtils analysis completed!\n")
cat("Results saved to:", args$output_dir, "\n")
if (args$run_emptydrops || args$run_barcoderanks) {
  cat("Cell barcodes saved for:\n")
  if (args$run_emptydrops) {
    cat("  - EmptyDrops FDR thresholds:", paste(fdr_thresholds, collapse=", "), "\n")
  }
  if (args$run_barcoderanks) {
    cat("  - BarcodeRanks knee and inflection points\n")
  }
}