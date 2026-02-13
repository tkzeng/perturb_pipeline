library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

result_dir <- "/Users/wkq953/Desktop/CBMR/Projects/Exchange/exchange_project/analysis/data/sequencing/results_eRZ58_40B_analysis"
split_pipe_GEX <- file.path(result_dir, "GEX_agg_samp_ana_summary.csv") %>% fread() 
split_pipe_guide <- fread("gRNA_agg_samp_ana_summary.csv") 
# perturb_pipe_raw <- fread("/Users/wkq953/Desktop/CBMR/Projects/Exchange/exchange_project/analysis/data/sequencing/eRZ58_20A/GEX_by_biological_sample.tsv")
perturb_pipe <- fread("qc_summary.csv")

split_stat_to_sample <- function(dt) {
  stopifnot("statistic" %in% names(dt))
  long <- melt(
    dt,
    id.vars = "statistic",
    variable.name = "biological_sample",
    value.name = "value"
  )
  
  wide <- dcast(
    long,
    biological_sample ~ statistic,
    value.var = "value"
  )
  
  # optional: keep sample_id as first column and sort
  setcolorder(wide, c("biological_sample", setdiff(names(wide), "biological_sample")))
  setorder(wide, biological_sample)
  
  return(wide)
}

split_pipe   <- split_stat_to_sample(split_pipe_GEX) %>%
  left_join(split_stat_to_sample(split_pipe_guide), by = "biological_sample", suffix = c("_GEX", "_guide"))

# > colnames(perturb_pipe)
# [1] "biological_sample"          "n_cells_per_sample"         "n_cells_with_guides"        "gex_umis_per_sample"        "gex_umis_per_cell_mean"    
# [6] "gex_umis_per_cell_median"   "genes_per_cell_mean"        "genes_per_cell_median"      "grna_umis_per_sample"       "grna_umis_per_cell_mean"   
# [11] "grna_umis_per_cell_median"  "grna_detections_per_sample" "grna_per_cell_mean"         "grna_per_cell_median"       "frac_cells_with_guides"    
# [16] "pct_mt_mean"                "pct_mt_median"              "pct_ribo_mean"              "pct_ribo_median"           

# QC metric -> extract corresponding columns in splite pipe

split_pipe_clean <- split_pipe %>%
  transmute(
     biological_sample = biological_sample,
     n_cells_per_sample = number_of_cells_GEX, 
     n_cells_with_guides  = gRNA_number_of_cells,
     gex_umis_per_sample        = number_of_tscp_GEX,
     gex_umis_per_cell_median   = `GRCm39-Cre-mCherry_median_tscp_per_cell`,
     genes_per_cell_median      = `GRCm39-Cre-mCherry_median_genes_per_cell`,
     grna_umis_per_sample       = number_of_tscp_guide,
     grna_umis_per_cell_median  = gRNA_median_tscp_per_cell,
     grna_detections_per_sample = detected_number_guides,
     frac_cells_with_guides     = 1 - frac_cells_with_0_crispr_guide
     )

metrics <- c(
  "biological_sample",
  "n_cells_per_sample",
  "n_cells_with_guides",
  "gex_umis_per_sample",
  "gex_umis_per_cell_median",
  "genes_per_cell_median",
  "grna_umis_per_sample",
  "grna_umis_per_cell_median",
  "grna_detections_per_sample",
  "frac_cells_with_guides"
)

merged <- split_pipe_clean %>%
  inner_join(perturb_pipe %>% select(all_of(c("biological_sample", metrics))),
             by = "biological_sample",
             suffix = c("_split", "_perturb"))

plot_dt <- merged %>%
  pivot_longer(
    cols = ends_with("_split"),
    names_to = "metric",
    values_to = "split_value"
  ) %>%
  mutate(metric = sub("_split$", "", metric)) %>%
  left_join(
    merged %>%
      pivot_longer(cols = ends_with("_perturb"),
                   names_to = "metric2",
                   values_to = "perturb_value") %>%
      mutate(metric = sub("_perturb$", "", metric2)) %>%
      select(biological_sample, metric, perturb_value),
    by = c("biological_sample", "metric")
  ) %>%
  filter(metric %in% metrics) %>%
  mutate(
    split_value = as.numeric(split_value),
    perturb_value = as.numeric(perturb_value)
  ) %>%
  filter(is.finite(split_value) & is.finite(perturb_value))


p <- ggplot(plot_dt, aes(x = split_value, y = perturb_value)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  facet_wrap(~ metric, scales = "free") +
  labs(x = "Split-pipe", y = "Perturb-pipe", title = "QC metrics: Split-pipe vs Perturb-pipe (per sample)") +
  theme_bw(base_size = 12)

print(p)
ggsave(file.path(result_dir,"/split_vs_perturb_scatter_by_metric.png"), p, width = 12, height = 8, dpi = 300)
