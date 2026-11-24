# Figure 3E: log2(CPM) vs mean deletion rate for miRNA T positions, per mouse tissue
# One PNG + CSV output per tissue group
#
# Inputs:  <ncrna_mm_dir>/<sample>.csv  or <sample>_annotated.csv
#          <data_dir>/samples_mm.csv  (cols: Sample, Treatment [NT/RBS], Group)
# Outputs per tissue: fig3E_<Group>_mirna_deletion_scatter.png + .csv

source("config.R")

library(dplyr)
library(readr)
library(ggplot2)
library(purrr)

# =============================================================================
# Configuration
# =============================================================================

MIN_DEPTH   <- 20          # CF112: MIN_DEPTH
NT_DEL_MAX  <- NT_cutoff   # CF112: NT_DEL_MAX  (config: 0.01)
RBS_DEL_MIN <- fig3e_repro_del_min  # CF112: RBS_DEL_MIN (config: 0.20)

SITE_KEYS <- c("ref", "start", "end", "base", "external_gene_name", "gene_biotype")

# =============================================================================
# Helper functions
# =============================================================================

add_cpm <- function(df) {
  total      <- sum(df$depth, na.rm = TRUE)
  df$cpm     <- 1e6 * df$depth / total
  df$log2cpm <- log2(df$cpm + 1)
  df
}

make_wide <- function(df_list, site_keys, prefix, extra_cols = character()) {
  imap(df_list, function(df, i) {
    keep <- c(site_keys, "depth", "del_rate", extra_cols)
    df %>%
      select(any_of(keep)) %>%
      rename(
        !!paste0("depth_",    prefix, "-", i) := depth,
        !!paste0("del_rate_", prefix, "-", i) := del_rate
      ) %>%
      {
        if ("log2cpm" %in% extra_cols)
          rename(., !!paste0("log2cpm_", prefix, "-", i) := log2cpm)
        else .
      }
  }) %>%
    reduce(full_join, by = site_keys)
}

load_pileup <- function(sample_name) {
  annotated <- file.path(ncrna_mm_dir, paste0(sample_name, "_annotated.csv"))
  plain     <- file.path(ncrna_mm_dir, paste0(sample_name, ".csv"))
  p <- if (file.exists(annotated)) annotated else plain
  if (!file.exists(p)) stop("Pileup file not found: ", plain)
  add_cpm(read_csv(p, show_col_types = FALSE))
}

# =============================================================================
# Load sample sheet
# =============================================================================

all_samples <- read_csv(file.path(data_dir, samples_mm_csv), show_col_types = FALSE)

groups <- fig3d_groups
if (is.null(groups)) groups <- unique(all_samples$Group)

# =============================================================================
# Main loop over tissue groups
# =============================================================================

for (GROUP_FILTER in groups) {

  message("\n========== Group: ", GROUP_FILTER, " ==========")

  samples     <- all_samples %>% filter(Group == GROUP_FILTER)
  nt_samples  <- samples %>% filter(Treatment == "NT")
  rbs_samples <- samples %>% filter(Treatment != "NT")

  if (nrow(nt_samples) == 0 || nrow(rbs_samples) == 0) {
    message("Skipping ", GROUP_FILTER, ": missing NT or RBS samples.")
    next
  }

  n_nt  <- nrow(nt_samples)
  n_rbs <- nrow(rbs_samples)

  nt_list  <- lapply(nt_samples$Sample,  load_pileup)
  rbs_list <- lapply(rbs_samples$Sample, load_pileup)

  nt_wide  <- make_wide(nt_list,  SITE_KEYS, "NT")
  rbs_wide <- make_wide(rbs_list, SITE_KEYS, "RBS", extra_cols = "log2cpm")
  combined <- full_join(nt_wide, rbs_wide, by = SITE_KEYS)

  nt_depth_cols   <- paste0("depth_NT-",    1:n_nt)
  nt_del_cols     <- paste0("del_rate_NT-", 1:n_nt)
  rbs_depth_cols  <- paste0("depth_RBS-",   1:n_rbs)
  rbs_del_cols    <- paste0("del_rate_RBS-",1:n_rbs)
  rbs_logcpm_cols <- paste0("log2cpm_RBS-", 1:n_rbs)

  # ---------------------------------------------------------------------------
  # Build plot_df: miRNA T positions, annotated only
  # ---------------------------------------------------------------------------

  plot_df <- combined %>%
    filter(base == "T") %>%
    mutate(
      mean_depth_NT  = rowMeans(pick(all_of(nt_depth_cols)),   na.rm = TRUE),
      mean_depth_RBS = rowMeans(pick(all_of(rbs_depth_cols)),  na.rm = TRUE),
      mean_del_NT    = rowMeans(pick(all_of(nt_del_cols)),     na.rm = TRUE),
      mean_del_RBS   = rowMeans(pick(all_of(rbs_del_cols)),    na.rm = TRUE),
      mean_log2cpm   = rowMeans(pick(all_of(rbs_logcpm_cols)), na.rm = TRUE),
      reproducible   = rowSums(across(all_of(rbs_del_cols), ~ !is.na(.) & . > RBS_DEL_MIN)) >= n_rbs
    ) %>%
    filter(
      mean_depth_NT  >= MIN_DEPTH,
      mean_depth_RBS >= MIN_DEPTH,
      mean_del_NT    <  NT_DEL_MAX,
      gene_biotype   == "miRNA",
      !is.na(external_gene_name)
    )

  if (nrow(plot_df) == 0) {
    message("No miRNA T positions passed filters for group: ", GROUP_FILTER)
    next
  }

  message(sum(plot_df$reproducible), " reproducible sites in ", GROUP_FILTER)

  # ---------------------------------------------------------------------------
  # Scatter plot
  # ---------------------------------------------------------------------------

  p <- ggplot(plot_df, aes(x = mean_log2cpm, y = mean_del_RBS)) +
    geom_hline(yintercept = RBS_DEL_MIN, linetype = "dashed", color = "grey40") +
    geom_point(aes(color = reproducible), alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey80"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    scale_y_continuous(limits = c(0, 0.8)) +
    labs(
      x = expression(log[2](CPM)),
      y = "Mean deletion rate (RBS)"
    ) +
    theme_classic(base_size = 12)

  out_png <- file.path(output_dir,
                       paste0("fig3E_", GROUP_FILTER, "_mirna_deletion_scatter.png"))
  out_csv <- sub("\\.png$", ".csv", out_png)

  ggsave(out_png, p, width = 5, height = 4, dpi = 300, bg = "white")
  write_csv(plot_df %>% select(ref, start, base, external_gene_name,
                               mean_del_NT, mean_del_RBS, mean_log2cpm, reproducible),
            out_csv)
  message("Saved: ", basename(out_png), " + CSV")
}
