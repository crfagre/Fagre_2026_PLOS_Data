# Figure 3B: miRNA deletion rate scatter (human HepG2)

source("config.R")

library(dplyr)
library(readr)
library(ggplot2)
library(purrr)
library(tidyr)

# =============================================================================
# Configuration
# =============================================================================

SAMPLE_DIR    <- ncrna_hs_dir
SAMPLE_SHEET  <- file.path(data_dir, samples_hs_csv)
GROUP_FILTER  <- fig3b_tissue     # NULL = use all samples
MIN_DEPTH     <- 20
NT_DEL_MAX    <- 0.01
REPRO_DEL_MIN <- 0.2

# =============================================================================
# Load sample sheet and split by treatment
# =============================================================================

samples <- read_csv(SAMPLE_SHEET, show_col_types = FALSE) %>%
  { if (!is.null(GROUP_FILTER)) filter(., Group == GROUP_FILTER) else . } %>%
  mutate(filename = paste0(Sample, "_annotated.csv"))

nt_samples  <- samples %>% filter(Treatment == "NT")
rbs_samples <- samples %>% filter(Treatment == "RBS")

stopifnot(nrow(nt_samples)  > 0)
stopifnot(nrow(rbs_samples) > 0)

n_nt  <- nrow(nt_samples)
n_rbs <- nrow(rbs_samples)

# =============================================================================
# Load pileup files and compute CPM
# =============================================================================

load_pileup_from <- function(filename, dir) {
  read_csv(file.path(dir, filename), show_col_types = FALSE)
}
load_pileup <- function(filename) load_pileup_from(filename, SAMPLE_DIR)

add_cpm <- function(df) {
  total      <- sum(df$depth, na.rm = TRUE)
  df$cpm     <- 1e6 * df$depth / total
  df$log2cpm <- log2(df$cpm + 1)
  df
}

nt_list  <- lapply(nt_samples$filename,  function(f) add_cpm(load_pileup(f)))
rbs_list <- lapply(rbs_samples$filename, function(f) add_cpm(load_pileup(f)))

# =============================================================================
# Merge replicates into wide format (one row per position)
# =============================================================================

SITE_KEYS <- c("ref", "start", "end", "base", "external_gene_name", "gene_biotype")

make_wide <- function(df_list, prefix, extra_cols = character()) {
  imap(df_list, function(df, i) {
    keep <- c(SITE_KEYS, "depth", "del_rate", extra_cols)
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
    reduce(full_join, by = SITE_KEYS)
}

nt_wide  <- make_wide(nt_list,  "NT")
rbs_wide <- make_wide(rbs_list, "RBS", extra_cols = "log2cpm")

combined <- full_join(nt_wide, rbs_wide, by = SITE_KEYS)

nt_depth_cols   <- paste0("depth_NT-",    1:n_nt)
nt_del_cols     <- paste0("del_rate_NT-", 1:n_nt)
rbs_depth_cols  <- paste0("depth_RBS-",   1:n_rbs)
rbs_del_cols    <- paste0("del_rate_RBS-",1:n_rbs)
rbs_logcpm_cols <- paste0("log2cpm_RBS-", 1:n_rbs)

# =============================================================================
# Scatter plot: log2(CPM) vs mean RBS deletion rate
# =============================================================================

paired_ttest_p <- function(rbs_vals, nt_vals) {
  tryCatch(
    t.test(rbs_vals, nt_vals, paired = TRUE, alternative = "greater")$p.value,
    error = function(e) NA_real_
  )
}

plot_df <- combined %>%
  filter(base == "T") %>%
  mutate(
    mean_depth_NT = rowMeans(pick(all_of(nt_depth_cols)),   na.rm = TRUE),
    mean_del_NT   = rowMeans(pick(all_of(nt_del_cols)),     na.rm = TRUE),
    mean_del_RBS  = rowMeans(pick(all_of(rbs_del_cols)),    na.rm = TRUE),
    mean_log2cpm  = rowMeans(pick(all_of(rbs_logcpm_cols)), na.rm = TRUE),
    reproducible  = rowSums(across(all_of(rbs_del_cols), ~ !is.na(.) & . > REPRO_DEL_MIN)) >= 2
  ) %>%
  filter(
    mean_depth_NT >= MIN_DEPTH,
    mean_del_NT   <  NT_DEL_MAX,
    gene_biotype  == "miRNA"
  ) %>%
  rowwise() %>%
  mutate(ttest_p = paired_ttest_p(
    c_across(all_of(rbs_del_cols)),
    c_across(all_of(nt_del_cols))
  )) %>%
  ungroup()

if (nrow(plot_df) == 0)
  stop("No miRNA T positions passed all filters.")

message(sum(plot_df$reproducible), " reproducible sites out of ",
        nrow(plot_df), " covered positions")

group_tag <- if (!is.null(GROUP_FILTER)) GROUP_FILTER else "all"

# --- Reproducibility-colored scatter -----------------------------------------

p <- ggplot(plot_df, aes(x = mean_log2cpm, y = mean_del_RBS)) +
  geom_hline(yintercept = REPRO_DEL_MIN, linetype = "dashed", color = "grey40") +
  geom_point(aes(color = reproducible), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey80"), guide = "none") +
  scale_y_continuous(limits = c(0, 0.8)) +
  labs(
    title = paste0(group_tag, " - miRNA"),
    x     = expression(log[2](CPM)),
    y     = "Mean deletion rate (RBS)"
  ) +
  theme_classic(base_size = 12)

out_png <- file.path(output_dir, paste0("fig3B_mirna_scatter_", group_tag, ".png"))
ggsave(out_png, p, width = 5, height = 4, dpi = 300, bg = "white")
message("Saved: ", out_png)
print(p)

# --- t-test-colored scatter --------------------------------------------------

p_ttest <- ggplot(plot_df, aes(x = mean_log2cpm, y = mean_del_RBS)) +
  geom_hline(yintercept = REPRO_DEL_MIN, linetype = "dashed", color = "grey40") +
  geom_point(aes(color = !is.na(ttest_p) & ttest_p < 0.01), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey80"), guide = "none") +
  scale_y_continuous(limits = c(0, 0.8)) +
  labs(
    title = paste0(group_tag, " - miRNA (paired t-test p<0.01)"),
    x     = expression(log[2](CPM)),
    y     = "Mean deletion rate (RBS)"
  ) +
  theme_classic(base_size = 12)

out_ttest_png <- file.path(output_dir, paste0("fig3B_mirna_scatter_ttest_", group_tag, ".png"))
ggsave(out_ttest_png, p_ttest, width = 5, height = 4, dpi = 300, bg = "white")
message("Saved: ", out_ttest_png)

# --- CSV export --------------------------------------------------------------

export_cols <- c("ref", "start", "base", "external_gene_name",
                 "mean_del_NT", "mean_del_RBS", "mean_log2cpm")
out_csv <- file.path(output_dir, paste0("fig3B_mirna_scatter_", group_tag, ".csv"))
write_csv(plot_df %>% select(any_of(export_cols)), out_csv)
message("Saved: ", out_csv)
