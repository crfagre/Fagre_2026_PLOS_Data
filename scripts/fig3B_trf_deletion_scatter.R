# Figure 3B: tRF deletion rate scatter (human HepG2)

source("config.R")

library(dplyr)
library(readr)
library(ggplot2)
library(purrr)
library(tidyr)

# =============================================================================
# Configuration
# =============================================================================

TRF_DIR       <- trf_hs_dir
SAMPLE_SHEET  <- file.path(data_dir, samples_hs_csv)
GROUP_FILTER  <- fig3b_tissue     # NULL = use all samples
MIN_DEPTH     <- 20
NT_DEL_MAX    <- 0.01
RBS_DEL_MIN   <- 0.20
REPRO_DEL_MIN <- 0.2

TRF_SITE_KEYS <- c("ref", "start", "end", "base")

# =============================================================================
# Load tRF FASTA (reference sequences for site calling)
# =============================================================================

read_fasta <- function(path) {
  lines      <- readLines(path)
  header_idx <- which(startsWith(lines, ">"))
  headers    <- sub("^>", "", lines[header_idx])
  seq_ends   <- c(header_idx[-1] - 1L, length(lines))
  seqs       <- mapply(
    function(s, e) paste(lines[(s + 1L):e], collapse = ""),
    header_idx, seq_ends
  )
  setNames(seqs, headers)
}

trf_fa <- read_fasta(file.path(ref_dir, trf_hs_fasta_file))

# =============================================================================
# Load sample sheet and split by treatment
# =============================================================================

load_pileup_from <- function(filename, dir) {
  read_csv(file.path(dir, filename), show_col_types = FALSE)
}

add_cpm <- function(df) {
  total      <- sum(df$depth, na.rm = TRUE)
  df$cpm     <- 1e6 * df$depth / total
  df$log2cpm <- log2(df$cpm + 1)
  df
}

trf_samples <- read_csv(SAMPLE_SHEET, show_col_types = FALSE) %>%
  { if (!is.null(GROUP_FILTER)) filter(., Group == GROUP_FILTER) else . } %>%
  mutate(filename = paste0(Sample, ".csv"))

trf_nt_samples  <- trf_samples %>% filter(Treatment == "NT")
trf_rbs_samples <- trf_samples %>% filter(Treatment == "RBS")

trf_nt_list  <- lapply(trf_nt_samples$filename,  function(f) add_cpm(load_pileup_from(f, TRF_DIR)))
trf_rbs_list <- lapply(trf_rbs_samples$filename, function(f) add_cpm(load_pileup_from(f, TRF_DIR)))

trf_n_nt  <- length(trf_nt_list)
trf_n_rbs <- length(trf_rbs_list)

# =============================================================================
# Merge replicates into wide format (one row per position)
# =============================================================================

make_wide_trf <- function(df_list, prefix, extra_cols = character()) {
  imap(df_list, function(df, i) {
    keep <- c(TRF_SITE_KEYS, "depth", "del_rate", extra_cols)
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
    reduce(full_join, by = TRF_SITE_KEYS)
}

trf_nt_wide  <- make_wide_trf(trf_nt_list,  "NT")
trf_rbs_wide <- make_wide_trf(trf_rbs_list, "RBS", extra_cols = "log2cpm")
trf_combined <- full_join(trf_nt_wide, trf_rbs_wide, by = TRF_SITE_KEYS)

trf_nt_depth_cols   <- paste0("depth_NT-",     1:trf_n_nt)
trf_nt_del_cols     <- paste0("del_rate_NT-",  1:trf_n_nt)
trf_rbs_depth_cols  <- paste0("depth_RBS-",    1:trf_n_rbs)
trf_rbs_del_cols    <- paste0("del_rate_RBS-", 1:trf_n_rbs)
trf_rbs_logcpm_cols <- paste0("log2cpm_RBS-",  1:trf_n_rbs)

# =============================================================================
# Site calling
# =============================================================================

trf_coverage_sites <- trf_combined %>%
  filter(base == "T") %>%
  filter(if_all(all_of(trf_nt_depth_cols),  ~ !is.na(.) & . >= MIN_DEPTH)) %>%
  filter(if_all(all_of(trf_rbs_depth_cols), ~ !is.na(.) & . >= MIN_DEPTH)) %>%
  mutate(
    passes_nt_del      = if_all(all_of(trf_nt_del_cols), ~ !is.na(.) & . < NT_DEL_MAX),
    passes_rbs_del     = rowSums(across(all_of(trf_rbs_del_cols), ~ !is.na(.) & . > RBS_DEL_MIN)) >= 2,
    psi_site           = passes_nt_del & passes_rbs_del,
    reference_sequence = trf_fa[ref]
  ) %>%
  select(-passes_nt_del, -passes_rbs_del) %>%
  rename(position = start) %>%
  select(-end)

group_tag <- if (!is.null(GROUP_FILTER)) GROUP_FILTER else "all"

trf_all_sites_file <- file.path(output_dir, paste0(group_tag, "_tRF_all_covered_sites.csv"))
write_csv(trf_coverage_sites, trf_all_sites_file)
message("Written ", nrow(trf_coverage_sites), " covered tRF T positions → ", trf_all_sites_file,
        " (", sum(trf_coverage_sites$psi_site), " positive hits)")

trf_called_sites <- filter(trf_coverage_sites, psi_site)
trf_out_file <- file.path(output_dir, paste0(group_tag, "_tRF_called_sites.csv"))
write_csv(trf_called_sites, trf_out_file)
message("Called ", nrow(trf_called_sites), " tRF sites → ", trf_out_file)

# =============================================================================
# Scatter plot: log2(CPM) vs mean RBS deletion rate
# =============================================================================

trf_plot_df <- trf_combined %>%
  filter(base == "T") %>%
  mutate(
    mean_depth_NT = rowMeans(pick(all_of(trf_nt_depth_cols)),   na.rm = TRUE),
    mean_del_NT   = rowMeans(pick(all_of(trf_nt_del_cols)),     na.rm = TRUE),
    mean_del_RBS  = rowMeans(pick(all_of(trf_rbs_del_cols)),    na.rm = TRUE),
    mean_log2cpm  = rowMeans(pick(all_of(trf_rbs_logcpm_cols)), na.rm = TRUE),
    reproducible  = rowSums(across(all_of(trf_rbs_del_cols), ~ !is.na(.) & . > REPRO_DEL_MIN)) >= 2
  ) %>%
  filter(
    mean_depth_NT >= MIN_DEPTH,
    mean_del_NT   <  NT_DEL_MAX
  )

if (nrow(trf_plot_df) == 0)
  stop("No tRF T positions passed all filters.")

message(sum(trf_plot_df$reproducible), " reproducible tRF sites out of ",
        nrow(trf_plot_df), " covered positions")

# --- Reproducibility-colored scatter (fig3B) ---------------------------------

trf_plot <- ggplot(trf_plot_df, aes(x = mean_log2cpm, y = mean_del_RBS)) +
  geom_hline(yintercept = REPRO_DEL_MIN, linetype = "dashed", color = "grey40") +
  geom_point(aes(color = reproducible), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey80"), guide = "none") +
  scale_y_continuous(limits = c(0, 0.8)) +
  labs(
    title = paste0(group_tag, " - tRF"),
    x     = expression(log[2](CPM)),
    y     = "Mean deletion rate (RBS)"
  ) +
  theme_classic(base_size = 12)

trf_plot_file <- file.path(output_dir, paste0("fig3B_trf_scatter_", group_tag, ".png"))
ggsave(trf_plot_file, trf_plot, width = 5, height = 4, dpi = 300, bg = "white")
message("Saved: ", trf_plot_file)
print(trf_plot)

# --- CSV export --------------------------------------------------------------

export_cols <- c("ref", "start", "base", "mean_del_NT", "mean_del_RBS", "mean_log2cpm")
out_csv <- file.path(output_dir, paste0("fig3B_trf_scatter_", group_tag, ".csv"))
write_csv(trf_plot_df %>% select(any_of(export_cols)), out_csv)
message("Saved: ", out_csv)
