# Figure 2D: Heatmap of deletion rates at specific called Ψ sites across tissue replicates
#
# Requires pipeline_site_calling.R to have been run first (call_all_groups = TRUE).
# Inputs:  <data_dir>/called_site_data.csv
#          <data_dir>/modification_sites_all_groups.csv
# Outputs: PNG heatmap + CSV of underlying matrix
#
# Edit config.R: data_dir, fig2d_positions, fig2d_end, fig2d_type

source("config.R")

library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(purrr)
library(ComplexHeatmap)

colors_custom <- c(
  colorRampPalette(c("white", "lightcoral"))(29),
  colorRampPalette(c("lightcoral", "darkred"))(70)
)

# -----------------------------------------------------------------------------
# Load pipeline outputs
# -----------------------------------------------------------------------------

called_site_data <- read_csv(file.path(pileup_dir, "called_site_data.csv"),
                              show_col_types = FALSE)
all_mod_sites    <- read_csv(file.path(pileup_dir, "modification_sites_all_groups.csv"),
                              show_col_types = FALSE)

# -----------------------------------------------------------------------------
# Build per-sample master table with pairwise t-test statistics
# -----------------------------------------------------------------------------

called_ids <- all_mod_sites %>%
  mutate(
    End_Type = if_else(str_detect(ref, "3p$"), "3p",
                       if_else(str_detect(ref, "5p$"), "5p", "other")),
    site_id  = paste(simplified_ref, End_Type, pos, sep = "_")
  ) %>%
  pull(site_id) %>% unique()

ttest_data <- called_site_data %>%
  mutate(site_id = paste(simplified_ref, End_Type, pos, sep = "_")) %>%
  filter(site_id %in% called_ids) %>%
  group_by(site_id, simplified_ref, End_Type, pos, Tissue, Sample_ID) %>%
  summarize(del_rate = mean(del_rate, na.rm = TRUE), .groups = "drop")

sample_rates <- ttest_data %>%
  select(site_id, Sample_ID, del_rate) %>%
  pivot_wider(names_from = Sample_ID, values_from = del_rate)

mean_rates <- ttest_data %>%
  group_by(site_id, simplified_ref, End_Type, pos, Tissue) %>%
  summarize(mean_del = mean(del_rate, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Tissue, values_from = mean_del, names_prefix = "Mean_")

master_table <- left_join(mean_rates, sample_rates, by = "site_id")

# Pairwise t-tests across tissue groups
groups <- unique(ttest_data$Tissue)
pairs  <- combn(groups, 2, simplify = FALSE)

for (pair in pairs) {
  g1 <- pair[1]; g2 <- pair[2]
  stats <- ttest_data %>%
    filter(Tissue %in% c(g1, g2)) %>%
    group_by(site_id) %>%
    summarize(
      p_val = {
        v1 <- del_rate[Tissue == g1]; v2 <- del_rate[Tissue == g2]
        if (sum(!is.na(v1)) >= 2 && sum(!is.na(v2)) >= 2 &&
            (var(v1, na.rm = TRUE) + var(v2, na.rm = TRUE)) > 0)
          t.test(v1, v2)$p.value else NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(!!paste0("padj_", g2, "_vs_", g1) := p.adjust(p_val, method = "BH")) %>%
    select(-p_val)
  master_table <- left_join(master_table, stats, by = "site_id")
}

write_csv(master_table, file.path(output_dir, "Master_Differential_Modification_Sites.csv"))
message("Master table → Master_Differential_Modification_Sites.csv")

# -----------------------------------------------------------------------------
# Heatmap function: plot deletion rates at specific positions across replicates
# Returns heatmap matrix invisibly; writes PNG and CSV.
# -----------------------------------------------------------------------------

known_cols  <- c("site_id", "simplified_ref", "End_Type", "pos")
stat_cols   <- grep("^(Mean|Diff|Log2FC|pval|padj)_", colnames(master_table), value = TRUE)
sample_cols <- setdiff(colnames(master_table), c(known_cols, stat_cols))

plot_position_heatmap <- function(master_df, target_pos, target_end = "3p",
                                   target_type = "cyto", top_n = 30) {

  pos_label <- paste(target_pos, collapse = "-")

  plot_data <- master_df %>%
    filter(pos %in% target_pos, End_Type == target_end) %>%
    mutate(Type = if_else(str_detect(simplified_ref, "mito|mt-"), "mito", "cyto")) %>%
    filter(Type == target_type)

  if (nrow(plot_data) == 0) {
    message("No called sites at position(s) ", pos_label, " in ",
            target_end, " ", target_type, " fragments.")
    return(invisible(NULL))
  }

  # Average across positions if multiple, keep lowest padj per tRF
  padj_cols <- grep("^padj_", colnames(plot_data), value = TRUE)
  plot_data <- plot_data %>%
    group_by(simplified_ref) %>%
    summarize(across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE)),
              across(all_of(padj_cols),   ~ min(.x,  na.rm = TRUE)),
              .groups = "drop") %>%
    rowwise() %>%
    mutate(min_padj = min(c_across(all_of(padj_cols)), na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(min_padj) %>%
    slice_head(n = top_n)

  mat <- plot_data %>%
    select(simplified_ref, all_of(sample_cols)) %>%
    column_to_rownames("simplified_ref") %>%
    as.matrix()

  # Order columns by group according to sample_order
  if (!is.null(sample_order)) {
    grp_map   <- setNames(ttest_data$Tissue, ttest_data$Sample_ID)
    col_order <- order(match(grp_map[colnames(mat)], sample_order))
    mat       <- mat[, col_order, drop = FALSE]
  }

  plot_height <- max(4, nrow(mat) * 0.25)
  tag     <- paste0("Pos", pos_label, "_", target_end, "_", target_type)
  out_png <- file.path(output_dir, paste0("fig2D_Heatmap_", tag, ".png"))
  out_csv <- sub("\\.png$", ".csv", out_png)

  png(out_png, width = 4, height = plot_height, units = "in", res = 300)
  p <- ComplexHeatmap::pheatmap(
    mat,
    cluster_rows   = FALSE, cluster_cols = FALSE,
    breaks         = seq(0, 1, length.out = length(colors_custom) + 1),
    color          = colors_custom,
    border_color   = "grey80",
    main           = paste0("Ψ deletion rates at pos", pos_label,
                             "\n(", target_end, " ", target_type, " tRFs)"),
    row_names_side = "left"
  )
  print(p)
  dev.off()

  write.csv(mat, out_csv)
  message("Saved: ", out_png, " + CSV")

  invisible(mat)
}

# -----------------------------------------------------------------------------
# Generate Figure 2E using parameters from config.R
# -----------------------------------------------------------------------------

plot_position_heatmap(master_table,
                      target_pos  = fig2d_positions,
                      target_end  = fig2d_end,
                      target_type = fig2d_type)
