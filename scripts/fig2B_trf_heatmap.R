# Figure 2C: Heatmap of deletion rates for 3p and 5p tRFs (mouse cerebellum)
#
# Requires pipeline_site_calling.R to have been run first.
# Inputs:  <data_dir>/means_df.csv
# Outputs: PNG heatmaps + CSV of underlying matrices
#
# Edit config.R: heatmap_tissue, heatmap_end, heatmap_type, heatmap_top_n, heatmap_aa

source("config.R")

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)
library(purrr)
library(ComplexHeatmap)

colors_custom <- c(
  colorRampPalette(c("white", "lightcoral"))(29),
  colorRampPalette(c("lightcoral", "darkred"))(70)
)

# -----------------------------------------------------------------------------
# Load pipeline output
# -----------------------------------------------------------------------------

means_df <- read_csv(file.path(pileup_dir, "means_df.csv"), show_col_types = FALSE)

# -----------------------------------------------------------------------------
# Heatmap function (uses local fragment coordinates; ComplexHeatmap)
# Returns the deletion-rate matrix invisibly; also writes PNG and CSV.
# -----------------------------------------------------------------------------

plot_trf_heatmaps <- function(df, tissue, end_type, frag_type = "cyto",
                               aa_filter = NULL, top_n = 50) {

  plot_data <- df %>% filter(Tissue == tissue, End_Type == end_type, Type == frag_type)
  if (!is.null(aa_filter)) plot_data <- plot_data %>% filter(AA %in% aa_filter)
  if (nrow(plot_data) == 0) {
    warning(sprintf(
      "No data for tissue='%s', end='%s', type='%s'.\n  Available Tissue values: %s\n  Available Type values:   %s",
      tissue, end_type, frag_type,
      paste(sort(unique(df$Tissue)), collapse = ", "),
      paste(sort(unique(df$Type)),   collapse = ", ")
    ))
    return(invisible(NULL))
  }

  create_matrix <- function(val_col) {
    mat <- plot_data %>%
      select(simplified_ref, start, !!sym(val_col)) %>%
      group_by(simplified_ref, start) %>%
      summarize(!!sym(val_col) := mean(!!sym(val_col), na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = start, values_from = !!sym(val_col)) %>%
      column_to_rownames("simplified_ref") %>%
      select(where(is.numeric)) %>%
      mutate(across(everything(), ~ if_else(is.na(.x), 0, .x)))
    as.matrix(mat[, as.character(sort(as.numeric(colnames(mat)))), drop = FALSE])
  }

  depth_mat <- create_matrix("mean_logcpm")
  dels_mat  <- create_matrix("mean_del_rate_cor")

  # Top N by depth at column closest to position 15
  sort_col <- as.character(min(15, max(as.numeric(colnames(depth_mat)))))
  if (sort_col %in% colnames(depth_mat)) {
    idx       <- order(depth_mat[, sort_col], decreasing = TRUE, na.last = TRUE)
    depth_mat <- depth_mat[idx, , drop = FALSE]
    dels_mat  <- dels_mat[rownames(depth_mat), , drop = FALSE]
  }
  if (!is.null(top_n) && nrow(depth_mat) > top_n) {
    depth_mat <- depth_mat[1:top_n, , drop = FALSE]
    dels_mat  <- dels_mat[1:top_n, , drop = FALSE]
  }

  # Sort alphabetically by isodecoder
  alpha_idx <- order(rownames(depth_mat))
  depth_mat <- depth_mat[alpha_idx, , drop = FALSE]
  dels_mat  <- dels_mat[rownames(depth_mat), , drop = FALSE]

  rownames(depth_mat) <- clean_labels(rownames(depth_mat))
  rownames(dels_mat)  <- clean_labels(rownames(dels_mat))

  # 1-based column tick labels (every 5)
  local_cols <- as.numeric(colnames(depth_mat)) + 1
  col_ticks  <- rep("", length(local_cols))
  col_ticks[which(local_cols %% 5 == 0)] <- local_cols[which(local_cols %% 5 == 0)]

  aa_tag  <- if (!is.null(aa_filter)) paste0("_", paste(aa_filter, collapse = "-")) else ""
  top_tag <- if (!is.null(top_n))     paste0("_top", top_n) else ""
  tag     <- paste0("tRF-", end_type, "_", frag_type, aa_tag, top_tag, "_", tissue)

  base_args <- list(cluster_cols = FALSE, cluster_rows = FALSE,
                    cellwidth = 10, cellheight = 10, fontsize_row = 10,
                    labels_col = col_ticks, angle_col = "0", na_col = "lightgrey")

  plot_height <- max(5, nrow(depth_mat) * 0.25)
  plot_width  <- 10

  # Coverage heatmap
  cov_png <- file.path(output_dir, paste0(tag, "_coverage.png"))
  png(cov_png, width = plot_width, height = plot_height, units = "in", res = 300)
  print(do.call(ComplexHeatmap::pheatmap,
                c(list(mat = depth_mat, breaks = seq(0, 12, length.out = 101),
                       main = paste(tissue, end_type, "Coverage")), base_args)))
  dev.off()

  # Deletion-rate heatmap
  del_png <- file.path(output_dir, paste0("fig2C_", tag, "_deletion_rates.png"))
  png(del_png, width = plot_width, height = plot_height, units = "in", res = 300)
  print(do.call(ComplexHeatmap::pheatmap,
                c(list(mat = dels_mat, color = colors_custom,
                       breaks = seq(0, 1, length.out = length(colors_custom) + 1),
                       main = paste(tissue, end_type, "Deletion Rates")), base_args)))
  dev.off()

  # CSV exports
  cov_csv <- sub("\\.png$", ".csv", cov_png)
  del_csv <- sub("\\.png$", ".csv", del_png)
  write.csv(depth_mat, cov_csv)
  write.csv(dels_mat,  del_csv)
  message("Saved: ", del_png, " + CSVs")

  invisible(dels_mat)
}

# -----------------------------------------------------------------------------
# Generate Figure 2C using parameters from config.R
# (call once per end type to match the figure panel)
# -----------------------------------------------------------------------------

plot_trf_heatmaps(means_df, heatmap_tissue, "3p", heatmap_type,
                  aa_filter = heatmap_aa, top_n = heatmap_top_n)
plot_trf_heatmaps(means_df, heatmap_tissue, "5p", heatmap_type,
                  aa_filter = heatmap_aa, top_n = heatmap_top_n)
