# Figure 1G: Heatmap of deletion rates for the top N expressed tRFs, grouped by isodecoders
#
# Requires pipeline_site_calling.R to have been run first.
# Edit config.R to set heatmap_tissue, heatmap_end, heatmap_type, heatmap_top_n, heatmap_aa.

source("config.R")

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(pheatmap)
library(tibble)
library(purrr)

colors_custom <- c(
  colorRampPalette(c("white", "lightcoral"))(29),
  colorRampPalette(c("lightcoral", "darkred"))(70)
)

#------------------------------------------------------------------------------
# Load pipeline output and position metadata
#------------------------------------------------------------------------------

means_df <- read_csv(file.path(pileup_dir, "means_df.csv"), show_col_types = FALSE)

ref_lookup   <- read_csv(file.path(ref_dir, trna_lookup_file),    show_col_types = FALSE)
sprinzl_meta <- read_csv(file.path(ref_dir, sprinzl_lookup_file), show_col_types = FALSE) %>%
  mutate(start = local_pos - 1) %>%
  select(ref, start, sprinzl_pos)

valid_positions <- sprinzl_meta %>%
  left_join(ref_lookup, by = c("ref" = "original_id")) %>%
  select(simplified_ref, sprinzl_pos) %>%
  distinct() %>%
  drop_na(simplified_ref, sprinzl_pos)

#------------------------------------------------------------------------------
# Heatmap function
#------------------------------------------------------------------------------

plot_trf_heatmaps <- function(df, tissue, end_type, frag_type = "cyto",
                               aa_filter = NULL, top_n = 50) {

  plot_data <- df %>%
    filter(Tissue == tissue, End_Type == end_type, Type == frag_type)
  if (!is.null(aa_filter)) plot_data <- plot_data %>% filter(AA %in% aa_filter)
  if (nrow(plot_data) == 0) { warning("No data for this filter combination."); return(invisible(NULL)) }

  create_matrix <- function(val_col) {
    if (end_type == "5p") {
      mat_data <- plot_data %>%
        select(simplified_ref, sprinzl_pos, !!sym(val_col)) %>%
        group_by(simplified_ref, sprinzl_pos) %>%
        summarize(!!sym(val_col) := mean(!!sym(val_col), na.rm = TRUE), .groups = "drop")
      mat <- valid_positions %>%
        filter(simplified_ref %in% plot_data$simplified_ref) %>%
        left_join(mat_data, by = c("simplified_ref", "sprinzl_pos")) %>%
        mutate(!!sym(val_col) := if_else(is.na(!!sym(val_col)), 0, !!sym(val_col))) %>%
        pivot_wider(names_from = sprinzl_pos, values_from = !!sym(val_col)) %>%
        column_to_rownames("simplified_ref") %>%
        select(where(is.numeric))
    } else {
      mat <- plot_data %>%
        select(simplified_ref, start, !!sym(val_col)) %>%
        group_by(simplified_ref, start) %>%
        summarize(!!sym(val_col) := mean(!!sym(val_col), na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = start, values_from = !!sym(val_col)) %>%
        column_to_rownames("simplified_ref") %>%
        select(where(is.numeric)) %>%
        mutate(across(everything(), ~ if_else(is.na(.x), 0, .x)))
    }
    mat[, as.character(sort(as.numeric(colnames(mat)))), drop = FALSE]
  }

  depth_mat <- create_matrix("mean_logcpm")
  dels_mat  <- create_matrix("mean_del_rate_cor")

  # Select top N by depth at column 15 (or last available)
  sort_col         <- min(15, ncol(depth_mat))
  depth_sorted_idx <- order(depth_mat[[sort_col]], decreasing = TRUE, na.last = TRUE)
  depth_mat        <- depth_mat[depth_sorted_idx, , drop = FALSE]
  dels_mat         <- dels_mat[rownames(depth_mat), , drop = FALSE]
  if (!is.null(top_n) && nrow(depth_mat) > top_n) {
    depth_mat <- depth_mat[1:top_n, , drop = FALSE]
    dels_mat  <- dels_mat[1:top_n, , drop = FALSE]
  }

  # Clean and sort row labels alphabetically by isodecoder
  rownames(depth_mat) <- clean_labels(rownames(depth_mat))
  rownames(dels_mat)  <- clean_labels(rownames(dels_mat))
  alpha_idx           <- order(rownames(depth_mat))
  depth_mat           <- depth_mat[alpha_idx, , drop = FALSE]
  dels_mat            <- dels_mat[rownames(depth_mat), , drop = FALSE]

  # Column tick labels (every 5 positions)
  n_cols    <- ncol(depth_mat)
  col_ticks <- rep("", n_cols)
  if (end_type == "3p") {
    pos_seq        <- -(n_cols):-1
    idx            <- which(pos_seq %% 5 == 0)
    col_ticks[idx] <- pos_seq[idx]
  } else {
    pos_seq        <- as.numeric(colnames(depth_mat))
    idx            <- which(pos_seq %% 5 == 0)
    col_ticks[idx] <- pos_seq[idx]
  }

  aa_tag  <- if (!is.null(aa_filter)) paste0("_", paste(aa_filter, collapse = "-")) else ""
  top_tag <- if (!is.null(top_n))     paste0("_top", top_n) else ""
  tag     <- paste0("tRF-", end_type, "_", frag_type, aa_tag, top_tag, "_", tissue)

  base_args <- list(
    cluster_cols = FALSE, cluster_rows = FALSE,
    cellwidth    = 10, cellheight = 10, fontsize_row = 10,
    labels_col   = col_ticks, angle_col = "0", na_col = "lightgrey"
  )

  do.call(pheatmap, c(list(
    mat      = depth_mat,
    breaks   = seq(0, 12, length.out = 101),
    filename = file.path(output_dir, paste0(tag, "_coverage.png"))
  ), base_args))

  out_del <- file.path(output_dir, paste0("fig1G_", tag, "_deletion_rates.png"))
  do.call(pheatmap, c(list(
    mat      = dels_mat,
    breaks   = seq(0, 1, length.out = 101),
    color    = colors_custom,
    filename = out_del
  ), base_args))

  write.csv(depth_mat, file.path(output_dir, paste0(tag, "_coverage.csv")))
  write.csv(dels_mat,  file.path(output_dir, paste0("fig1G_", tag, "_deletion_rates.csv")))
  message("Saved: ", out_del, " + CSVs")
}

#------------------------------------------------------------------------------
# Generate Figure 1G using parameters from config.R
#------------------------------------------------------------------------------

plot_trf_heatmaps(
  df        = means_df,
  tissue    = heatmap_tissue,
  end_type  = heatmap_end,
  frag_type = heatmap_type,
  aa_filter = heatmap_aa,
  top_n     = heatmap_top_n
)
