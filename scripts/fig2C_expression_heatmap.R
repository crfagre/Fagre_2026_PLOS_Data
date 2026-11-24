# Figure 2C: Clustered heatmap of normalised tRF abundance across tissues
#
# Reads idxstats files directly (no pipeline dependency).
# Inputs:  <sample_data_dir>/<sample_id>/idxstats_trna.txt
#          <data_dir>/samples.csv
# Outputs: PNG heatmaps + CSV of underlying log2-CPM matrix
#
# Edit config.R: data_dir, sample_data_folder, fig2c_cpm_filter

source("config.R")

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)
library(purrr)
library(ComplexHeatmap)

# -----------------------------------------------------------------------------
# Load idxstats count tables
# -----------------------------------------------------------------------------

samples_meta <- read_csv(file.path(data_dir, samples_csv), show_col_types = FALSE) %>%
  filter(Treatment == "NT")

counts_list <- lapply(1:nrow(samples_meta), function(i) {
  sid <- samples_meta$Sample[i]
  p   <- file.path(sample_data_dir, sid, "idxstats_trna.txt")
  if (!file.exists(p)) {
    message("MISSING idxstats for ", sid, ". Checked: ", p)
    return(NULL)
  }
  read_tsv(p, col_names = c("ref", "length", "aligned", "unaligned"),
           show_col_types = FALSE) %>%
    select(ref, aligned) %>%
    rename(!!sid := aligned)
})

counts_list <- compact(counts_list)
if (length(counts_list) == 0) stop("No idxstats files found.")

raw_counts <- reduce(counts_list, full_join, by = "ref") %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

# -----------------------------------------------------------------------------
# Clean, collapse grouped variants, filter mito
# -----------------------------------------------------------------------------

expr_df <- raw_counts %>%
  filter(str_detect(ref, "tRNA-"),
         !str_detect(ref, "^mt-"),
         !str_detect(ref, "mito")) %>%
  mutate(
    End_Type       = if_else(str_detect(ref, "3p$"), "3p",
                     if_else(str_detect(ref, "5p$"), "5p", "other")),
    simplified_ref = clean_labels(ref)
  ) %>%
  filter(End_Type %in% c("3p", "5p")) %>%
  select(-ref) %>%
  group_by(simplified_ref, End_Type) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop")

# -----------------------------------------------------------------------------
# Heatmap function: CPM normalise, log2 transform, cluster, plot
# Returns log2-CPM matrix invisibly; writes PNG and CSV.
# -----------------------------------------------------------------------------

plot_expression_heatmap <- function(df, target_end) {

  mat_df <- df %>%
    filter(End_Type == target_end) %>%
    select(-End_Type) %>%
    column_to_rownames("simplified_ref")

  if (nrow(mat_df) == 0) return(invisible(NULL))

  mat     <- as.matrix(mat_df)
  cpm_mat <- sweep(mat, 2, colSums(mat), FUN = "/") * 1e6
  cpm_mat <- cpm_mat[rowSums(cpm_mat) > fig2c_cpm_filter, , drop = FALSE]

  if (nrow(cpm_mat) == 0) {
    message("No ", target_end, " tRFs passed CPM filter > ", fig2c_cpm_filter)
    return(invisible(NULL))
  }

  plot_mat <- log2(cpm_mat + 1)   # tRFs as rows, samples as columns

  # Order columns by group according to sample_order
  if (!is.null(sample_order)) {
    grp_map    <- setNames(samples_meta$Group, samples_meta$Sample)
    col_order  <- order(match(grp_map[colnames(plot_mat)], sample_order))
    plot_mat   <- plot_mat[, col_order, drop = FALSE]
    do_cluster_cols <- FALSE
  } else {
    do_cluster_cols <- TRUE
  }

  plot_width  <- 3
  plot_height <- 5

  out_png <- file.path(output_dir, paste0("fig2C_Expression_", target_end, ".png"))
  out_csv <- sub("\\.png$", ".csv", out_png)

  png(out_png, width = plot_width, height = plot_height, units = "in", res = 300)
  p <- ComplexHeatmap::pheatmap(
    plot_mat,
    show_rownames  = FALSE,
    show_colnames  = TRUE,
    cluster_cols   = do_cluster_cols,
    cluster_rows   = TRUE,
    main           = paste0(target_end, " tRF Expression (Log2 CPM)"),
    angle_col      = "45"
  )
  print(p)
  dev.off()

  write.csv(plot_mat, out_csv)
  message("Saved: ", out_png, " + CSV")

  invisible(plot_mat)
}

# -----------------------------------------------------------------------------
# Generate Figure 2D
# -----------------------------------------------------------------------------

plot_expression_heatmap(expr_df, "5p")
plot_expression_heatmap(expr_df, "3p")
