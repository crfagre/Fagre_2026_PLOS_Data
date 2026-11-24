# Figure S3A: Clustered heatmap of Pearson correlation (r) of tRNA fragments across mouse tissues
# Input: idxstats_trna.txt per NT mouse sample; uses clean_labels from helpers.R

source("config.R")

library(pheatmap)
library(purrr)

# =============================================================================
# Configuration
# =============================================================================

IDXSTATS_FILE <- "idxstats_trna.txt"
CPM_FILTER    <- 50     # matching R_plots_fig2.Rmd: rowSums > 50

# =============================================================================
# Load idxstats_trna.txt for all NT mouse samples
# =============================================================================

nt_samples <- read.csv(file.path(data_dir, samples_mm_csv)) %>%
  filter(Treatment == "NT")

counts_list <- lapply(seq_len(nrow(nt_samples)), function(i) {
  sample <- nt_samples$Sample[i]
  path   <- file.path(sample_data_dir, sample, IDXSTATS_FILE)
  if (!file.exists(path)) {
    message("MISSING: ", path)
    return(NULL)
  }
  read.table(path, header = FALSE, sep = "\t",
             col.names = c("ref", "length", "aligned", "unaligned")) %>%
    select(ref, aligned) %>%
    setNames(c("ref", sample))
})

counts_list <- Filter(Negate(is.null), counts_list)
if (length(counts_list) == 0) stop("No idxstats_trna.txt files found in ", sample_data_dir)

raw_counts <- reduce(counts_list, function(a, b) merge(a, b, by = "ref", all = TRUE))
raw_counts[is.na(raw_counts)] <- 0

# =============================================================================
# Clean, filter, and collapse by simplified_ref + End_Type
# =============================================================================

expr_df <- raw_counts %>%
  filter(
    grepl("tRNA-", ref),
    !grepl("^mt-", ref),
    !grepl("mito", ref, ignore.case = TRUE),
    ref != "*"
  ) %>%
  mutate(
    End_Type       = case_when(
      grepl("3p$", ref) ~ "3p",
      grepl("5p$", ref) ~ "5p",
      TRUE              ~ "other"
    ),
    simplified_ref = clean_labels(ref)
  ) %>%
  filter(End_Type %in% c("3p", "5p")) %>%
  select(-ref) %>%
  group_by(simplified_ref, End_Type) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")

# =============================================================================
# CPM normalize, filter, log2 transform, and compute Pearson correlation
# =============================================================================

corr_data <- expr_df %>%
  mutate(feature_id = paste(simplified_ref, End_Type, sep = "_")) %>%
  select(-simplified_ref, -End_Type) %>%
  column_to_rownames("feature_id")

mat      <- as.matrix(corr_data)
lib_sizes <- colSums(mat)
cpm_mat  <- sweep(mat, 2, lib_sizes, FUN = "/") * 1e6
cpm_mat  <- cpm_mat[rowSums(cpm_mat) > CPM_FILTER, , drop = FALSE]
log_cpm  <- log2(cpm_mat + 1)

cor_matrix <- cor(log_cpm, method = "pearson")

# =============================================================================
# Plot and save
# =============================================================================

out_png <- file.path(output_dir, "figS3A_trf_correlation_heatmap.png")
out_csv <- file.path(output_dir, "figS3A_trf_correlation_heatmap_data.csv")
write.csv(as.data.frame(cor_matrix), out_csv)
message("Saved: ", basename(out_csv))

pheatmap(
  cor_matrix,
  display_numbers  = TRUE,
  number_color     = "black",
  fontsize_number  = 8,
  color            = colorRampPalette(c("white", "lightcoral", "darkred"))(100),
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,
  main             = "Sample-to-Sample Pearson Correlation (tRF Expression)",
  filename         = out_png,
  width            = 7,
  height           = 6
)

message("Saved: ", out_png)
