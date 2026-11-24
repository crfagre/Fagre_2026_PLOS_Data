# Figure S3C: tRF coverage heatmap — 5p fragments, mouse cerebellum NT
# Mouse equivalent of figS2C; uses helpers from helpers.R

source("config.R")

library(pheatmap)
library(RColorBrewer)

# =============================================================================
# Configuration
# =============================================================================

TRF_TYPE      <- "5p"
CONDITION     <- "NT"
TISSUE        <- "Cerebellum"
MIN_DEPTH_TRF <- trf_min_depth

col_depth <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)

# =============================================================================
# Load tRF pileup CSVs for Cerebellum NT samples (mouse)
# =============================================================================

nt_samples <- read.csv(file.path(data_dir, samples_mm_csv)) %>%
  filter(Treatment == CONDITION, Group == TISSUE)

trfs_all <- bind_rows(lapply(seq_len(nrow(nt_samples)), function(i) {
  row  <- nt_samples[i, ]
  path <- file.path(trf_mm_dir, paste0(row$Sample, ".csv"))
  df   <- load_trfs(path)
  if (is.null(df)) return(NULL)
  df %>% mutate(sample_id = row$Sample)
}))

if (nrow(trfs_all) == 0) stop("No tRF data loaded from ", trf_mm_dir)

# All-base filter for depth/coverage visualisation
trfs_depth_all <- trfs_all %>%
  filter(depth >= MIN_DEPTH_TRF) %>%
  group_by(sample_id) %>%
  mutate(cpm = depth / sum(depth, na.rm = TRUE) * 1e6) %>%
  ungroup()

# =============================================================================
# Build depth matrix (5p, one row per isoacceptor family)
# =============================================================================

sub_d <- trfs_depth_all %>% filter(trf_type == TRF_TYPE)
if (nrow(sub_d) < 5) stop("Too few 5p depth rows after filtering.")

mat_d_abs <- make_trf_matrix_ref(sub_d, "cpm")
if (is.null(mat_d_abs) || nrow(mat_d_abs) < 2)
  stop("Could not build depth matrix.")

mat_log_abs <- log2(mat_d_abs + 1)

isoacc_per_row <- parse_trfs_family(rownames(mat_log_abs))

sort_col  <- "13"
if (sort_col %in% colnames(mat_log_abs)) {
  sort_vals <- mat_log_abs[, sort_col]; sort_vals[is.na(sort_vals)] <- 0
} else {
  message("Sort column 13 absent — falling back to row means")
  sort_vals <- rowMeans(mat_log_abs, na.rm = TRUE)
}
names(sort_vals) <- rownames(mat_log_abs)

all_df <- data.frame(ref_row = rownames(mat_log_abs), ind_depth = sort_vals,
                     isoacc = isoacc_per_row, stringsAsFactors = FALSE)
top_df           <- all_df[order(-all_df$ind_depth), ]
top_df           <- top_df[!duplicated(top_df$isoacc), ]
grp_avg          <- tapply(top_df$ind_depth, top_df$isoacc, mean, na.rm = TRUE)
top_df$grp_depth <- grp_avg[top_df$isoacc]
top_df           <- top_df[order(-top_df$grp_depth, -top_df$ind_depth), ]

mat_log_abs <- mat_log_abs[top_df$ref_row, , drop = FALSE]
rownames(mat_log_abs) <- sub("^tRNA-", "", top_df$isoacc)

# =============================================================================
# Plot and save
# =============================================================================

n_rows <- nrow(mat_log_abs)
n_cols <- ncol(mat_log_abs)
out_pdf <- file.path(output_dir, "figS3C_trf_coverage_heatmap_mm_5p.pdf")
out_csv <- file.path(output_dir, "figS3C_trf_coverage_heatmap_mm_5p_data.csv")
write.csv(as.data.frame(mat_log_abs), out_csv)
message("Saved: ", basename(out_csv))

pdf(out_pdf, width = max(4, n_cols * 0.18 + 3), height = max(4, n_rows * 0.18 + 2))
pheatmap(mat_log_abs,
         color        = col_depth,
         cluster_rows = FALSE, cluster_cols = FALSE,
         na_col       = "grey90",
         cellheight   = 8, cellwidth = 10,
         main         = paste0(TISSUE, " ", CONDITION, " (5p)  log2(CPM+1) depth"),
         fontsize = 10, fontsize_row = 7)
dev.off()
message("Saved: ", out_pdf)
