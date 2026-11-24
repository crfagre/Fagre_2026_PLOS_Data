# Figure S3E: tRF deletion rate heatmap — 5p fragments, mouse cerebellum RBS
# Mouse equivalent (5p) of figS2E; uses helpers from helpers.R

source("config.R")

library(pheatmap)

# =============================================================================
# Configuration
# =============================================================================

TRF_TYPE      <- "5p"
CONDITION     <- "RBS"
TISSUE        <- "Cerebellum"
MIN_DEPTH_TRF <- trf_min_depth

col_del <- colorRampPalette(c("white", "#c0392b"))(100)

# =============================================================================
# Load tRF pileup CSVs for Cerebellum RBS samples (mouse)
# =============================================================================

bs_samples <- read.csv(file.path(data_dir, samples_mm_csv)) %>%
  filter(Treatment == CONDITION, Group == TISSUE)

trfs_all <- bind_rows(lapply(seq_len(nrow(bs_samples)), function(i) {
  row  <- bs_samples[i, ]
  path <- file.path(trf_mm_dir, paste0(row$Sample, ".csv"))
  df   <- load_trfs(path)
  if (is.null(df)) return(NULL)
  df %>% mutate(sample_id = row$Sample)
}))

if (nrow(trfs_all) == 0) stop("No tRF data loaded from ", trf_mm_dir)

# T-only filter for del_rate analysis
trfs_filt <- trfs_all %>%
  filter(toupper(.data$base) == "T", depth >= MIN_DEPTH_TRF) %>%
  group_by(sample_id) %>%
  mutate(cpm = depth / sum(depth, na.rm = TRUE) * 1e6) %>%
  ungroup()

# =============================================================================
# Build Sprinzl coordinate mapping
# =============================================================================

sprinzl_res   <- build_sprinzl_lookup(hs_modomics_tsv, unique(trfs_filt$trf_family))
sprinzl_lut   <- sprinzl_res$lookup
SPRINZL_ORDER <- sprinzl_res$sprinzl_order

if (is.null(sprinzl_lut) || nrow(sprinzl_lut) == 0) {
  message("Sprinzl mapping unavailable — using raw fragment positions.")
  trfs_filt$sprinzl <- as.character(trfs_filt$start)
  SPRINZL_ORDER     <- as.character(sort(unique(trfs_filt$start)))
} else {
  message("Sprinzl map built: ", length(unique(sprinzl_lut$trf_family)), " tRF families mapped.")
  trfs_filt <- trfs_filt %>%
    left_join(sprinzl_lut, by = c("trf_family", "trf_type", "start"))
}

# =============================================================================
# Build deletion rate matrix (5p, Sprinzl coordinates)
# =============================================================================

sub_r <- trfs_filt %>% filter(trf_type == TRF_TYPE)
if (nrow(sub_r) < 5) stop("Too few 5p T-position rows after filtering.")

mat_r <- make_trf_matrix(sub_r, "del_rate", SPRINZL_ORDER)
if (is.null(mat_r) || nrow(mat_r) < 2)
  stop("Could not build deletion rate matrix.")

# =============================================================================
# Plot and save
# =============================================================================

out_pdf <- file.path(output_dir, "figS3E_trf_deletion_heatmap_mm_5p.pdf")

pdf(out_pdf, width = 14, height = max(6, nrow(mat_r) * 0.28))
write.csv(mat_r, sub(".pdf", ".csv", out_pdf), na = "")
pheatmap(mat_r,
         color         = col_del,
         breaks        = seq(0, 1, length.out = 101),
         cluster_rows  = FALSE, cluster_cols = FALSE,
         na_col        = "grey90",
         main          = paste0(TISSUE, " ", CONDITION, " (5p)  avg deletion rate"),
         fontsize = 8, fontsize_row = 7)
dev.off()
message("Saved: ", out_pdf)
