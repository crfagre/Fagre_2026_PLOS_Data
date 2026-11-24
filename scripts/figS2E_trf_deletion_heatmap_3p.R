# Figure S2E: tRF deletion rate heatmap — 3p fragments, human WT RBS
# Uses helpers: load_trfs, build_sprinzl_lookup, make_trf_matrix

source("config.R")

library(pheatmap)

# =============================================================================
# Configuration
# =============================================================================

TRF_TYPE      <- "3p"
CONDITION     <- "RBS"    # RBS condition for deletion rate heatmap
MIN_DEPTH_TRF <- trf_min_depth

col_del <- colorRampPalette(c("white", "#c0392b"))(100)

# =============================================================================
# Load tRF pileup CSVs for RBS samples (human)
# =============================================================================

rbs_samples <- read.csv(file.path(data_dir, samples_hs_csv)) %>%
  filter(Treatment == CONDITION)

trfs_all <- bind_rows(lapply(seq_len(nrow(rbs_samples)), function(i) {
  row  <- rbs_samples[i, ]
  path <- file.path(trf_hs_dir, paste0(row$Sample, ".csv"))
  df   <- load_trfs(path)
  if (is.null(df)) return(NULL)
  df %>% mutate(sample_id = row$Sample)
}))

if (nrow(trfs_all) == 0) stop("No tRF data loaded from ", trf_hs_dir)

# T-only filter for del_rate analysis (pseudouridine detection)
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
# Build deletion rate matrix (3p, Sprinzl coordinates)
# =============================================================================

sub_r <- trfs_filt %>% filter(trf_type == TRF_TYPE)
if (nrow(sub_r) < 5) stop("Too few 3p T-position rows after filtering.")

mat_r <- make_trf_matrix(sub_r, "del_rate", SPRINZL_ORDER)
if (is.null(mat_r) || nrow(mat_r) < 2)
  stop("Could not build deletion rate matrix for ", TRF_TYPE)

# =============================================================================
# Plot and save
# =============================================================================

out_pdf <- file.path(output_dir, "figS2E_trf_deletion_heatmap_3p.pdf")

pdf(out_pdf, width = 14, height = max(6, nrow(mat_r) * 0.28))
write.csv(mat_r, sub(".pdf", ".csv", out_pdf), na = "")
pheatmap(mat_r,
         color         = col_del,
         breaks        = seq(0, 1, length.out = 101),
         cluster_rows  = FALSE, cluster_cols = FALSE,
         na_col        = "grey90",
         main          = paste0("WT ", CONDITION, " (3p)  avg deletion rate"),
         fontsize = 8, fontsize_row = 7)
dev.off()
message("Saved: ", out_pdf)
