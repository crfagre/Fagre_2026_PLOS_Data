# Figure S4C: Clustered heatmap of Pearson correlation (r) of mapped miRNAs across mouse tissue replicates
# Uses samples_mm_csv (NT only) and sample_data_dir

source("config.R")

library(pheatmap)
library(purrr)

# =============================================================================
# Configuration
# =============================================================================

IDXSTATS_FILE <- "idxstats.txt"   # contains mmu- miRNA names

# =============================================================================
# Load idxstats.txt for all NT mouse samples
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
if (length(counts_list) == 0) stop("No idxstats.txt files found in ", sample_data_dir)

raw_counts <- reduce(counts_list, function(a, b) merge(a, b, by = "ref", all = TRUE))
raw_counts[is.na(raw_counts)] <- 0

# =============================================================================
# Filter for mmu- miRNA refs
# =============================================================================

mirna_df <- raw_counts[grepl("^mmu", raw_counts$ref), ]
rownames(mirna_df) <- mirna_df[, 1]
mirna_df <- mirna_df[, -1]

# =============================================================================
# CPM normalize, filter zeros, log2 transform, Pearson correlation
# =============================================================================

# Remove refs with zero counts in any sample (all(x > 0))
nonzero_filter <- apply(mirna_df, 1, function(x) all(x > 0))
mirna_df       <- mirna_df[nonzero_filter, , drop = FALSE]

lib_sizes    <- colSums(mirna_df)
cpm_matrix   <- sweep(mirna_df, 2, lib_sizes, FUN = "/") * 1e6
log2_cpm_mat <- log2(cpm_matrix)

cor_matrix <- cor(log2_cpm_mat, method = "pearson")

# =============================================================================
# Plot and save
# =============================================================================

out_png <- file.path(output_dir, "figS4C_mirna_correlation_heatmap.png")
out_csv <- file.path(output_dir, "figS4C_mirna_correlation_heatmap_data.csv")
write.csv(as.data.frame(cor_matrix), out_csv)
message("Saved: ", basename(out_csv))

pheatmap(
  cor_matrix,
  display_numbers = TRUE,
  number_color    = "black",
  fontsize_number = 8,
  color           = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  main            = "Sample-to-Sample Pearson Correlation (miRNA Expression)",
  filename        = out_png,
  width           = 7,
  height          = 6
)

message("Saved: ", out_png)
