# Figure 3D: miRNA expression heatmap across mouse tissues

source("config.R")

library(pheatmap)

# =============================================================================
# Configuration
# =============================================================================

IDXSTATS_FILE <- fig3d_idxstats_file   # "idxstats.txt" — contains mmu- miRNA names
CPM_FILTER    <- fig3d_cpm_filter      # 200

# =============================================================================
# Build sample_list and sample_names from samples_mm.csv (NT only)
# =============================================================================

samples_mm   <- read.csv(file.path(data_dir, samples_mm_csv))
nt_samples   <- samples_mm[samples_mm$Treatment == "NT", ]
sample_list  <- file.path(sample_data_dir, nt_samples$Sample, IDXSTATS_FILE)
sample_names <- nt_samples$Sample

# =============================================================================
# =============================================================================

for(i in c(1:length(sample_list))){
  if(i==1){
    counts_df <- read.csv(sample_list[1], header = FALSE, sep = "\t",
                          col.names = c("ref","length","aligned","unaligned"))[,c("ref","aligned")]
    colnames(counts_df)[colnames(counts_df) == "aligned"] <- sample_names[1]
  } else {
    df_tmp <- read.csv(sample_list[i], header = FALSE, sep = "\t",
                       col.names = c("ref","length","aligned","unaligned"))[,c("ref","aligned")]
    colnames(df_tmp)[colnames(df_tmp) == "aligned"] <- sample_names[i]
    counts_df <- merge(counts_df, df_tmp, by=c("ref"))
  }
}

miRNA_df <- counts_df[grepl("^mmu", counts_df$ref), ]
rownames(miRNA_df) <- miRNA_df[, 1]
miRNA_df <- miRNA_df[, -1]

lib_sizes <- colSums(miRNA_df)
miRNA_df <- sweep(miRNA_df, 2, lib_sizes, FUN = "/") * 1e6

miRNA_df <- miRNA_df[rowSums(miRNA_df) > CPM_FILTER, ]

miRNA_df_log <- log2(miRNA_df)
miRNA_df_log[is.na(miRNA_df_log)] <- 0
miRNA_df_log[is.nan(as.matrix(miRNA_df_log))] <- 0
miRNA_df_log[is.infinite(as.matrix(miRNA_df_log))] <- 0

miRNA_df_log <- miRNA_df_log[rowSums(miRNA_df) != 0, ]

# =============================================================================
# Plot and save
# =============================================================================

out_png <- file.path(output_dir, "fig3D_mirna_expression_heatmap.png")
out_csv <- file.path(output_dir, "fig3D_mirna_expression_heatmap_data.csv")
write.csv(as.data.frame(miRNA_df_log), out_csv)
message("Saved: ", basename(out_csv))

png(out_png, width = 10, height = 8, units = "in", res = 300)
pheatmap(
  t(as.matrix(miRNA_df_log)),
  show_rownames = TRUE,
)
dev.off()
message("Saved: ", out_png)

pheatmap(
  t(as.matrix(miRNA_df_log)),
  show_rownames = TRUE,
)
