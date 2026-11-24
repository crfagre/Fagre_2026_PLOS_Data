# Figure 4B/C: DESeq2 pipeline → CDF of Log2FC binned by TargetScan site type
#
# Inputs:  <deseq_data_dir>/<fig4_conditions_file>        — sample sheet (ID + Sample cols)
#          <deseq_data_dir>/<ID>/<fig4_counts_file>       — HTSeq count files
#          <ref_dir>/<fig4_gene_site_map_file>            — pre-built RefSeq → gene + site map
#            Columns: transcript_id, gene_name, miRNA_site_category, miRNA_site_simple
# Outputs per comparison (tag = sanitised label):
#   fig4_deseq_<tag>.csv           — DESeq2 results with gene + site annotation
#   fig4_mirna_cdf_<tag>.png       — CDF plot with zoomed inset
#   fig4_mirna_cdf_<tag>_data.csv  — per-gene LFC + site annotation used for plot
#   fig4_mirna_cdf_<tag>_ks.csv    — KS test results vs "None"
# Shared output:
#   fig4_normalized_counts.csv     — DESeq2 size-factor-normalized counts (all samples)
#
# Edit config.R: fig4_conditions_file, fig4_counts_file, fig4_gene_site_map_file,
#   fig4_ref_condition, fig4_min_counts, fig4_min_samples, fig4_comparisons

source("config.R")

library(DESeq2)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

# =============================================================================
# 1. Load sample sheet and build DESeq2 object
# =============================================================================

conditions_path <- file.path(deseq_data_dir, fig4_conditions_file)
if (!file.exists(conditions_path))
  stop("Conditions file not found: ", conditions_path)

conditions <- read.table(conditions_path, header = TRUE)

for (col in c("ID", "Sample")) {
  if (!col %in% colnames(conditions))
    stop("Column '", col, "' not found in ", fig4_conditions_file,
         ".\n  Available: ", paste(colnames(conditions), collapse = ", "))
}

sampleTable <- data.frame(
  sampleName       = conditions[["ID"]],
  sampleFiles      = file.path(conditions[["ID"]], fig4_counts_file),
  sampleCondition  = conditions[["Sample"]],
  stringsAsFactors = FALSE
)

missing_files <- sampleTable$sampleFiles[
  !file.exists(file.path(deseq_data_dir, sampleTable$sampleFiles))]
if (length(missing_files) > 0)
  stop("Missing count file(s):\n",
       paste(file.path(deseq_data_dir, missing_files), collapse = "\n"))

message("Building DESeq2 object from ", nrow(sampleTable), " samples...")

dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable,
  directory   = deseq_data_dir,
  design      = ~sampleCondition
)

keep <- rowSums(counts(dds) >= fig4_min_counts) >= fig4_min_samples
dds  <- dds[keep, ]
message(sum(keep), " genes retained after count filter (>= ", fig4_min_counts,
        " counts in >= ", fig4_min_samples, " samples)")

dds$sampleCondition <- relevel(dds$sampleCondition, ref = fig4_ref_condition)
dds <- DESeq(dds)
message("DESeq2 complete.")

# Export normalized counts (once, shared across all comparisons)
norm_df <- as.data.frame(counts(dds, normalized = TRUE))
write_csv(cbind(transcript_id = rownames(norm_df), norm_df),
          file.path(output_dir, "fig4_normalized_counts.csv"))
message("Saved: fig4_normalized_counts.csv")

# =============================================================================
# 2. Load pre-built gene + site annotation map
# =============================================================================

map_path <- file.path(ref_dir, fig4_gene_site_map_file)
if (!file.exists(map_path))
  stop("Gene-site map not found: ", map_path,
       "\n  Build it once from biomaRt + TargetScan and save with columns: ",
       "transcript_id, gene_name, miRNA_site_category, miRNA_site_simple")

gene_site_map <- read.csv(map_path, na.strings = c("", "NA"))

for (col in c("transcript_id", "gene_name", "miRNA_site_category", "miRNA_site_simple")) {
  if (!col %in% colnames(gene_site_map))
    stop("Column '", col, "' not found in ", fig4_gene_site_map_file)
}

message(nrow(gene_site_map), " entries loaded from gene-site map.")

# =============================================================================
# 3. Per-comparison: extract results → annotate → CDF + exports
# =============================================================================

run_comparison <- function(label, contrast) {

  message("\n=== Processing: ", label, " ===")

  res    <- results(dds, contrast = contrast)
  res_df <- as.data.frame(res)

  # Annotate from pre-built map
  res_df$gene_name           <- gene_site_map$gene_name[
    match(rownames(res_df), gene_site_map$transcript_id)]
  res_df$miRNA_site_category <- gene_site_map$miRNA_site_category[
    match(rownames(res_df), gene_site_map$transcript_id)]
  res_df$miRNA_site_simple   <- gene_site_map$miRNA_site_simple[
    match(rownames(res_df), gene_site_map$transcript_id)]
  res_df$miRNA_site_simple[is.na(res_df$miRNA_site_simple)] <- "None"

  tag       <- gsub("[^A-Za-z0-9]", "_", label)
  out_deseq <- file.path(output_dir, paste0("fig4_deseq_", tag, ".csv"))
  write_csv(cbind(transcript_id = rownames(res_df), res_df), out_deseq)
  message("Saved DESeq results: ", basename(out_deseq))

  fc_df_clean <- res_df %>%
    filter(!is.na(log2FoldChange), !is.na(miRNA_site_simple))

  if (nrow(fc_df_clean) == 0) {
    warning("No rows with valid LFC for: ", label, " — skipping plots.")
    return(invisible(NULL))
  }

  # Shared color scale and ECDF layers (reused in main plot and inset)
  site_colors <- c("None"    = "grey",
                   "8mer"    = "red",
                   "7mer-M8" = "green",
                   "7mer-A1" = "blue",
                   "6mer"    = "orange")

  ecdf_layers <- list(
    stat_ecdf(data = subset(fc_df_clean, miRNA_site_simple != "None"),
              linewidth = 1.2, alpha = 0.8),
    stat_ecdf(data = subset(fc_df_clean, miRNA_site_simple == "None"),
              linewidth = 1.2, alpha = 0.8),
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6),
    scale_color_manual(values = site_colors, breaks = names(site_colors))
  )

  # Main plot — zoomed to ±0.5; ECDF computed from full data
  cdf_plot <- ggplot(fc_df_clean, aes(x = log2FoldChange, color = miRNA_site_simple)) +
    ecdf_layers +
    coord_cartesian(xlim = c(-0.5, 0.5)) +
    labs(x     = bquote(Log[2] ~ "(Fold-change)"),
         y     = "Cumulative Probability",
         color = "Site type") +
    theme_minimal() +
    theme(legend.title      = element_text(size = 18),
          legend.text       = element_text(size = 16),
          legend.key.height = unit(1.8, "lines"),
          axis.title        = element_text(size = 16),
          axis.text         = element_text(size = 14))

  # Inset — zoomed to crossover region
  cdf_inset <- ggplot(fc_df_clean, aes(x = log2FoldChange, color = miRNA_site_simple)) +
    ecdf_layers +
    coord_cartesian(xlim = c(-0.025, 0.025), ylim = c(0.4, 0.6)) +
    scale_x_continuous(breaks = c(-0.025, 0.025)) +
    scale_y_continuous(breaks = c(0.4, 0.6)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(legend.position  = "none",
          panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
          axis.text        = element_text(size = 10),
          panel.grid.minor = element_blank())

  cdf_composite <- cdf_plot +
    inset_element(cdf_inset, left = 0.02, bottom = 0.54, right = 0.40, top = 0.98)

  # KS tests: each site bin vs "None"
  none_values <- fc_df_clean$log2FoldChange[fc_df_clean$miRNA_site_simple == "None"]

  ks_results <- bind_rows(lapply(c("6mer", "7mer-A1", "7mer-M8", "8mer"), function(cat) {
    cat_values <- fc_df_clean$log2FoldChange[fc_df_clean$miRNA_site_simple == cat]
    if (length(cat_values) == 0)
      return(data.frame(site = cat, n = 0L, ks_D = NA_real_, p_value = NA_real_))
    ks <- ks.test(none_values, cat_values)
    data.frame(site = cat, n = length(cat_values), ks_D = ks$statistic, p_value = ks$p.value)
  }))

  message("KS tests (None vs site bin):")
  for (i in seq_len(nrow(ks_results))) {
    message(sprintf("  None vs %-8s  n = %4d  p = %s  D = %.3f",
                    ks_results$site[i], ks_results$n[i],
                    format.pval(ks_results$p_value[i], digits = 3),
                    ks_results$ks_D[i]))
  }

  fc_df_clean %>%
    group_by(miRNA_site_simple) %>%
    summarise(n             = n(),
              median_log2FC = round(median(log2FoldChange), 3),
              mean_log2FC   = round(mean(log2FoldChange),   3),
              q25           = round(quantile(log2FoldChange, 0.25), 3),
              q75           = round(quantile(log2FoldChange, 0.75), 3),
              .groups = "drop") %>%
    print()

  out_png  <- file.path(output_dir, paste0("fig4_mirna_cdf_", tag, ".png"))
  out_data <- file.path(output_dir, paste0("fig4_mirna_cdf_", tag, "_data.csv"))
  out_ks   <- file.path(output_dir, paste0("fig4_mirna_cdf_", tag, "_ks.csv"))

  ggsave(out_png, plot = cdf_composite, width = 8, height = 6, dpi = 300, bg = "white")
  write_csv(cbind(transcript_id = rownames(fc_df_clean), fc_df_clean), out_data)
  write_csv(ks_results, out_ks)
  message("Saved: ", basename(out_png), " + data CSV + KS CSV")
  print(cdf_composite)
}

for (label in names(fig4_comparisons)) {
  run_comparison(label, fig4_comparisons[[label]])
}
