# Figure 4B/C: Cumulative distribution of Log2FC binned by TargetScan site type
#              One CDF panel per comparison defined in fig4_comparisons.
#              KS test p-values compare each site bin to transcripts with no site.
#
# Inputs:  annotated DESeq2 CSVs in <data_dir> — output of CF153_analysis.Rmd
#            save_annotated chunk; must contain log2FoldChange and
#            miRNA_site_category columns
# Outputs: per-comparison PNG, data CSV, and KS CSV in <output_dir>
#
# Edit config.R: fig4_comparisons, fig4_lfc_col, fig4_site_category_col

source("config.R")

library(dplyr)
library(readr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Helper: run one CDF panel and export
# -----------------------------------------------------------------------------

run_cdf <- function(label, file, lfc_col, site_cat_col) {

  fc_df <- read_csv(file.path(data_dir, file), show_col_types = FALSE)

  for (col in c(lfc_col, site_cat_col)) {
    if (!col %in% colnames(fc_df))
      stop("Column '", col, "' not found in ", file,
           ".\n  Available: ", paste(colnames(fc_df), collapse = ", "))
  }

  # Simplify to strongest site type (8mer > 7mer-M8 > 7mer-A1 > 6mer > None)
  fc_df <- fc_df %>%
    mutate(miRNA_site_simple = case_when(
      is.na(.data[[site_cat_col]])                                        ~ "None",
      grepl("8mer",    .data[[site_cat_col]], ignore.case = TRUE)        ~ "8mer",
      grepl("7mer-M8", .data[[site_cat_col]], ignore.case = TRUE)        ~ "7mer-M8",
      grepl("7mer-A1", .data[[site_cat_col]], ignore.case = TRUE)        ~ "7mer-A1",
      grepl("6mer",    .data[[site_cat_col]], ignore.case = TRUE)        ~ "6mer",
      TRUE                                                                 ~ "None"
    ))

  fc_df_clean <- fc_df %>%
    filter(!is.na(.data[[lfc_col]]), !is.na(miRNA_site_simple)) %>%
    rename(log2FoldChange = !!lfc_col)

  # CDF plot
  cdf_plot <- ggplot(fc_df_clean, aes(x = log2FoldChange, color = miRNA_site_simple)) +
    stat_ecdf(data = subset(fc_df_clean, miRNA_site_simple != "None"),
              linewidth = 1.2, alpha = 0.8) +
    stat_ecdf(data = subset(fc_df_clean, miRNA_site_simple == "None"),
              linewidth = 1.2, alpha = 0.8) +
    scale_color_manual(
      values = c("None"    = "grey",
                 "7mer-A1" = "blue",
                 "7mer-M8" = "green",
                 "8mer"    = "red",
                 "6mer"    = "orange"),
      breaks = c("None", "7mer-A1", "7mer-M8", "8mer", "6mer")
    ) +
    labs(title = paste("CDF of Log2FC —", label),
         x     = "Log2 Fold Change",
         y     = "Cumulative Probability",
         color = "Site type") +
    theme_minimal() +
    theme(plot.title   = element_text(size = 14, hjust = 0.5),
          legend.title = element_text(size = 18),
          legend.text  = element_text(size = 16),
          axis.title   = element_text(size = 16),
          axis.text    = element_text(size = 14)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6)

  # KS tests vs "None"
  none_values <- fc_df_clean$log2FoldChange[fc_df_clean$miRNA_site_simple == "None"]

  ks_results <- bind_rows(lapply(c("6mer", "7mer-A1", "7mer-M8", "8mer"), function(cat) {
    cat_values <- fc_df_clean$log2FoldChange[fc_df_clean$miRNA_site_simple == cat]
    if (length(cat_values) == 0)
      return(data.frame(site = cat, n = 0L, ks_D = NA_real_, p_value = NA_real_))
    ks <- ks.test(none_values, cat_values)
    data.frame(site = cat, n = length(cat_values), ks_D = ks$statistic, p_value = ks$p.value)
  }))

  message("\n--- KS tests: ", label, " ---")
  for (i in seq_len(nrow(ks_results))) {
    cat(sprintf("None vs %-8s  p = %s  D = %.3f\n",
                ks_results$site[i],
                format.pval(ks_results$p_value[i], digits = 3),
                ks_results$ks_D[i]))
  }

  # Summary statistics
  summary_stats <- fc_df_clean %>%
    group_by(miRNA_site_simple) %>%
    summarise(n             = n(),
              median_log2FC = round(median(log2FoldChange), 3),
              mean_log2FC   = round(mean(log2FoldChange),   3),
              q25           = round(quantile(log2FoldChange, 0.25), 3),
              q75           = round(quantile(log2FoldChange, 0.75), 3),
              .groups = "drop")
  print(summary_stats)

  # Export
  tag      <- gsub("[^A-Za-z0-9]", "_", label)
  out_png  <- file.path(output_dir, paste0("fig4_mirna_cdf_", tag, ".png"))
  out_data <- file.path(output_dir, paste0("fig4_mirna_cdf_", tag, "_data.csv"))
  out_ks   <- file.path(output_dir, paste0("fig4_mirna_cdf_", tag, "_ks.csv"))

  ggsave(out_png, plot = cdf_plot, width = 8, height = 6, dpi = 300, bg = "white")
  write_csv(fc_df_clean, out_data)
  write_csv(ks_results,  out_ks)
  message("Saved: ", out_png, " + data CSV + KS CSV")
  print(cdf_plot)
}

# -----------------------------------------------------------------------------
# Run all comparisons defined in config.R
# -----------------------------------------------------------------------------

for (label in names(fig4_comparisons)) {
  run_cdf(label        = label,
          file         = fig4_comparisons[[label]],
          lfc_col      = fig4_lfc_col,
          site_cat_col = fig4_site_category_col)
}
