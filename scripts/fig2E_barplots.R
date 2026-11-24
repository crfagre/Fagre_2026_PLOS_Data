# Figure 2E: Barplots of deletion rates at specific tRF sites, one-way ANOVA + Tukey HSD
#
# Requires pipeline_site_calling.R to have been run first (call_all_groups = TRUE).
# Inputs:  <data_dir>/called_site_data.csv
# Outputs: PNG barplot + CSV of underlying data
#
# Edit config.R: data_dir, fig2e_ref, fig2e_pos, fig2e_end

source("config.R")

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Load pipeline output
# -----------------------------------------------------------------------------

called_site_data <- read_csv(file.path(pileup_dir, "called_site_data.csv"),
                              show_col_types = FALSE)

# -----------------------------------------------------------------------------
# Barplot + ANOVA/Tukey function
# Returns the site data invisibly; writes PNG and CSV.
# -----------------------------------------------------------------------------

plot_site_barplot <- function(df, target_ref, target_pos, target_end,
                               min_depth_filter = 50) {

  site_data <- df %>%
    filter(
      str_detect(ref, target_ref),
      str_detect(ref, paste0("_", target_end, "$")),
      pos == target_pos
    ) %>%
    filter(depth >= min_depth_filter)

  if (nrow(site_data) == 0) {
    message("No data found for ref='", target_ref, "', pos=", target_pos,
            ", end=", target_end)
    return(invisible(NULL))
  }

  # One-way ANOVA
  aov_fit <- aov(del_rate ~ Tissue, data = site_data)
  aov_p   <- summary(aov_fit)[[1]][["Pr(>F)"]][1]
  message("ANOVA p-value: ", signif(aov_p, 3))

  # Tukey HSD
  tukey       <- TukeyHSD(aov_fit)$Tissue
  tukey_df    <- as.data.frame(tukey) %>%
    tibble::rownames_to_column("comparison") %>%
    mutate(sig = case_when(
      `p adj` < 0.01 ~ "**",
      `p adj` < 0.05 ~ "*",
      TRUE            ~ "ns"
    ))

  # Summary stats per tissue
  summary_df <- site_data %>%
    group_by(Tissue) %>%
    summarize(mean_del = mean(del_rate, na.rm = TRUE),
              sem_del  = sd(del_rate, na.rm = TRUE) / sqrt(n()),
              .groups = "drop")

  tissue_levels <- if (!is.null(sample_order)) sample_order else sort(unique(site_data$Tissue))
  summary_df    <- summary_df %>% mutate(Tissue = factor(Tissue, levels = tissue_levels))
  site_data     <- site_data   %>% mutate(Tissue = factor(Tissue, levels = tissue_levels))

  # Build significance brackets for significant Tukey pairs
  sig_pairs <- tukey_df %>% filter(sig != "ns")
  brackets  <- list()
  if (nrow(sig_pairs) > 0) {
    y_max   <- max(site_data$del_rate, na.rm = TRUE)
    y_step  <- (1 - y_max) / (nrow(sig_pairs) + 1)
    for (k in seq_len(nrow(sig_pairs))) {
      parts <- strsplit(sig_pairs$comparison[k], "-")[[1]]
      brackets[[k]] <- list(
        x1 = parts[1], x2 = parts[2],
        y  = y_max + k * y_step,
        label = sig_pairs$sig[k]
      )
    }
  }

  p <- ggplot(summary_df, aes(x = Tissue, y = mean_del, fill = Tissue)) +
    geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean_del - sem_del, ymax = mean_del + sem_del),
                  width = 0.2) +
    geom_jitter(data = site_data, aes(x = Tissue, y = del_rate),
                width = 0.1, size = 1.5, alpha = 0.7, inherit.aes = FALSE) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.15))) +
    labs(
      title    = paste0(target_ref, " pos.", target_pos, " (", target_end, ")"),
      subtitle = paste0("One-way ANOVA p = ", signif(aov_p, 3)),
      x        = NULL,
      y        = "Deletion rate (BS)"
    ) +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(size = 11, face = "bold"))

  # Add significance brackets
  for (br in brackets) {
    x1_num <- which(levels(summary_df$Tissue) == br$x1)
    x2_num <- which(levels(summary_df$Tissue) == br$x2)
    p <- p +
      annotate("segment", x = x1_num, xend = x2_num, y = br$y, yend = br$y) +
      annotate("segment", x = x1_num, xend = x1_num, y = br$y, yend = br$y - 0.015) +
      annotate("segment", x = x2_num, xend = x2_num, y = br$y, yend = br$y - 0.015) +
      annotate("text", x = (x1_num + x2_num) / 2, y = br$y + 0.02,
               label = br$label, size = 4)
  }

  ref_clean <- str_replace_all(target_ref, "[^A-Za-z0-9]", "_")
  tag       <- paste0(ref_clean, "_pos", target_pos, "_", target_end)
  out_png   <- file.path(output_dir, paste0("fig2E_barplot_", tag, ".png"))
  out_csv   <- sub("\\.png$", ".csv", out_png)

  ggsave(out_png, plot = p, width = 4, height = 5, dpi = 300, bg = "white")
  write_csv(site_data, out_csv)
  write_csv(tukey_df,  sub("\\.csv$", "_tukey.csv", out_csv))
  message("Saved: ", out_png, " + data CSVs")

  invisible(site_data)
}

# -----------------------------------------------------------------------------
# Generate Figure 2F using parameters from config.R
# -----------------------------------------------------------------------------

plot_site_barplot(called_site_data,
                  target_ref = fig2e_ref,
                  target_pos = fig2e_pos,
                  target_end = fig2e_end)
