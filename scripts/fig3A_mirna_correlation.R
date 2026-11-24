# Figure 3A: Correlation of miRNA counts from untreated vs bisulfite-treated samples
#
# Edit config.R: fig3a_nt_counts_file, fig3a_bs_counts_file,
#                fig3a_cpm_threshold

source("config.R")

library(dplyr)
library(readr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Load count tables and compute log2 CPM
# -----------------------------------------------------------------------------

nt <- read.table(fig3a_nt_counts_file)
bs <- read.table(fig3a_bs_counts_file)

nt <- nt[grepl("^hsa-mir", nt[, 1]), ]
bs <- bs[grepl("^hsa-mir", bs[, 1]), ]

plot_df <- merge(nt[, c(1, 3)], bs[, c(1, 3)], by = 1) %>%
  setNames(c("ref", "NT", "BS")) %>%
  as.data.frame() %>%
  mutate(
    logNT = log2(1e6 * NT / sum(NT)),
    logBS = log2(1e6 * BS / sum(BS))
  )

r_val <- round(cor(plot_df$logNT, plot_df$logBS, use = "complete.obs"), 3)

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

cpm_thresh <- log2(fig3a_cpm_threshold)

p <- ggplot(plot_df, aes(x = logNT, y = logBS)) +
  geom_point(aes(color = logNT > cpm_thresh & logBS > cpm_thresh), alpha = 0.7) +
  geom_vline(xintercept = cpm_thresh, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = cpm_thresh, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("grey60", "steelblue"),
                     labels = c("Other", "High expression"),
                     name   = NULL) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = paste0("R = ", r_val), size = 4) +
  labs(
    x = expression("Untreated log"[2] * " CPM"),
    y = expression("Bisulfite log"[2] * " CPM")
  ) +
  theme_classic()

out_file <- file.path(output_dir, "fig3A_mirna_correlation.png")
ggsave(out_file, plot = p, width = 6, height = 5, dpi = 300, bg = "white")
write_csv(plot_df, file.path(output_dir, "fig3A_mirna_correlation.csv"))
message("Saved: ", out_file, " + fig3A_mirna_correlation.csv")
print(p)
