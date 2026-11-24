# Figure 1B: Correlation of untreated vs. bisulfite-treated sample counts
#
# Edit config.R to set nt_counts_file, bs_counts_file, counts_prefix.

source("config.R")

library(dplyr)
library(readr)
library(ggplot2)

#------------------------------------------------------------------------------
# Load count tables and compute log2 CPM
#------------------------------------------------------------------------------

nt <- read.table(nt_counts_file)
bs <- read.table(bs_counts_file)

nt <- nt[grepl(counts_prefix, nt[, 1]), ]
bs <- bs[grepl(counts_prefix, bs[, 1]), ]

plot_df <- merge(nt[, c(1, 3)], bs[, c(1, 3)], by = 1) %>%
  setNames(c("ref", "NT", "BS")) %>%
  as.data.frame() %>%
  mutate(
    logNT = log2(1e6 * NT / sum(NT)),
    logBS = log2(1e6 * BS / sum(BS))
  )

r_val <- round(cor(plot_df$logNT, plot_df$logBS, use = "complete.obs"), 3)

#------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------

p <- ggplot(plot_df, aes(x = logNT, y = logBS)) +
  geom_point(aes(color = logNT > log2(100) & logBS > log2(100)), alpha = 0.7) +
  geom_vline(xintercept = log2(100), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = log2(100), linetype = "dashed", color = "grey50") +
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

out_file <- file.path(output_dir, "fig1B_correlation.png")
ggsave(out_file, plot = p, width = 6, height = 5, dpi = 300, bg = "white")
write_csv(plot_df, file.path(output_dir, "fig1B_correlation.csv"))
message("Saved: ", out_file, " + fig1B_correlation.csv")
print(p)
