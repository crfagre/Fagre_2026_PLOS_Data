# Figure S4A: miRNA replicate correlation scatter — human WT NT (ncrna idxstats)
# miRNA equivalent of figS2A; uses load_idxstats from helpers.R

source("config.R")

library(dplyr)
library(readr)
library(ggplot2)

# =============================================================================
# Configuration
# =============================================================================

SAMPLE_SHEET <- file.path(data_dir, samples_hs_csv)
RNA_TYPE     <- "ncrna"

# =============================================================================
# Build miRNA ENST whitelist
# =============================================================================

if (file.exists(hs_mirna_annot_path)) {
  mirna_refs <- read_csv(hs_mirna_annot_path, show_col_types = FALSE) %>%
    filter(gene_biotype == "miRNA") %>%
    pull(ref) %>%
    unique()
  message("miRNA whitelist: ", length(mirna_refs), " ENST IDs from ", hs_mirna_annot_path)
} else {
  mirna_refs <- NULL
  warning("miRNA annotation file not found — ncRNA plots will include all transcripts.")
}

# =============================================================================
# Helper: load one idxstats_ncrna.txt file
# =============================================================================

load_idxstats <- function(path, rna_type = c("trna", "ncrna")) {
  rna_type <- match.arg(rna_type)
  if (!file.exists(path)) return(NULL)
  df <- read_tsv(path, col_names = c("ref", "length", "mapped", "unmapped"),
                 show_col_types = FALSE)
  df <- df[df$ref != "*", ]
  df <- df[!grepl("mito", df$ref, ignore.case = TRUE), ]
  if (rna_type == "trna")
    df <- df[!grepl("^tRNA-aaa[/_]", df$ref, ignore.case = TRUE), ]
  total <- sum(df$mapped, na.rm = TRUE)
  df$cpm <- if (total > 0) df$mapped / total * 1e6 else 0
  list(data = df[, c("ref", "cpm")], total_mapped = total)
}

# =============================================================================
# Load NT replicates
# =============================================================================

nt_samples <- read_csv(SAMPLE_SHEET, show_col_types = FALSE) %>%
  filter(Treatment == "NT")

rep_names <- nt_samples$Sample

rep_results <- lapply(seq_along(rep_names), function(i) {
  path <- file.path(hs_sample_data_dir, rep_names[i],
                    paste0("idxstats_", RNA_TYPE, ".txt"))
  res  <- load_idxstats(path, RNA_TYPE)
  if (is.null(res) || nrow(res$data) == 0) return(NULL)
  if (!is.null(mirna_refs))
    res$data <- res$data[res$data$ref %in% mirna_refs, ]
  if (nrow(res$data) == 0) return(NULL)
  list(data         = setNames(res$data, c("ref", paste0("rep", i))),
       total_mapped = res$total_mapped)
})

valid <- which(!sapply(rep_results, is.null))
if (length(valid) < 2)
  stop("Fewer than 2 NT replicates found in ", hs_sample_data_dir)

# =============================================================================
# Pairwise correlation plots (ncrna branch)
# =============================================================================

pairs   <- combn(valid, 2, simplify = FALSE)

# Write merged CPM table (all valid replicates) as underlying data
merged_cpm <- Reduce(function(a, b) inner_join(a, b, by = "ref"),
                     lapply(valid, function(i) rep_results[[i]]$data))
out_csv <- file.path(output_dir, "figS4A_mirna_replicate_scatter_data.csv")
write.csv(merged_cpm, out_csv, row.names = FALSE)
message("Saved: ", basename(out_csv))

out_pdf <- file.path(output_dir, "figS4A_mirna_replicate_scatter.pdf")
pdf(out_pdf, width = 7, height = 6.5)

for (pair in pairs) {
  r1 <- pair[1]; r2 <- pair[2]
  joined <- inner_join(rep_results[[r1]]$data,
                       rep_results[[r2]]$data, by = "ref")
  colnames(joined)[2:3] <- c("x", "y")

  thresh_x  <- 100 / rep_results[[r1]]$total_mapped * 1e6
  thresh_y  <- 100 / rep_results[[r2]]$total_mapped * 1e6
  lthresh_x <- log2(thresh_x + 1)
  lthresh_y <- log2(thresh_y + 1)

  joined$above_thresh <- joined$x >= thresh_x & joined$y >= thresh_y

  lx     <- log2(joined$x + 1)
  ly     <- log2(joined$y + 1)
  r2_val <- round(cor(lx, ly, use = "complete.obs")^2, 3)

  xlab <- bquote("Rep" ~ .(r1) ~ log[2](CPM + 1))
  ylab <- bquote("Rep" ~ .(r2) ~ log[2](CPM + 1))

  joined$color_grp <- ifelse(joined$above_thresh, "high", "low")
  n_above          <- sum(joined$above_thresh)
  stats_label      <- paste0("R² = ", r2_val, "\nmiRNA: ", n_above)
  y_mid            <- mean(range(ly, na.rm = TRUE))

  p <- ggplot(joined, aes(x = log2(x + 1), y = log2(y + 1),
                          color = color_grp)) +
    scale_color_manual(values = c("high" = "blue", "low" = "#DEDEDE"),
                       guide = "none") +
    geom_abline(slope = 1, intercept = 0,
                color = "grey75", linetype = "dashed", linewidth = 0.4) +
    geom_point(alpha = 0.65, size = 1.5) +
    geom_vline(xintercept = lthresh_x,
               linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = lthresh_y,
               linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_smooth(data = subset(joined, above_thresh),
                method = "lm", formula = y ~ x,
                se = FALSE, color = "grey30", linewidth = 0.6) +
    annotate("text", x = Inf, y = y_mid,
             hjust = -0.15, vjust = 0.5,
             label = stats_label, size = 5,
             color = "grey20", lineheight = 1.6) +
    coord_cartesian(clip = "off") +
    labs(x = xlab, y = ylab,
         title = paste0("WT NT  rep", r1, " vs rep", r2, "  (miRNA)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "bottom",
          plot.margin     = margin(t = 5, r = 90, b = 5, l = 5, unit = "pt"))

  print(p)
}

dev.off()
message("Saved: ", out_pdf)
