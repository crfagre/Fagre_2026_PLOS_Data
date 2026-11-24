# Figure 1D: Venn diagram of tRNA Ψ site overlap (BID-Seq vs. BACS) + overlap metrics
#
# Requires pipeline_site_calling.R to have been run first.
# Edit config.R to set target_group, data_dir, bacs_sites_file.

source("config.R")

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggVennDiagram)
library(ggplot2)

#------------------------------------------------------------------------------
# Load pipeline output and BACS dataset
#------------------------------------------------------------------------------

sites_file <- file.path(pileup_dir, paste0("modification_sites_", target_group, ".csv"))
mod_sites  <- read_csv(sites_file, show_col_types = FALSE)

orthogonal_data <- read.csv(file.path(ref_dir, bacs_sites_file), header = TRUE) %>%
  mutate(match_name = str_replace_all(name, "_", "-"))

#------------------------------------------------------------------------------
# Expand to individual tRNA + position keys
#------------------------------------------------------------------------------

my_expanded_sites <- mod_sites %>%
  mutate(individual_trna = str_split(ref, "/")) %>%
  unnest(individual_trna) %>%
  mutate(
    individual_trna = str_remove(individual_trna, "_[35]p$"),
    site_id         = paste(individual_trna, mature_trna_pos, sep = "_")
  ) %>%
  pull(site_id) %>%
  unique()

bacs_expanded_sites <- orthogonal_data %>%
  mutate(site_id = paste(match_name, pos, sep = "_")) %>%
  pull(site_id) %>%
  unique()

#------------------------------------------------------------------------------
# Venn diagram
#------------------------------------------------------------------------------

venn_list <- list(
  "BID-Seq"      = my_expanded_sites,
  "BACS Dataset" = bacs_expanded_sites
)

p_venn <- ggVennDiagram(venn_list, label_alpha = 0, set_color = "black") +
  scale_fill_gradient(low = "white", high = "dodgerblue") +
  theme(legend.position = "none") +
  labs(title = paste("tRNA Ψ Site Overlap —", target_group)) +
  scale_x_continuous(expand = expansion(mult = 0.2))

out_file <- file.path(output_dir, paste0("fig1D_venn_", target_group, ".png"))
ggsave(out_file, plot = p_venn, width = 6, height = 5, bg = "white", dpi = 300)
message("Saved: ", out_file)
print(p_venn)

#------------------------------------------------------------------------------
# BACS overlap metrics
#------------------------------------------------------------------------------

n_my   <- length(my_expanded_sites)
n_bacs <- length(bacs_expanded_sites)
n_both <- length(intersect(my_expanded_sites, bacs_expanded_sites))

cat("\n=== BACS OVERLAP METRICS (", target_group, ") ===\n", sep = "")
cat("BID-Seq sites:                 ", n_my,               "\n")
cat("BACS sites:                    ", n_bacs,              "\n")
cat("Overlapping:                   ", n_both,              "\n")
cat("BID-Seq only:                  ", n_my   - n_both,     "\n")
cat("BACS only:                     ", n_bacs - n_both,     "\n")
cat("Fraction of BID-Seq in BACS:   ",
    round(n_both / n_my * 100, 1), "%\n")

#------------------------------------------------------------------------------
# CSV exports: expanded site lists + overlap summary
#------------------------------------------------------------------------------

overlap  <- intersect(my_expanded_sites, bacs_expanded_sites)
only_bid <- setdiff(my_expanded_sites,   bacs_expanded_sites)
only_bac <- setdiff(bacs_expanded_sites, my_expanded_sites)

sites_csv <- bind_rows(
  data.frame(site_id = my_expanded_sites,    set = "BID-Seq"),
  data.frame(site_id = bacs_expanded_sites,  set = "BACS"),
  data.frame(site_id = overlap,              set = "Overlap")
)

overlap_summary <- data.frame(
  group             = target_group,
  n_bidseq          = n_my,
  n_bacs            = n_bacs,
  n_overlap         = n_both,
  n_bidseq_only     = n_my   - n_both,
  n_bacs_only       = n_bacs - n_both,
  pct_bidseq_in_bacs = round(n_both / n_my * 100, 1)
)

out_prefix <- file.path(output_dir, paste0("fig1D_", target_group))
write_csv(sites_csv,       paste0(out_prefix, "_site_lists.csv"))
write_csv(overlap_summary, paste0(out_prefix, "_overlap_summary.csv"))
message("Saved CSVs: ", out_prefix, "_site_lists.csv + _overlap_summary.csv")
