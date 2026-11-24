# Site-calling pipeline for BID-Seq tRNA pseudouridine analysis

# Edit config.R to control all behaviour.

source("config.R")

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)
library(purrr)

load_pileup <- function(pileup_dir, samples_path) {
  if (!file.exists(samples_path)) stop("samples CSV not found: ", samples_path)
  message("Reading samples from: ", samples_path)
  meta <- read_csv(samples_path, show_col_types = FALSE) %>%
    mutate(Treatment = if_else(Treatment == "NT", "NT", "BS"))
  map_dfr(1:nrow(meta), function(i) {
    s    <- meta$Sample[i]
    cond <- meta$Treatment[i]
    grp  <- meta$Group[i]
    p    <- file.path(pileup_dir, paste0(s, ".csv"))
    if (!file.exists(p)) { message("MISSING FILE: ", p); return(NULL) }
    read_csv(p, show_col_types = FALSE) %>%
      filter(!is.na(depth)) %>%
      mutate(cpm = 1e6 * depth / sum(depth, na.rm = TRUE), logcpm = log2(cpm),
             Tissue = grp, Condition = cond, Sample_ID = s) %>%
      select(ref, start, base, depth, del_rate, logcpm, Tissue, Condition, Sample_ID)
  })
}

# =============================================================================
# PART 1 — Ψ site calling
# =============================================================================

# Reference data (Sprinzl and FASTA needed only when use_sprinzl / use_bacs)
if (use_sprinzl) {
  sprinzl_meta <- read_csv(file.path(ref_dir, sprinzl_lookup_file),
                            show_col_types = FALSE) %>%
    mutate(start = local_pos - 1) %>% select(ref, start, sprinzl_pos)

  fasta_raw <- tibble(line = readLines(file.path(ref_dir, fasta_file))) %>%
    mutate(is_header = str_starts(line, ">"),
           trna_id   = if_else(is_header, str_remove(line, ">"), NA_character_)) %>%
    fill(trna_id, .direction = "down") %>% filter(!is_header) %>%
    group_by(trna_id) %>%
    summarize(seq = paste(line, collapse = ""), .groups = "drop") %>%
    mutate(full_len = nchar(seq))
}

if (use_bacs) {
  orthogonal_data <- read.csv(file.path(ref_dir, bacs_sites_file), header = TRUE) %>%
    mutate(match_name = str_replace_all(name, "_", "-"))
}

all_data <- load_pileup(pileup_dir, file.path(data_dir, samples_csv))
if (nrow(all_data) == 0) stop("No pileup files loaded. Check MISSING FILE messages above for the exact paths being searched.")

# Fragment lengths — used to compute pos from 3' end for 3p fragments
frag_lens <- all_data %>%
  group_by(ref) %>%
  summarize(frag_len = max(start) + 1L, .groups = "drop")

groups_to_call <- unique(all_data$Tissue)
mod_sites_all  <- list()

for (grp in groups_to_call) {

  sites <- all_data %>%
    filter(Tissue == grp) %>%
    group_by(ref, start, base) %>%
    summarize(
      n_nt        = sum(Condition == "NT"),
      n_bs        = sum(Condition == "BS"),
      all_nt_pass = n_nt > 0 && all(del_rate[Condition == "NT"] < NT_cutoff  & depth[Condition == "NT"] > min_depth),
      all_bs_pass = n_bs > 0 && all(del_rate[Condition == "BS"] > BS_min_del & depth[Condition == "BS"] > min_depth),
      .groups = "drop"
    ) %>%
    filter(base == "T", all_nt_pass, all_bs_pass)

  if (use_sprinzl) {
    sites <- sites %>%
      mutate(
        first_trna_full = map_chr(str_split(ref, "/"), 1),
        first_trna_base = str_remove(first_trna_full, "_[35]p$"),
        is_5p           = str_detect(first_trna_full, "5p$")
      ) %>%
      left_join(fasta_raw, by = c("first_trna_base" = "trna_id")) %>%
      mutate(local_pos_1 = start + 1,
             mature_trna_pos = if_else(is_5p, local_pos_1, full_len - 35 + local_pos_1)) %>%
      left_join(sprinzl_meta, by = c("ref", "start")) %>%
      drop_na(sprinzl_pos)
  } else {
    sites <- sites %>%
      left_join(frag_lens, by = "ref") %>%
      mutate(pos = if_else(str_detect(ref, "3p$"), frag_len - start, start + 1L)) %>%
      select(-frag_len)
  }

  if (use_bacs) {
    sites <- sites %>%
      mutate(Orthogonal_Match = map2_lgl(ref, mature_trna_pos, ~ {
        trnas <- str_remove_all(str_split(.x, "/")[[1]], "_[35]p$")
        any(orthogonal_data$match_name %in% trnas & orthogonal_data$pos == .y)
      }))
  }

  sites <- sites %>%
    mutate(
      simplified_ref = clean_labels(ref),
      Tissue         = grp,
      End_Type       = if_else(str_detect(ref, "3p$"), "3p",
                               if_else(str_detect(ref, "5p$"), "5p", "other"))
    )

  if (use_sprinzl) {
    sites <- sites %>%
      select(Tissue, ref, simplified_ref, End_Type, mature_trna_pos, sprinzl_pos,
             base, n_nt, n_bs, any_of("Orthogonal_Match"))
  } else {
    sites <- sites %>%
      select(Tissue, ref, simplified_ref, End_Type, pos, base, n_nt, n_bs)
  }

  sites <- sites %>% arrange(simplified_ref, if (use_sprinzl) mature_trna_pos else pos)

  out <- file.path(pileup_dir, paste0("modification_sites_", grp, ".csv"))
  write_csv(sites, out)
  message(nrow(sites), " Ψ sites called for '", grp, "' → ", out)

  mod_sites_all[[grp]] <- sites
}

all_mod_sites <- bind_rows(mod_sites_all)
write_csv(all_mod_sites,
          file.path(pileup_dir, "modification_sites_all_groups.csv"))
message("Combined site calls → modification_sites_all_groups.csv")

# Export raw per-sample BS deletion rates at called sites (used by fig2D, fig2E)
# pos convention: 1-based from 5' end for 5p fragments, from 3' end for 3p fragments
site_keys <- all_mod_sites %>%
  mutate(key = paste(ref, pos, sep = "_")) %>%
  pull(key) %>% unique()

called_site_data <- all_data %>%
  filter(Condition == "BS", base == "T") %>%
  left_join(frag_lens, by = "ref") %>%
  mutate(
    End_Type = if_else(str_detect(ref, "3p$"), "3p",
                       if_else(str_detect(ref, "5p$"), "5p", "other")),
    pos      = if_else(End_Type == "3p", frag_len - start, start + 1L),
    key      = paste(ref, pos, sep = "_")
  ) %>%
  filter(key %in% site_keys) %>%
  mutate(simplified_ref = clean_labels(ref)) %>%
  select(simplified_ref, ref, pos, End_Type, Tissue, Sample_ID, depth, del_rate)

write_csv(called_site_data, file.path(pileup_dir, "called_site_data.csv"))
message("Per-sample called-site rates → called_site_data.csv")

# =============================================================================
# PART 2 — tRF position summary for heatmap scripts (fig1G / fig2C)
# =============================================================================

trfs_data <- load_pileup(pileup_dir, file.path(data_dir, samples_csv))
if (nrow(trfs_data) == 0) stop("No pileup files loaded (Part 2). Check MISSING FILE messages above.")

summary_trfs <- trfs_data %>%
  group_by(ref, start, base, Tissue) %>%
  summarize(
    mean_logcpm        = mean(logcpm[Condition == "NT"], na.rm = TRUE),
    mean_del_nt        = mean(del_rate[Condition == "NT"], na.rm = TRUE),
    mean_del_bs        = mean(del_rate[Condition == "BS"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mean_del_rate_cor = if_else(base == "T", mean_del_bs, 0))

if (use_trna_lookup) {
  ref_lookup   <- read_csv(file.path(ref_dir, trna_lookup_file), show_col_types = FALSE)
  sprinzl_trfs <- read_csv(file.path(ref_dir, sprinzl_lookup_file), show_col_types = FALSE) %>%
    mutate(start = local_pos - 1) %>% select(ref, start, sprinzl_pos)
  means_df <- summary_trfs %>%
    left_join(ref_lookup,   by = c("ref" = "original_id")) %>%
    left_join(sprinzl_trfs, by = c("ref", "start")) %>%
    drop_na(simplified_ref)
} else {
  means_df <- summary_trfs %>%
    mutate(simplified_ref = clean_labels(ref))
}

means_df <- means_df %>%
  mutate(
    Type     = if_else(str_detect(ref, "mito"), "mito", "cyto"),
    End_Type = if_else(str_detect(ref, "3p$"), "3p",
                       if_else(str_detect(ref, "5p$"), "5p", "other")),
    AA       = str_extract(ref, "(?<=tRNA-)[A-Za-z]+")
  )

out_means <- file.path(pileup_dir, "means_df.csv")
write_csv(means_df, out_means)
message("tRF position summary → ", out_means)
