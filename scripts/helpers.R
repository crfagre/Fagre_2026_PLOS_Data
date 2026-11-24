library(purrr)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)

clean_labels <- function(labels) {
  map_chr(labels, function(label) {
    parts  <- str_split(label, "/")[[1]]
    parsed <- str_match(parts, "(?:mt-)?tRNA-([A-Za-z]+)-([A-Za-z0-9]+)-([0-9]+).*_[35]p")
    if (any(is.na(parsed[, 1]))) return(label)
    prefix <- if (str_detect(parts[1], "^mt-")) "mt-" else ""
    aa <- unique(parsed[, 2])
    ac <- unique(parsed[, 3])
    aa_ac_str <- if (length(aa) == 1) {
      paste0(prefix, aa, "-", paste(ac, collapse = "/"))
    } else {
      paste(unique(paste0(prefix, aa, "-", parsed[, 3])), collapse = "/")
    }
    paste0(aa_ac_str, "-", paste(sort(as.numeric(unique(parsed[, 4]))), collapse = "/"))
  })
}

# =============================================================================
# tRF helper functions
# =============================================================================

# Extract isoacceptor family from slash-collapsed tRF ref string.
# e.g. "tRNA-Cys-GCA-10_5p/tRNA-Cys-GCA-11_5p" -> "tRNA-Cys-GCA"
parse_trfs_family <- function(ref) {
  first  <- str_extract(ref, "^[^/]+")
  strip  <- str_remove(first, "_[35]p$")
  family <- str_remove(strip, "-[0-9]+$")
  family
}

# Extract 5p / 3p fragment type from trfs ref string.
parse_trfs_type <- function(ref) {
  first    <- str_extract(ref, "^[^/]+")
  trf_type <- str_extract(first, "[35]p$")
  if_else(is.na(trf_type), "unknown", trf_type)
}

# Load a tRF pileup CSV and add trf_family + trf_type columns.
load_trfs <- function(path) {
  if (!file.exists(path)) {
    warning("tRF file not found: ", path)
    return(NULL)
  }
  df <- read_csv(path, show_col_types = FALSE)
  df <- df[!grepl("_mito|^tRNA-aaa[_/]", df$ref, ignore.case = TRUE), ]
  df$trf_family <- parse_trfs_family(df$ref)
  df$trf_type   <- parse_trfs_type(df$ref)
  df
}

# Parse hs_modomics_tRNA_alignment.tsv and build a data frame mapping
# (trf_family, trf_type, start) -> Sprinzl position label.
# Returns: list(lookup = tibble | NULL, sprinzl_order = character vector)
build_sprinzl_lookup <- function(path, trf_families, frag_len = 35L) {
  if (!file.exists(path)) {
    warning("Modomics file not found: ", path, " — Sprinzl mapping disabled.")
    return(list(lookup = NULL, sprinzl_order = character(0)))
  }
  raw <- suppressMessages(
    read_tsv(path, col_names = FALSE, show_col_types = FALSE, progress = FALSE)
  )
  n_cols <- ncol(raw)
  sprinzl_labels <- unlist(raw[3, 9:n_cols], use.names = FALSE)
  sprinzl_labels[is.na(sprinzl_labels) | sprinzl_labels == ""] <- NA_character_
  sprinzl_order  <- sprinzl_labels[!is.na(sprinzl_labels)]

  entries <- lapply(5:nrow(raw), function(i) {
    aa     <- as.character(raw[[3]][i])
    ac_raw <- as.character(raw[[4]][i])
    chars  <- unlist(raw[i, 9:n_cols], use.names = FALSE)
    chars[is.na(chars)] <- "-"
    non_gap        <- which(chars != "-")
    abs_to_sprinzl <- sprinzl_labels[non_gap]
    total_len      <- length(abs_to_sprinzl)
    ac <- toupper(ac_raw)
    ac <- chartr("UI", "TA", ac)
    ac <- gsub("[^ACGT]", "N", ac)
    list(aa = aa, ac = ac, abs_to_sprinzl = abs_to_sprinzl, total_len = total_len)
  })

  unique_fams <- unique(trf_families[!is.na(trf_families)])
  lookup <- bind_rows(lapply(unique_fams, function(fam) {
    m <- str_match(fam, "^tRNA-([A-Za-z]+)-([A-Za-z]+)$")
    if (is.na(m[1])) return(NULL)
    aa_q <- m[2]; ac_q <- m[3]

    entry <- NULL
    for (e in entries) {
      if (e$aa == aa_q && e$ac == ac_q) { entry <- e; break }
    }
    if (is.null(entry)) {
      for (e in entries) {
        if (e$aa != aa_q || nchar(e$ac) != nchar(ac_q)) next
        e_ch <- strsplit(e$ac, "")[[1]]; q_ch <- strsplit(ac_q, "")[[1]]
        if (all(e_ch == q_ch | e_ch == "N")) { entry <- e; break }
      }
    }
    if (is.null(entry)) {
      for (e in entries) { if (e$aa == aa_q) { entry <- e; break } }
    }
    if (is.null(entry)) return(NULL)

    s5  <- seq_len(min(frag_len, entry$total_len)) - 1L
    sp5 <- entry$abs_to_sprinzl[s5 + 1L]

    s3      <- seq_len(frag_len) - 1L
    sp3_base <- 76L - frag_len + 1L
    sp3     <- as.character(sp3_base + s3)

    bind_rows(
      tibble(trf_family = fam, trf_type = "5p", start = as.integer(s5),  sprinzl = sp5),
      tibble(trf_family = fam, trf_type = "3p", start = as.integer(s3), sprinzl = sp3)
    )
  }))

  list(lookup = lookup, sprinzl_order = sprinzl_order)
}

# Build position matrix using Sprinzl labels as columns.
make_trf_matrix <- function(data, value_col, col_order = NULL) {
  data <- data %>% filter(!is.na(sprinzl))
  if (nrow(data) == 0) return(NULL)
  mat <- data %>%
    group_by(trf_family, sprinzl) %>%
    summarise(val = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = sprinzl, values_from = val) %>%
    column_to_rownames("trf_family") %>%
    as.matrix()
  if (!is.null(col_order)) {
    keep <- intersect(col_order, colnames(mat))
    if (length(keep) > 0) mat <- mat[, keep, drop = FALSE]
  }
  mat
}

# Build position matrix using slash-collapsed ref as row key (absolute positions).
make_trf_matrix_ref <- function(data, value_col) {
  if (nrow(data) == 0) return(NULL)
  mat <- data %>%
    group_by(ref, start) %>%
    summarise(val = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = start, values_from = val) %>%
    column_to_rownames("ref") %>%
    as.matrix()
  mat[, order(as.integer(colnames(mat))), drop = FALSE]
}

# Cluster rows on NA-imputed copy, return matrix reordered with NAs intact.
cluster_rows_na <- function(mat) {
  if (nrow(mat) < 2) return(mat)
  imp <- mat; imp[is.na(imp)] <- 0
  ord <- tryCatch(hclust(dist(imp))$order, error = function(e) seq_len(nrow(mat)))
  mat[ord, , drop = FALSE]
}
