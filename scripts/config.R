# =============================================================================
# config.R — user-specified parameters for all analysis scripts
# See README.md for a full description of each parameter.
# =============================================================================

# --- Directories -------------------------------------------------------------
data_dir           <- "/Users/christianfagre/Yale University/Gilbert Lab/mirna_psi_manuscript/11162026_Compiled_Data"
trf_hs_folder      <- "trf_hs"       # human tRF pileup CSVs
trf_mm_folder      <- "trf_mm"       # mouse tRF pileup CSVs
ncrna_hs_folder    <- "ncrna_hs"     # human ncRNA pileup CSVs (fig3B/C)
ncrna_mm_folder    <- "ncrna_mm"     # mouse ncRNA pileup CSVs (fig3D)
ref_folder         <- "refs"         # reference files
sample_data_folder <- "sample_data"  # per-sample pipeline outputs (<sample>/<file>)
hs_sample_data_dir <- "/Users/christianfagre/Yale University/Gilbert Lab/Data/CF092"  # human (WT) sample folders
deseq_folder       <- "deseq"        # DESeq2 inputs: Conditions.txt + <sample>/gene_counts.txt
output_folder      <- "figures"      # all figure PNGs and CSVs

# Derived full paths — do not edit
trf_hs_dir      <- file.path(data_dir, trf_hs_folder)
trf_mm_dir      <- file.path(data_dir, trf_mm_folder)
ncrna_hs_dir    <- file.path(data_dir, ncrna_hs_folder)
ncrna_mm_dir    <- file.path(data_dir, ncrna_mm_folder)
ref_dir         <- file.path(data_dir, ref_folder)
sample_data_dir <- file.path(data_dir, sample_data_folder)
deseq_data_dir  <- file.path(data_dir, deseq_folder)
output_dir      <- file.path(data_dir, output_folder)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Sample sheets -----------------------------------------------------------
# Each CSV: columns Sample, Treatment (NT or RBS), Group
samples_hs_csv <- "samples_hs.csv"   # human samples (fig3A/B/C, fig4)
samples_mm_csv <- "samples_mm.csv"   # mouse samples (fig3D, fig1/2 tRF)

# Order of groups in all plots — set to NULL for alphabetical
sample_order <- c("Cerebellum", "Lung", "Liver")

# --- Reference files (filenames only; looked up in ref_dir) ------------------
sprinzl_lookup_file  <- "sprinzl_lookup.csv"
fasta_file           <- "hg38-mature-tRNAs_clean.fa"
hs_modomics_tsv      <- "/Users/christianfagre/Yale University/Gilbert Lab/Data/CF092/refs/hs_modomics_tRNA_alignment.tsv"  # Sprinzl coordinate map
hs_mirna_annot_path  <- "/Users/christianfagre/Yale University/Gilbert Lab/Data/CF092/ncrna/WT_1_NT_annotated.csv"  # miRNA ENST whitelist (gene_biotype == "miRNA")
bacs_sites_file     <- "BACS_sites.csv"
trna_lookup_file    <- "trna_lookup.csv"
trf_hs_fasta_file   <- "hg38_tRFs.fa"
trf_mm_fasta_file   <- "mm39_tRFs.fa"

# --- Site-calling thresholds -------------------------------------------------
NT_cutoff  <- 0.01   # max deletion rate in untreated to pass
BS_min_del <- 0.10   # min deletion rate in bisulfite to call a site (all replicates must exceed)
min_depth  <- 50     # minimum read depth at a position (NT samples)

# --- Coordinate output for site calls ----------------------------------------
# TRUE:  map to Sprinzl + mature tRNA positions (requires fasta_file; fig1 mode)
# FALSE: use local 1-based position only (fig2 mode)
use_sprinzl <- FALSE

# --- BACS orthogonal check ---------------------------------------------------
# TRUE: add Orthogonal_Match column (fig1 mode); FALSE: skip (fig2 mode)
use_bacs <- FALSE

# --- tRF summary: simplified_ref computation --------------------------------
# TRUE:  use trna_lookup.csv (pre-computed mapping; fig1 mode)
# FALSE: use clean_labels() regex function (fig2 mode)
use_trna_lookup <- FALSE

source("helpers.R")

# =============================================================================
# Figure-specific parameters
# =============================================================================

# --- Figure 1B: correlation scatter plot -------------------------------------
nt_counts_file <- "../CF092/CF092-1_counts.txt"
bs_counts_file <- "../CF092/CF092-13_counts.txt"
counts_prefix  <- "^hsa-mir"   # regex to filter feature names

# --- Figure 1G / 2C: tRF deletion-rate heatmap -------------------------------
heatmap_tissue <- "Cerebellum"  # tissue/condition to plot
heatmap_end    <- "3p"       # "3p" or "5p"
heatmap_type   <- "cyto"     # "cyto" or "mito"
heatmap_top_n  <- 50
heatmap_aa     <- NULL       # e.g. c("Val","Gly") or NULL for all

# --- tRF pileup heatmaps (figS2C/D/E, figS3C/D/E/F) -------------------------
trf_min_depth <- 10   # min read depth per position for tRF coverage/deletion heatmaps

# --- Figure 2C: tRF expression heatmap (idxstats) ----------------------------
fig2c_cpm_filter <- 200   # min total CPM across all samples to retain a tRF

# --- Figure 2D: deletion heatmap at specific modification sites --------------
fig2d_end       <- "5p" # 5p or 3p modified fragment
fig2d_positions <- c(27,28)   # vector of 1-based local positions. For 3p fragment
                           # sites, position is relative to 3' end
                           # (e.g. for position 55, enter 22)
fig2d_type      <- "cyto"

# --- Figure 2E: barplots + ANOVA at a specific site --------------------------
fig2e_ref <- "tRNA-Ala-AGC-2"   # regex to match tRNA ref names
fig2e_pos <- 17                  # 1-based local position
fig2e_end <- "5p"                # "5p" or "3p"

# --- Figure 3A: miRNA count correlation scatter ------------------------------
fig3a_nt_counts_file <- "../CF092/CF092-1_counts.txt"   # untreated count table
fig3a_bs_counts_file <- "../CF092/CF092-13_counts.txt"  # bisulfite count table
fig3a_cpm_threshold  <- 100          # CPM threshold for quadrant highlight

# --- Figure 3B: miRNA deletion rate scatter (HepG2) -------------------------
# Uses ncrna_hs_dir and samples_hs_csv
fig3b_tissue       <- NULL   # NULL = use all samples in samples_hs.csv
fig3b_min_depth    <- 20     # CF092/scripts/R_plots_pileup_ncrna.R: MIN_DEPTH
fig3b_repro_del_min <- 0.2   # CF092/scripts/R_plots_pileup_ncrna.R: REPRO_DEL_MIN

# --- Figure 3D: miRNA expression heatmap across mouse tissues (idxstats) -----
fig3d_idxstats_file <- "idxstats.txt"   # per-sample idxstats file (mmu- miRNA names)
fig3d_cpm_filter    <- 200              # min total CPM across samples (rowSums > 200)

# --- Figure 3E: miRNA deletion rate scatter, mouse tissues ------------------
# Groups (Group column values in samples_mm.csv) to include; order = facet order
fig3d_groups        <- c("Liver", "Lung", "Cerebellum")
fig3e_repro_del_min <- 0.2    # CF112/scripts/R_plots_pileup_ncrna.R: RBS_DEL_MIN

# --- Figure 4B/C: DESeq2 pipeline → CDF of Log2FC binned by TargetScan site type
fig4_conditions_file   <- "Conditions.txt"         # sample sheet in deseq_data_dir (cols: ID, Sample)
fig4_counts_file       <- "gene_counts.txt"        # HTSeq count filename inside each sample dir
fig4_ref_condition     <- "Scramble-U"             # reference level for DESeq2
fig4_min_counts        <- 10                       # min counts per gene per sample
fig4_min_samples       <- 3                        # min samples satisfying fig4_min_counts
fig4_gene_site_map_file <- "refseq_gene_site_map.csv"  # pre-built map in ref_dir
                                                        # cols: transcript_id, gene_name,
                                                        #       miRNA_site_category, miRNA_site_simple
# Named list of comparisons: label → DESeq2 contrast vector c(factor, numerator, denominator)
fig4_comparisons <- list(
  "miR-3068-U vs Scramble-U"  = c("sampleCondition", "miR3068-U",   "Scramble-U"),
  "miR-3068-U vs miR-3068-Y"  = c("sampleCondition", "miR3068-17Y", "miR3068-U")
)

# --- Figure S5: CDF of Log2FC for miR3068-5p-Ψ vs Scramble-Ψ -----------------
figS5_comparison <- c("sampleCondition", "miR3068-17Y", "Scramble-Y")
