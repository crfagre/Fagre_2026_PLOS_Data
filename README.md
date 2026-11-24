# Transcriptome-wide profiling reveals rare pseudouridylation of mammalian microRNAs
# Christian Fagre, Yixuan Chen, and Wendy Gilbert
# Supplemental Data and Code


All figure data is located in data/, in respective figure directories
('e.g. Fig 1B data is in 'data/Fig 1/B'). Each figure plot can be reproduced by
running the scripts in this repository.

## Code Overview (Located in scripts/)

After alignment, pileup files are generated using 'count_dels_pileup.py' on
untreated and bisulfite-treated samples:


python count_dels_pileup_1.2.py \
     -n [bam file] \
     -f [fasta reference] \
     -o [output prefix (optional)] \
     -t [threads (optional)] \
     -w [window size (optional)] \
     -d [minimum depth (optional)]
     
Pseudouridine sites are called with 'call_sites.R'



**Figure 1:**
```
1. call_sites.R   # Calls psi sites 
2. fig1B_correlation_plot.R  # Figure 1B — NT vs BS correlation scatter
3. fig1D_BACS_overlap.R      # Figure 1D — BID-Seq / BACS Venn + overlap metrics
4. fig1G_trf_heatmap.R       # Figure 1G — tRF deletion-rate heatmap
```

**Figure 2 (all groups, local coordinates, no BACS):**
```
1. pipeline_site_calling.R    # call Ψ sites for all groups
2. fig2B_trf_heatmap.R       # Figure 2B — tRF coverage + deletion heatmap
3. fig2C_expression_heatmap.R # Figure 2C — tRF expression heatmap (idxstats)
4. fig2D_site_heatmap.R      # Figure 2D — differential Ψ site heatmap
5. fig2E_barplots.R          # Figure 2E — barplots with ANOVA/Tukey
```

> **Note:** Fig2B, 2D and 2E require the pipeline to have been run first. Fig2C reads idxstats files directly and has no pipeline dependency.

**Figure 3 (miRNA Ψ identification):**
```
1. fig3A_mirna_correlation.R           # Figure 3A — miRNA NT vs BS count scatter
2. fig3B_mirna_deletion_scatter.R      # Figure 3B — deletion rate vs depth (HepG2)
3. fig3C_mirna_expression_heatmap.R    # Figure 3C — miRNA expression heatmap (idxstats)
4. fig3D_mirna_deletion_scatter_tissues.R  # Figure 3D — deletion rate vs depth, mouse tissues
```

> **Note:** Fig3B and 3D require miRNA pileup CSVs in `mirna_pileup_dir`. Fig3A and 3C have no pipeline dependency.

**Figure 4 (miR-3068-5p modification impacts repression):**
```
1. fig4BC_deseq_cdf.R  # Figure 4B/C — DESeq2 → CDF Log2FC, one panel per comparison
```

> **Note:** Complete pipeline — no pre-processed CSVs needed. Requires HTSeq count files, a `Conditions.txt` sample sheet, and a TargetScan TSV in `ref_dir` Add or remove comparisons in `fig4_comparisons` in `config.R` to change which panels are produced.

---

## Parameters (`config.R`)

### Directories

| Parameter | Description |
|-----------|-------------|
| `data_dir` | Root data directory (absolute path). |
| `pileup_folder` | Subfolder under `data_dir` containing tRF pileup CSVs. |
| `ref_folder` | Subfolder under `data_dir` containing all reference files (TargetScan, Sprinzl lookup, BACS, etc.). |
| `sample_data_folder` | Subfolder under `data_dir` containing per-sample pipeline outputs, one subdirectory per sample named by `Sample` from `samples.csv` (e.g. `sample_data/Lvr-1-NT/idxstats_trna.txt`). |
| `mirna_pileup_folder` | Subfolder under `data_dir` containing miRNA pileup CSVs (fig3B/D). |
| `deseq_folder` | Subfolder under `data_dir` containing DESeq2 inputs: `Conditions.txt` and one subdirectory per sample holding `gene_counts.txt` (fig4). |
| `output_folder` | Subfolder under `data_dir` where all figure PNGs and CSVs are written. Created automatically if it does not exist. |

`pileup_dir`, `ref_dir`, `sample_data_dir`, `mirna_pileup_dir`, `deseq_data_dir`, and `output_dir` are derived automatically in `config.R`. Scripts use those; do not edit them directly.

### Sample sheet

| Parameter | Description |
|-----------|-------------|
| `samples_csv` | Filename of a CSV with a header row and three columns: `Sample` (file base name), `Treatment` (`NT` or anything else for bisulfite), `Group` (tissue / condition label). Must exist inside `data_dir` (not inside `pileup_dir`, since the same sheet applies to all RNA classes). |

The `Group` column groups replicates and is what `target_group` / `heatmap_tissue` match against.

### Reference files

| Parameter | Description |
|-----------|-------------|
| `sprinzl_lookup_file` | CSV mapping `ref`, `local_pos` → `sprinzl_pos` |
| `fasta_file` | FASTA of mature tRNA sequences, used to compute mature-tRNA positions |
| `bacs_sites_file` | BACS orthogonal Ψ-site table (CSV with columns `name`, `pos`) |
| `trna_lookup_file` | CSV mapping `original_id` → `simplified_ref` for tRF fragments |

All reference files are looked up in `ref_dir`.

### Site-calling thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NT_cutoff` | `0.01` | Maximum deletion rate in untreated samples to pass filtering |
| `BS_min_del` | `0.10` | Minimum deletion rate in bisulfite-treated samples to call a site |
| `min_depth` | `50` | Minimum read depth required at a position in both conditions |

The pipeline always calls sites for every group present in `samples.csv`.

### Pipeline mode flags

| Parameter | Description |
|-----------|-------------|
| `use_sprinzl` | `TRUE`: map local positions to Sprinzl + mature-tRNA coordinates; `FALSE`: keep local 1-based positions |
| `use_bacs` | `TRUE`: annotate sites with BACS orthogonal match; `FALSE`: skip BACS check |
| `use_trna_lookup` | `TRUE`: use `trna_lookup_file` to simplify tRF reference names; `FALSE`: use `clean_labels()` regex |

### Figure 1B — correlation scatter plot

| Parameter | Description |
|-----------|-------------|
| `nt_counts_file` | Path to the untreated read-count table (tab-separated, columns: name, length, count, …) |
| `bs_counts_file` | Path to the bisulfite-treated read-count table (same format) |
| `counts_prefix` | Regex to filter feature names, e.g. `"^hsa-mir"` for human miRNA or `"^mmu-tRNA"` for mouse tRNA |

### Figure 1G / 2C — tRF deletion-rate heatmap

| Parameter | Default | Description |
|-----------|---------|-------------|
| `heatmap_tissue` | `"shEmpty"` | Value in the `Tissue` column of `means_df.csv` to plot |
| `heatmap_end` | `"3p"` | Fragment end to show: `"3p"` or `"5p"` |
| `heatmap_type` | `"cyto"` | Fragment type: `"cyto"` (cytoplasmic) or `"mito"` (mitochondrial) |
| `heatmap_top_n` | `100` | Number of top-expressed tRFs to include (ranked by read depth at position 15) |
| `heatmap_aa` | `NULL` | Optional character vector of amino acids to restrict to, e.g. `c("Val", "Gly")`. `NULL` includes all. |

### Figure 2C — tRF expression heatmap

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fig2c_cpm_filter` | `200` | Minimum sum CPM across all samples for a tRF to be included in the heatmap |

Idxstats files are read from `<sample_data_dir>/<sample_id>/idxstats_trna.txt`, where `sample_id` comes from the `Sample` column of `samples.csv` (NT samples only).

### Figure 2D — differential Ψ site heatmap

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fig2d_positions` | `c(14)` | Integer vector of local positions to filter on (can span multiple positions) |
| `fig2d_end` | `"3p"` | Fragment end to show: `"3p"` or `"5p"` |
| `fig2d_type` | `"cyto"` | Fragment type: `"cyto"` or `"mito"` |

### Figure 2E — barplots

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fig2e_ref` | `"tRNA-Ala-AGC-2"` | tRF reference name (partial match) for the barplot |
| `fig2e_pos` | `17` | Local position to plot |
| `fig2e_end` | `"5p"` | Fragment end: `"3p"` or `"5p"` |

### Figure 3A — miRNA count correlation scatter

| Parameter | Description |
|-----------|-------------|
| `fig3a_nt_counts_file` | Path to the untreated read-count table (tab-separated, columns: name, length, count, …) |
| `fig3a_bs_counts_file` | Path to the bisulfite-treated read-count table (same format) |
| `fig3a_cpm_threshold` | CPM value used as the quadrant highlight threshold (dashed lines + colour split) |

### Figure 3B — miRNA deletion rate scatter (HepG2)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mirna_pileup_folder` | `"ncrna"` | Subfolder under `data_dir` containing miRNA pileup CSVs (same format as tRF pileups) |
| `fig3b_tissue` | `NULL` | Tissue to restrict to; `NULL` uses all samples in `samples.csv` |

U positions are filtered for `mean NT del_rate > NT_cutoff` and `mean BS depth ≥ min_depth`. Sites with `mean BS del_rate ≥ BS_min_del` are highlighted in red.

### Figure 3C — miRNA expression heatmap

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fig3c_idxstats_file` | `"idxstats_ncrna.txt"` | Per-sample idxstats filename (looked up in `<sample_data_dir>/<sample_id>/`) |
| `fig3c_cpm_filter` | `50` | Minimum total CPM across all samples for a miRNA to be included |

### Figure 3D — miRNA deletion rate scatter, mouse tissues

Uses `mirna_pileup_folder` and `samples.csv`. Produces a faceted plot with one panel per tissue; panels are ordered by `sample_order`.

### Figure 4B/C — DESeq2 pipeline → CDF of Log2FC by TargetScan site type

`fig4BC_deseq_cdf.R` is a complete pipeline. It runs DESeq2 from raw HTSeq count files, converts RefSeq IDs to gene names via biomaRt, categorises genes by predicted TargetScan site strength, and produces CDF plots + all output CSVs. No pre-processed input files are required.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fig4_conditions_file` | `"Conditions.txt"` | Sample sheet in `deseq_data_dir`; must have columns `ID` (matches a subdirectory name there) and `Sample` (condition label) |
| `fig4_counts_file` | `"gene_counts.txt"` | HTSeq count filename inside each sample's subdirectory of `deseq_data_dir` |
| `fig4_ref_condition` | `"Scramble-Y"` | Reference level for DESeq2 (passed to `relevel()`) |
| `fig4_min_counts` | `10` | Minimum counts per gene per sample for the row-filter |
| `fig4_min_samples` | `3` | Minimum number of samples satisfying `fig4_min_counts` |
| `fig4_gene_site_map_file` | `"refseq_gene_site_map.csv"` | Pre-built annotation map in `ref_dir`; must have columns `transcript_id`, `gene_name`, `miRNA_site_category`, `miRNA_site_simple`. Build once from biomaRt + TargetScan and save — avoids repeated API calls. |
| `fig4_comparisons` | see config | Named list: label → DESeq2 contrast vector `c(factor, numerator, denominator)`. Each entry produces one set of output files. |

Site annotation is loaded from the pre-built map, not queried at run time. KS tests vs "None" and summary statistics are printed and exported per comparison. DESeq2 runs once and is reused for all comparisons.

---

## Outputs

### Pipeline

| File | Description |
|------|-------------|
| `<pileup_dir>/modification_sites_<group>.csv` | Called Ψ sites for `target_group` (fig1 mode) |
| `<pileup_dir>/modification_sites_all_groups.csv` | Called Ψ sites for all groups (fig2 mode) |
| `<pileup_dir>/called_site_data.csv` | Raw per-sample BS deletion rates at all called sites. `pos` is 1-based from the 5' end for 5p fragments and from the 3' end for 3p fragments (position 1 = CCA). |
| `<pileup_dir>/means_df.csv` | Per-position mean coverage and deletion rates across all samples |

### Figure 1

| Script | Output files |
|--------|-------------|
| `fig1B_correlation_plot.R` | `<output_dir>/fig1B_correlation.png` — scatter plot |
| `fig1B_correlation_plot.R` | `<output_dir>/fig1B_correlation.csv` — log2 CPM values for both samples |
| `fig1D_venn_diagram.R` | `<output_dir>/fig1D_venn_<group>.png` — Venn diagram |
| `fig1D_venn_diagram.R` | `<output_dir>/fig1D_<group>_site_lists.csv` — BID-Seq, BACS, and overlap site IDs |
| `fig1D_venn_diagram.R` | `<output_dir>/fig1D_<group>_overlap_summary.csv` — counts per set |
| `fig1G_trf_heatmap.R` | `<output_dir>/fig1G_tRF-<end>_<type>_top<n>_<tissue>_deletion_rates.png` — heatmap |
| `fig1G_trf_heatmap.R` | `<output_dir>/tRF-<end>_<type>_top<n>_<tissue>_coverage.csv` — coverage matrix |
| `fig1G_trf_heatmap.R` | `<output_dir>/fig1G_tRF-<end>_<type>_top<n>_<tissue>_deletion_rates.csv` — deletion rate matrix |

### Figure 2

| Script | Output files |
|--------|-------------|
| `fig2B_trf_heatmap.R` | `<output_dir>/tRF-<end>_<type>_<tissue>_coverage.png` — coverage heatmap |
| `fig2B_trf_heatmap.R` | `<output_dir>/tRF-<end>_<type>_<tissue>_coverage.csv` — coverage matrix |
| `fig2B_trf_heatmap.R` | `<output_dir>/fig2C_tRF-<end>_<type>_<tissue>_deletion_rates.png` — deletion rate heatmap |
| `fig2B_trf_heatmap.R` | `<output_dir>/fig2C_tRF-<end>_<type>_<tissue>_deletion_rates.csv` — deletion rate matrix |
| `fig2C_expression_heatmap.R` | `<output_dir>/fig2C_Expression_<end>.png` — expression heatmap |
| `fig2C_expression_heatmap.R` | `<output_dir>/fig2C_Expression_<end>.csv` — log2-CPM matrix (tRFs × samples) |
| `fig2D_site_heatmap.R` | `<output_dir>/Master_Differential_Modification_Sites.csv` — per-site means + pairwise BH-adjusted t-tests |
| `fig2D_site_heatmap.R` | `<output_dir>/fig2D_Heatmap_Pos<pos>_<end>_<type>.png` — heatmap |
| `fig2D_site_heatmap.R` | `<output_dir>/fig2D_Heatmap_Pos<pos>_<end>_<type>.csv` — deletion rate matrix |
| `fig2E_barplots.R` | `<output_dir>/fig2E_barplot_<ref>_pos<n>_<end>.png` — barplot |
| `fig2E_barplots.R` | `<output_dir>/fig2E_barplot_<ref>_pos<n>_<end>.csv` — per-sample deletion rates |
| `fig2E_barplots.R` | `<output_dir>/fig2E_barplot_<ref>_pos<n>_<end>_tukey.csv` — Tukey HSD results |

### Figure 3

| Script | Output files |
|--------|-------------|
| `fig3A_mirna_correlation.R` | `<output_dir>/fig3A_mirna_correlation.png` — scatter plot |
| `fig3A_mirna_correlation.R` | `<output_dir>/fig3A_mirna_correlation.csv` — log2 CPM values |
| `fig3B_mirna_deletion_scatter.R` | `<output_dir>/fig3B_mirna_deletion_scatter_<tissue>.png` — scatter plot |
| `fig3B_mirna_deletion_scatter.R` | `<output_dir>/fig3B_mirna_deletion_scatter_<tissue>.csv` — per-position summary |
| `fig3C_mirna_expression_heatmap.R` | `<output_dir>/fig3C_mirna_expression_heatmap.png` — clustered heatmap |
| `fig3C_mirna_expression_heatmap.R` | `<output_dir>/fig3C_mirna_expression_heatmap.csv` — log2-CPM matrix |
| `fig3D_mirna_deletion_scatter_tissues.R` | `<output_dir>/fig3D_mirna_deletion_scatter_tissues.png` — faceted scatter |
| `fig3D_mirna_deletion_scatter_tissues.R` | `<output_dir>/fig3D_mirna_deletion_scatter_tissues.csv` — per-position summary |

### Figure 4

| Script | Output files |
|--------|-------------|
| `fig4BC_deseq_cdf.R` | `<output_dir>/fig4_normalized_counts.csv` — DESeq2 size-factor-normalized counts (all samples) |
| `fig4BC_deseq_cdf.R` | `<output_dir>/fig4_deseq_<label>.csv` — DESeq2 results with gene + site annotation (one per comparison) |
| `fig4BC_deseq_cdf.R` | `<output_dir>/fig4_mirna_cdf_<label>.png` — CDF plot with zoomed inset (one per comparison) |
| `fig4BC_deseq_cdf.R` | `<output_dir>/fig4_mirna_cdf_<label>_data.csv` — per-gene Log2FC + site bin used for plot |
| `fig4BC_deseq_cdf.R` | `<output_dir>/fig4_mirna_cdf_<label>_ks.csv` — KS test D + p-value per site bin vs "None" |

---

## Input file formats

**`samples.csv`** (comma-separated, header required):
```
Sample	Treatment	Group
Lvr-1-NT	NT	Liver
Lvr-1-BS	RBS	Liver
Lvr-2-NT	NT	Liver
Lvr-2-BS	RBS	Liver
```
Each row points to `<Sample>.csv` in `pileup_dir`. `Treatment` is normalised internally: `NT` stays `NT`, anything else (e.g. `RBS`, `BS`) is treated as bisulfite. `Group` is the tissue/condition label used for grouping.

**Per-sample pileup CSVs** (named `<sample_id>-NT.csv` and `<sample_id>-BS.csv`):
```
ref,start,base,depth,del_rate
tRNA-Val-TAC-1_3p,0,A,1234,0.001
...
```

**`BACS_sites.csv`** (fig1 only):
```
name,pos
tRNA_Val_TAC_1,27
...
```

**`idxstats_trna.txt`** (fig2D; SAMtools idxstats format):
```
tRNA-Val-TAC-1_3p  73  1234  0
...
```
