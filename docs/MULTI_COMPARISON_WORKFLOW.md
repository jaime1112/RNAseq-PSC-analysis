# Multi-Comparison Workflow Guide

**Date:** April 2026

---

## Overview

This workflow compares differentially expressed (DE) gene lists across several
analyses and visualizes where they overlap. It lets you:

1. Compare DE genes across multiple comparisons
2. Filter for a specific pathway (e.g. WNT signaling)
3. Visualize overlaps with UpSet plots and Venn diagrams
4. Extract genes from a chosen overlap pattern (shared by all, unique to one, or a custom subset)
5. Build targeted heatmaps for selected gene sets and samples

**You normally do not need to write any of this by hand.** The workflow is driven
by two notebooks — see [USER_GUIDE.md](../USER_GUIDE.md) for step-by-step
instructions:

- `run_upset_venn.Rmd` — overlap analysis (UpSet plots, Venn diagrams, overlap exports)
- `run_heatmap.Rmd` — flexible heatmaps of selected genes/samples

The rest of this document is a **reference for the underlying helper functions**,
for when you want to script a custom analysis directly.

---

## Function Reference

### Gene set operations (`scripts/utils/gene_set_operations.R`)

#### `load_de_gene_lists()`
Load significant-gene lists from one or more DE result files (`allGenes.xlsx`).

- `results_files` — named vector of paths to Excel/RDS results files
- `lfc_threshold` — log2 fold-change cutoff (default: 1)
- `padj_threshold` — adjusted p-value cutoff (default: 0.05)
- `direction` — `"up"`, `"down"`, or `"both"` (default: `"both"`)
- `gene_filter` — optional vector of genes to keep (e.g. a pathway)

#### `find_overlap_genes()`
Return genes matching an overlap pattern.

- `gene_lists` — named list of gene vectors
- `pattern` — `"all"` (in every list), `"any"` (in at least one), `"unique_to"`
  (only one list — also give `comparison_name`), or a vector like
  `c("List1", "List2")` for a specific combination
- `comparison_name` — required for the `"unique_to"` pattern

#### `generate_upset_plot()` / `generate_venn_diagram()`
Draw the overlaps. UpSet plots scale to any number of sets; Venn diagrams are
limited to 2–4 lists (larger inputs are skipped automatically).

#### `export_gene_lists()`
Write one Excel workbook with a sheet per list plus the all-way and every
pairwise overlap.

#### `load_pathway_genes()`
Read a pathway gene list (one column of symbols) from a file such as
`config/WntSignalingPathway.txt`.

### Heatmap helpers (`scripts/utils/heatmap_helpers.R`)

#### `generate_flexible_heatmap()`
Draw a heatmap with flexible gene/sample selection. Key parameters:

- `count_data` — normalized count matrix (a DE run's `normalized_counts.xlsx`)
- `gene_list` / `sample_list` — genes / samples to include (optional)
- `metadata` — sample metadata for annotations
- `cluster_rows` / `cluster_cols` — enable/disable clustering
- `sort_rows_alpha` — order gene rows alphabetically (overrides row clustering)
- `color_scheme` — `"viridis"`, `"blue_red"`, or `"blue_yellow"`
- `scale_method` — `"0to1"`, `"row"`, `"column"`, or `"none"`
- `annotation_col` — metadata columns to show as color bars

#### `get_samples_by_criteria()`
Select samples by metadata, e.g.
`get_samples_by_criteria(metadata, list(CellType = "GroupA"))`.

---

## Scripting Example

Finding WNT genes shared across three comparisons, then plotting them. This is
the same analysis the two notebooks perform, shown here as a script for
reference.

```r
source("scripts/utils/color_helpers.R")
source("scripts/utils/gene_set_operations.R")
source("scripts/utils/heatmap_helpers.R")

# 1. Load the pathway to filter on
wnt_genes <- load_pathway_genes("config/WntSignalingPathway.txt")

# 2. Load upregulated WNT genes from each comparison's allGenes.xlsx
comparisons <- c(
  "Comp1" = "figures/comparison1/allGenes.xlsx",
  "Comp2" = "figures/comparison2/allGenes.xlsx",
  "Comp3" = "figures/comparison3/allGenes.xlsx"
)
gene_lists <- load_de_gene_lists(
  results_files  = comparisons,
  lfc_threshold  = 1,
  padj_threshold = 0.05,
  direction      = "up",
  gene_filter    = wnt_genes
)

# 3. Visualize and export the overlaps
generate_upset_plot(gene_lists, "UpSet_WNT_up.pdf",
                    title = "Upregulated WNT genes across comparisons")
export_gene_lists(gene_lists, "WNT_up_overlaps.xlsx", include_overlaps = TRUE)

# 4. Pull out the genes shared by all three comparisons
universal_wnt <- find_overlap_genes(gene_lists, pattern = "all")

# 5. Heatmap of those shared genes
counts   <- load_normalized_counts("figures/combined/normalized_counts.xlsx")
metadata <- read.table("config/sample_metadata_full.txt",
                       header = TRUE, sep = "\t", row.names = 1)
generate_flexible_heatmap(
  count_data     = counts,
  gene_list      = universal_wnt,
  metadata       = metadata,
  output_file    = "Heatmap_Universal_WNT.pdf",
  title          = "WNT genes shared across all comparisons",
  color_scheme   = "viridis",
  scale_method   = "0to1",
  annotation_col = c("CellType", "Dataset")
)
```

---

## Tips

**UpSet vs Venn** — use a Venn diagram for a quick look at 2–4 lists; use an
UpSet plot for 3+ lists (especially 5+), where it stays readable.

**Heatmap scaling**
- `"0to1"` divides each gene by its maximum (highest sample = 1) — good for comparing patterns
- `"row"` is a per-gene z-score — good for emphasizing relative changes
- `"none"` shows the normalized counts as-is

**Fixed gene order** — turn clustering off, or set `sort_rows_alpha = TRUE`, when
you want a predictable (alphabetical) gene order instead of a clustered one.

---

## Output Files

A typical run produces:

- `UpSet_*.pdf` — all intersection patterns
- `Venn_*.pdf` — 2–4 way overlaps
- `gene_lists_with_overlaps.xlsx` — per-list sheets plus overlap sheets
- `gene_list_summary.xlsx` — summary counts
- one or more heatmap PDFs for the gene/sample subsets you choose

---

## Related Documentation

- [USER_GUIDE.md](../USER_GUIDE.md) — step-by-step notebook instructions
- `docs/PCA_DUAL_GROUPING.md` — PCA and correlation-heatmap QC
- `docs/ORGANISM_SELECTION.md` — organism database selection

---

**Last Updated:** April 2026
