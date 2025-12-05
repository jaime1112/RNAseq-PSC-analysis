# Multi-Comparison Workflow Guide

**Date:** December 4, 2025
**Feature:** Flexible gene set comparison across multiple analyses with targeted visualizations

---

## Overview

This workflow enables you to:
1. Compare differentially expressed genes across multiple analyses
2. Filter for specific pathways (e.g., WNT signaling)
3. Visualize overlaps with UpSet plots and Venn diagrams
4. Extract genes from specific overlap patterns (universal, unique, custom)
5. Create targeted heatmaps for selected gene sets and samples

---

## Key Use Cases

### Use Case 1: Finding Universal WNT Genes Across Conditions

Compare upregulated WNT genes in 3 different experimental conditions to find genes that are consistently upregulated.

### Use Case 2: Condition-Specific Responses

Identify genes that are differentially expressed uniquely in one condition but not others.

### Use Case 3: Targeted Heatmaps

Create heatmaps showing expression patterns of specific gene subsets (e.g., genes common to all conditions) in selected samples.

---

## New Utility Functions

### Gene Set Operations (`scripts/utils/gene_set_operations.R`)

#### 1. `load_de_gene_lists()`
Load gene lists from multiple DE analysis results.

**Parameters:**
- `results_files` - Named vector of paths to Excel/RDS results files
- `lfc_threshold` - Log2 fold change cutoff (default: 1)
- `padj_threshold` - Adjusted p-value cutoff (default: 0.05)
- `direction` - "up", "down", or "both" (default: "both")
- `gene_filter` - Optional vector of genes to filter for

**Example:**
```r
# Compare 3 different analyses, get upregulated WNT genes
wnt_genes <- load_pathway_genes("config/WntSignalingPathway.txt")

gene_lists <- load_de_gene_lists(
  results_files = c(
    "Condition_A" = "figures/conditionA/allGenes.xlsx",
    "Condition_B" = "figures/conditionB/allGenes.xlsx",
    "Condition_C" = "figures/conditionC/allGenes.xlsx"
  ),
  lfc_threshold = 1,
  padj_threshold = 0.05,
  direction = "up",
  gene_filter = wnt_genes
)
```

#### 2. `find_overlap_genes()`
Find genes in specific overlap patterns.

**Parameters:**
- `gene_lists` - Named list of gene vectors
- `pattern` - Overlap pattern to find:
  - `"all"` - Genes in ALL lists
  - `"any"` - Genes in ANY list
  - `"unique_to"` - Genes unique to one list (specify `comparison_name`)
  - `c("List1", "List2")` - Genes in specific combination
- `comparison_name` - For "unique_to" pattern

**Examples:**
```r
# Find genes in all comparisons
universal_genes <- find_overlap_genes(gene_lists, pattern = "all")

# Find genes unique to Condition_A
unique_to_a <- find_overlap_genes(
  gene_lists,
  pattern = "unique_to",
  comparison_name = "Condition_A"
)

# Find genes in both A and B but not necessarily C
a_and_b <- find_overlap_genes(
  gene_lists,
  pattern = c("Condition_A", "Condition_B")
)
```

#### 3. `generate_upset_plot()`
Create UpSet plot showing all overlaps.

```r
generate_upset_plot(
  gene_lists = gene_lists,
  output_file = "UpSet_comparisons.pdf",
  title = "Gene Set Overlaps",
  n_intersects = 20,
  width = 12,
  height = 7
)
```

#### 4. `generate_venn_diagram()`
Create Venn diagram (2-4 way comparisons).

```r
# Works with 2-4 lists
generate_venn_diagram(
  gene_lists = gene_lists[1:3],  # First 3 comparisons
  output_file = "Venn_comparisons.pdf",
  title = "3-way Comparison",
  width = 10,
  height = 10
)
```

#### 5. `export_gene_lists()`
Export gene lists to Excel with overlap information.

```r
export_gene_lists(
  gene_lists = gene_lists,
  output_file = "gene_lists_with_overlaps.xlsx",
  include_overlaps = TRUE
)
# Creates sheets for each list plus "All_Overlaps" and pairwise overlaps
```

### Heatmap Helpers (`scripts/utils/heatmap_helpers.R`)

#### 1. `generate_flexible_heatmap()`
Create heatmap with gene and sample selection.

**Key Parameters:**
- `count_data` - Normalized count matrix
- `gene_list` - Vector of genes to include (optional)
- `sample_list` - Vector of samples to include (optional)
- `metadata` - Sample metadata for annotations
- `cluster_rows/cluster_cols` - Enable/disable clustering
- `color_scheme` - "viridis", "blue_red", "blue_yellow"
- `scale_method` - "0to1", "row", "column", "none"
- `annotation_col` - Metadata columns to display

**Example:**
```r
# Heatmap of universal genes in human samples only
generate_flexible_heatmap(
  count_data = normalized_counts,
  gene_list = universal_genes,
  sample_list = human_samples,
  metadata = metadata,
  output_file = "Heatmap_Universal_Genes_Human.pdf",
  title = "Universal Genes in Human Samples",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color_scheme = "viridis",
  scale_method = "0to1",
  annotation_col = c("CellType", "Dataset"),
  width = 10,
  height = 8
)
```

#### 2. `get_samples_by_criteria()`
Select samples based on metadata.

```r
# Get only pig samples
pig_samples <- get_samples_by_criteria(
  metadata = metadata,
  criteria = list(CellType = "pPSC")
)

# Get specific datasets
dataset_samples <- get_samples_by_criteria(
  metadata = metadata,
  criteria = list(Dataset = c("piPSC_ours", "pESC_Choi2024"))
)
```

#### 3. `generate_multiple_heatmaps()`
Batch create heatmaps for gene list collection.

```r
# Create one heatmap per gene list
generate_multiple_heatmaps(
  count_data = counts,
  gene_list_collection = gene_lists,
  output_dir = "figures/heatmaps",
  base_filename = "Heatmap_Category",
  metadata = metadata,
  annotation_col = c("CellType", "Dataset"),
  color_scheme = "viridis",
  scale_method = "0to1"
)
```

---

## Complete Workflow Example

### Scenario: Compare WNT Genes Across 3 Conditions

```r
# === SETUP ===
library(tidyverse)
library(VennDiagram)
library(UpSetR)
library(pheatmap)

source("scripts/utils/gene_set_operations.R")
source("scripts/utils/heatmap_helpers.R")

# === 1. LOAD PATHWAY GENES ===
wnt_genes <- load_pathway_genes(
  "config/WntSignalingPathway.txt",
  gene_column = "GeneID"
)

# === 2. LOAD DE RESULTS FROM MULTIPLE COMPARISONS ===
# Method A: From different analysis results
comparisons <- c(
  "Treatment_vs_Control" = "figures/treatment_vs_control/allGenes.xlsx",
  "TimePoint1_vs_Baseline" = "figures/timepoint1/allGenes.xlsx",
  "TimePoint2_vs_Baseline" = "figures/timepoint2/allGenes.xlsx"
)

gene_lists_up <- load_de_gene_lists(
  results_files = comparisons,
  lfc_threshold = 1,
  padj_threshold = 0.05,
  direction = "up",
  gene_filter = wnt_genes
)

gene_lists_down <- load_de_gene_lists(
  results_files = comparisons,
  lfc_threshold = 1,
  padj_threshold = 0.05,
  direction = "down",
  gene_filter = wnt_genes
)

# Method B: From single analysis with categories
# (see run_upset_venn.R for implementation)

# === 3. VISUALIZE OVERLAPS ===
# UpSet plot
generate_upset_plot(
  gene_lists = gene_lists_up,
  output_file = "UpSet_WNT_upregulated.pdf",
  title = "Upregulated WNT Genes Across Conditions"
)

# Venn diagram
generate_venn_diagram(
  gene_lists = gene_lists_up,
  output_file = "Venn_WNT_upregulated.pdf",
  title = "Upregulated WNT Genes"
)

# === 4. EXTRACT SPECIFIC OVERLAPS ===
# Universal genes (in all conditions)
universal_up <- find_overlap_genes(gene_lists_up, pattern = "all")

# Unique to treatment
unique_treatment <- find_overlap_genes(
  gene_lists_up,
  pattern = "unique_to",
  comparison_name = "Treatment_vs_Control"
)

# In treatment AND timepoint1
treatment_time1 <- find_overlap_genes(
  gene_lists_up,
  pattern = c("Treatment_vs_Control", "TimePoint1_vs_Baseline")
)

# === 5. EXPORT GENE LISTS ===
export_gene_lists(
  gene_lists = gene_lists_up,
  output_file = "WNT_upregulated_all_comparisons.xlsx",
  include_overlaps = TRUE
)

# === 6. LOAD NORMALIZED COUNTS ===
counts <- load_normalized_counts("figures/combined/normalized_counts.xlsx")
metadata <- read.table("config/sample_metadata.txt", header = TRUE, row.names = 1)

# === 7. CREATE TARGETED HEATMAPS ===
# Heatmap 1: All upregulated WNT genes
all_up_wnt <- unique(unlist(gene_lists_up))
generate_flexible_heatmap(
  count_data = counts,
  gene_list = all_up_wnt,
  metadata = metadata,
  output_file = "Heatmap_All_Up_WNT.pdf",
  title = "All Upregulated WNT Genes",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color_scheme = "viridis",
  annotation_col = c("CellType", "Treatment")
)

# Heatmap 2: Universal genes only
if (length(universal_up) > 0) {
  generate_flexible_heatmap(
    count_data = counts,
    gene_list = universal_up,
    metadata = metadata,
    output_file = "Heatmap_Universal_Up_WNT.pdf",
    title = "Universally Upregulated WNT Genes",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color_scheme = "blue_red",
    scale_method = "row"  # Z-score for better visualization
  )
}

# Heatmap 3: Treatment-unique genes in treated samples only
treated_samples <- get_samples_by_criteria(
  metadata = metadata,
  criteria = list(Treatment = "Treated")
)

generate_flexible_heatmap(
  count_data = counts,
  gene_list = unique_treatment,
  sample_list = treated_samples,
  metadata = metadata,
  output_file = "Heatmap_Treatment_Unique_TreatedSamples.pdf",
  title = "Treatment-Unique Genes in Treated Samples"
)
```

---

##Real-World Example: Original Workflow

Your original workflow of finding WNT genes across 3 comparisons:

```r
# === LOAD 3 COMPARISON RESULTS ===
comparisons <- c(
  "hPSC_vs_pPSC" = "figures/hPSC_vs_pPSC/allGenes.xlsx",
  "hiPSC_vs_piPSC" = "figures/hiPSC_vs_piPSC/allGenes.xlsx",
  "hESC_vs_pESC" = "figures/hESC_vs_pESC/allGenes.xlsx"
)

# Load WNT pathway genes
wnt_genes <- load_pathway_genes("config/WntSignalingPathway.txt")

# Get upregulated WNT genes from each comparison
gene_lists <- load_de_gene_lists(
  results_files = comparisons,
  lfc_threshold = 1,
  padj_threshold = 0.05,
  direction = "up",
  gene_filter = wnt_genes
)

# === VISUALIZE WITH UPSET PLOT ===
generate_upset_plot(
  gene_lists = gene_lists,
  output_file = "UpSet_WNT_Upregulated_3Comparisons.pdf",
  title = "Upregulated WNT Genes: 3 Cross-Species Comparisons"
)

# === FIND UNIVERSAL GENES ===
universal_wnt_up <- find_overlap_genes(gene_lists, pattern = "all")
cat("Universal upregulated WNT genes:", length(universal_wnt_up), "\n")
print(universal_wnt_up)

# === SAVE UNIVERSAL GENES ===
write.table(
  data.frame(Gene = universal_wnt_up),
  "Universal_Upregulated_WNT_Genes.txt",
  row.names = FALSE,
  quote = FALSE
)

# === CREATE HEATMAP OF UNIVERSAL GENES ===
counts <- load_normalized_counts("figures/combined/normalized_counts.xlsx")
metadata <- read.table("config/sample_metadata.txt", header = TRUE, row.names = 1)

generate_flexible_heatmap(
  count_data = counts,
  gene_list = universal_wnt_up,
  metadata = metadata,
  output_file = "Heatmap_Universal_WNT_Upregulated.pdf",
  title = "Universal WNT Genes (Upregulated in All Comparisons)",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color_scheme = "viridis",
  scale_method = "0to1",
  annotation_col = c("CellType", "Species"),
  width = 12,
  height = 10
)
```

---

## Advanced Features

### 1. Multiple Pathway Analysis

Compare different pathways across conditions:

```r
# Load multiple pathways
pathways <- list(
  WNT = load_pathway_genes("config/WntSignalingPathway.txt"),
  Pluripotency = c("POU5F1", "NANOG", "SOX2", ...),
  Metabolism = load_pathway_genes("config/MetabolicPathway.txt")
)

# Create gene lists for each pathway
for (pathway_name in names(pathways)) {
  gene_lists[[pathway_name]] <- load_de_gene_lists(
    results_files = comparisons,
    direction = "up",
    gene_filter = pathways[[pathway_name]]
  )
}
```

### 2. Sample Subset Comparisons

Create heatmaps for different sample subsets:

```r
# Define sample groups
sample_groups <- list(
  Human = get_samples_by_criteria(metadata, list(Species = "Human")),
  Pig = get_samples_by_criteria(metadata, list(Species = "Pig")),
  PSC = get_samples_by_criteria(metadata, list(CellType = "PSC")),
  ESC = get_samples_by_criteria(metadata, list(CellType = "ESC"))
)

# Create heatmap for each sample group
for (group_name in names(sample_groups)) {
  generate_flexible_heatmap(
    count_data = counts,
    gene_list = universal_genes,
    sample_list = sample_groups[[group_name]],
    output_file = paste0("Heatmap_Universal_", group_name, ".pdf"),
    title = paste("Universal Genes in", group_name, "Samples")
  )
}
```

### 3. Regulation Direction Comparison

Compare up vs down regulation patterns:

```r
# Get both up and down lists
gene_lists_combined <- c(
  load_de_gene_lists(comparisons, direction = "up", gene_filter = wnt_genes),
  load_de_gene_lists(comparisons, direction = "down", gene_filter = wnt_genes)
)

# Rename for clarity
names(gene_lists_combined) <- c(
  paste0(names(comparisons), "_Up"),
  paste0(names(comparisons), "_Down")
)

# Create comprehensive UpSet plot
generate_upset_plot(
  gene_lists = gene_lists_combined,
  output_file = "UpSet_WNT_UpDown_AllComparisons.pdf",
  n_intersects = 30
)
```

---

## Tips and Best Practices

### 1. Choosing UpSet vs Venn

- **Use Venn diagrams for:**
  - 2-4 comparisons
  - Simple overlap visualization
  - Quick visual assessment

- **Use UpSet plots for:**
  - 3+ comparisons (especially 5+)
  - Detailed overlap patterns
  - Quantitative comparison

### 2. Heatmap Scaling

Different scaling methods for different purposes:

- **`scale_method = "0to1"`**: Best for comparing expression levels across genes
- **`scale_method = "row"`**: Best for comparing expression patterns (Z-score)
- **`scale_method = "column"`**: Best for sample-to-sample comparison
- **`scale_method = "none"`**: Raw normalized counts

### 3. Clustering Decisions

- **Cluster rows (genes):** Shows which genes have similar expression patterns
- **Cluster cols (samples):** Shows which samples are similar
- **No clustering:** Preserves original gene/sample order (useful for comparing specific orderings)

### 4. Gene List Size Considerations

- **Large gene lists (>100 genes):** Consider reducing font size or filtering to top genes
- **Small gene lists (<10 genes):** May not need clustering
- **Very different sizes:** Consider normalizing or using log-scale visualizations

---

## Output Files

A typical multi-comparison workflow creates:

**Overlap Visualizations:**
- `UpSet_*.pdf` - Shows all intersection patterns
- `Venn_*.pdf` - 2-4 way visual overlaps

**Gene Lists:**
- `*_gene_lists_with_overlaps.xlsx` - All lists + overlap sheets
- `gene_list_summary.xlsx` - Statistics table

**Heatmaps:**
- `Heatmap_All_*.pdf` - All genes combined
- `Heatmap_Universal_*.pdf` - Only universal genes
- `Heatmap_Unique_*.pdf` - Condition-specific genes
- `Heatmap_Category_*.pdf` - One per gene list

---

## Implementation Scripts

The workflows described in this document are implemented in the main pipeline scripts:
- `run_upset_venn.R` - For gene set overlap analysis (UpSet plots and Venn diagrams)
- `run_heatmap.R` - For flexible heatmap generation with gene/sample filtering

Refer to the main [README.md](../README.md) for usage instructions.

---

## Related Documentation

- `docs/PCA_DUAL_GROUPING.md` - PCA visualization options
- `docs/ORGANISM_SELECTION.md` - Organism database selection
- `README.md` - Main project documentation

---

**Status:** ✅ Tested and Production-Ready
**Last Updated:** December 4, 2025
