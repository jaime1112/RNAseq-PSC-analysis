---
editor_options: 
  markdown: 
    wrap: 72
---

# RNAseq Analysis Pipeline - User Guide

This guide provides step-by-step instructions for running the three main
analysis scripts in this pipeline.

## Table of Contents

1.  [Quick Start](#quick-start)
2.  [Script 1: Differential Expression
    Analysis](#script-1-differential-expression-analysis)
3.  [Script 2: UpSet Plots and Venn
    Diagrams](#script-2-upset-plots-and-venn-diagrams)
4.  [Script 3: Heatmap Generation](#script-3-heatmap-generation)
5.  [Configuration Files](#configuration-files)
6.  [Output Files](#output-files)
7.  [Troubleshooting](#troubleshooting)

------------------------------------------------------------------------

## Quick Start {#quick-start}

**Prerequisites:** - R (version 4.0 or higher) - **RStudio** (highly
recommended - download from
<https://posit.co/download/rstudio-desktop/>) - Git (optional, for
version control)

**RStudio Project Setup (Recommended):**

1.  In RStudio: `File` → `New Project` → `Existing Directory`
2.  Browse to your RNAseq folder
3.  Click "Create Project"
4.  This will create an `.Rproj` file and set the working directory
    automatically

**Installation:**

Open RStudio and run these commands in the Console:

``` r
# Install rstudioapi (for RStudio integration)
install.packages("rstudioapi")

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "org.Hs.eg.db", "org.Ss.eg.db",
                       "clusterProfiler", "enrichplot"))

# Install CRAN packages
install.packages(c("tidyverse", "readxl", "writexl", "ggplot2",
                   "ggrepel", "ggpubr", "pheatmap", "viridis",
                   "RColorBrewer", "VennDiagram", "UpSetR"))
```

------------------------------------------------------------------------

## Script 1: Differential Expression Analysis {#script-1-differential-expression-analysis}

**Purpose:** Run a complete differential expression analysis comparing
two groups (e.g., pPSC vs hPSC).

**File:** `run_differential_expression.R`

### How to Run in RStudio

1.  **Open the script** in RStudio: `File` → `Open File` → select
    `run_differential_expression.R`

2.  **Configure your analysis** by editing the parameters at the top of
    the script:

``` r
# Analysis name (used for output directory)
ANALYSIS_NAME <- "pPSC_vs_hPSC"

# Data configuration file
DATA_CONFIG <- "config/data_sources_full.csv"

# Sample metadata file
METADATA_FILE <- "config/sample_metadata_full.txt"

# Experimental design
BASELINE_LEVEL <- "hPSC"      # Reference group
COMPARISON_LEVEL <- "pPSC"    # Comparison group
DESIGN_FORMULA <- "~CellType" # DESeq2 design formula
FACTOR_NAME <- "CellType"     # Factor to test

# Organism for annotation and enrichment
# Options: "human"/"Hs"/"hsa" OR "pig"/"Ss"/"ssc"
ORGANISM <- "human"

# Analysis options
PERFORM_ENRICHMENT <- TRUE    # Run GO/KEGG enrichment
RUN_PATHWAY_ANALYSIS <- TRUE  # Filter pathway-specific genes
```

3.  **Run the script:** Click the **"Source"** button (top-right of
    script pane) or press `Ctrl/Cmd+Shift+S`

    > **Note:** You can also run it from Terminal with
    > `Rscript run_differential_expression.R` if preferred

4.  **Monitor progress:** Watch the R Console for progress messages and
    any warnings/errors

### What It Does

The script performs a complete analysis pipeline: 1. Loads and merges
count data from multiple sources 2. Creates DESeq2 dataset with your
experimental design 3. Runs differential expression analysis 4.
Generates PCA plots (by CellType and Dataset) 5. Extracts and annotates
results with gene symbols and Entrez IDs 6. Creates MA plots (overall
and pathway-specific) 7. (Optional) Performs GO/KEGG enrichment analysis
8. (Optional) Filters pathway-specific genes (e.g., WNT pathway)

### Output Files {#output-files}

All files are saved to `figures/[ANALYSIS_NAME]/`: - `allGenes.xlsx` -
Complete results for all genes - `sigGenes.xlsx` - Significant genes
only (padj \< 0.05) - `normalized_counts.xlsx` - DESeq2 normalized count
matrix - `PCA_by_CellType.pdf` - PCA plot colored by cell type -
`PCA_by_Dataset.pdf` - PCA plot colored by dataset - `MAplot.pdf` - Main
MA plot with top 20 genes labeled - `MAplot_WNTligands.pdf` - WNT
ligands highlighted - `MAplot_WNTreceptors.pdf` - WNT receptors
highlighted - `MAplot_Pluripotency.pdf` - Pluripotency genes
highlighted - `MAplot_WNTlinked.pdf` - WNT pathway genes labeled -
`GO_combined.xlsx` - GO enrichment results (if enabled) -
`KEGG_combined.xlsx` - KEGG enrichment results (if enabled) -
`GO_dotplot.pdf` - GO enrichment visualization (if enabled) -
`KEGG_dotplot.pdf` - KEGG enrichment visualization (if enabled) -
`WntGenes.xlsx` - WNT pathway genes only (if enabled)

------------------------------------------------------------------------

## Script 2: UpSet Plots and Venn Diagrams {#script-2-upset-plots-and-venn-diagrams}

**Purpose:** Compare gene lists across multiple comparisons or
categories to find overlaps.

**File:** `run_upset_venn.R`

### How to Run in RStudio

1.  **Open the script** in RStudio: `File` → `Open File` → select
    `run_upset_venn.R`

2.  **Choose your analysis mode** by editing the parameters:

**Method 1: Compare Multiple DE Results**

``` r
# Analysis name
ANALYSIS_NAME <- "multi_comparison"

# Method 1: Load from multiple comparison results
USE_MULTIPLE_COMPARISONS <- TRUE

# Paths to DE results files
COMPARISON_RESULTS <- c(
  "Comparison1" = "figures/comp1/allGenes.xlsx",
  "Comparison2" = "figures/comp2/allGenes.xlsx",
  "Comparison3" = "figures/comp3/allGenes.xlsx"
)

# Gene filtering
DIRECTION <- "both"  # "up", "down", or "both"
LFC_THRESHOLD <- 1
PADJ_THRESHOLD <- 0.05
```

**Method 2: Create Categories from Single Result**

``` r
# Analysis name
ANALYSIS_NAME <- "WNT_multi_comparison"

# Method 2: Load from single result file
USE_SINGLE_RESULT <- TRUE

# Path to single DE results file
SINGLE_RESULT_FILE <- "figures/full_test_pPSC_vs_hPSC/allGenes.xlsx"

# This will create categories:
# - High_Up (log2FC > 2)
# - Moderate_Up (1 < log2FC <= 2)
# - High_Down (log2FC < -2)
# - Moderate_Down (-2 <= log2FC < -1)
```

**Optional Pathway Filtering:**

``` r
# Pathway filtering (optional)
FILTER_PATHWAY <- TRUE
PATHWAY_FILE <- "config/WntSignalingPathway.txt"
```

3.  **Run the script:** Click the **"Source"** button or press
    `Ctrl/Cmd+Shift+S`

    > **Note:** You can also run from Terminal with
    > `Rscript run_upset_venn.R` if preferred

### What It Does

1.  Loads gene lists from either multiple comparisons or single result
2.  Filters genes by log2FC and adjusted p-value thresholds
3.  (Optional) Filters for pathway-specific genes
4.  Generates UpSet plot showing all overlaps
5.  Generates Venn diagram (if 2-4 lists)
6.  Finds and exports overlap genes
7.  Creates summary statistics

### Output Files

All files are saved to `figures/[ANALYSIS_NAME]/`: -
`UpSet_gene_overlaps.pdf` - UpSet plot showing all intersections -
`Venn_gene_overlaps.pdf` - Venn diagram (only if 2-4 lists) -
`gene_lists_with_overlaps.xlsx` - All gene lists plus overlap sets -
`gene_list_summary.xlsx` - Summary table with gene counts

------------------------------------------------------------------------

## Script 3: Heatmap Generation {#script-3-heatmap-generation}

**Purpose:** Create customizable heatmaps with flexible gene and sample
selection.

**File:** `run_heatmap.R`

### How to Run in RStudio

1.  **Open the script** in RStudio: `File` → `Open File` → select
    `run_heatmap.R`

2.  **Configure your heatmap** by editing the parameters:

``` r
# Analysis name (used for output filename)
ANALYSIS_NAME <- "WNT_heatmap"

# Output directory
OUTPUT_DIR <- "figures"

# Input files
COUNTS_FILE <- "figures/full_test_pPSC_vs_hPSC/normalized_counts.xlsx"
METADATA_FILE <- "config/sample_metadata_full.txt"
```

2.  **Choose gene selection method** (pick ONE):

``` r
# Option 1: Use all genes
USE_ALL_GENES <- FALSE

# Option 2: Use pathway genes
USE_PATHWAY_GENES <- TRUE
PATHWAY_FILE <- "config/WntSignalingPathway.txt"

# Option 3: Use specific gene list from file
USE_GENE_LIST_FILE <- FALSE
GENE_LIST_FILE <- "gene_list.txt"  # One gene per line

# Option 4: Use custom vector
USE_CUSTOM_GENES <- FALSE
CUSTOM_GENES <- c("WNT1", "WNT3A", "CTNNB1", "AXIN2")
```

3.  **Choose sample selection** (optional):

``` r
# Sample selection
# Leave NULL to use all samples
SAMPLE_CRITERIA <- NULL

# Examples:
# list(CellType = "hPSC")  # Only human samples
# list(CellType = "pPSC")  # Only pig samples
# list(Dataset = c("hESC_ours", "hESC_Selmi"))  # Specific datasets
```

4.  **Customize appearance:**

``` r
# Heatmap parameters
CLUSTER_ROWS <- TRUE
CLUSTER_COLS <- TRUE
COLOR_SCHEME <- "viridis"  # "viridis", "blue_red", or "blue_yellow"
SCALE_METHOD <- "0to1"     # "0to1", "row", "column", or "none"
ANNOTATION_COLS <- c("CellType", "Dataset")  # Metadata columns to display

# Plot dimensions
PLOT_WIDTH <- 12
PLOT_HEIGHT <- 10
```

5.  **Run the script:** Click the **"Source"** button or press
    `Ctrl/Cmd+Shift+S`

    > **Note:** You can also run from Terminal with
    > `Rscript run_heatmap.R` if preferred

### What It Does

1.  Loads normalized counts and metadata
2.  Selects genes based on your chosen method
3.  Filters samples based on criteria (if specified)
4.  Generates heatmap with clustering and annotations
5.  Saves as PDF

### Output Files

Output saved to `[OUTPUT_DIR]/`: - `[ANALYSIS_NAME].pdf` - Heatmap
visualization

------------------------------------------------------------------------

## Configuration Files {#configuration-files}

### Data Sources Configuration (`config/data_sources_full.csv`)

Format:

``` csv
dataset_name,filepath,remove_pattern
hESC_ours,data/hESC_ours/ReadsPerGene.txt,
hESC_Selmi,data/hESC_Selmi/counts.txt,ENSG[0-9]+
```

-   `dataset_name` - Name for this dataset
-   `filepath` - Path to count matrix file
-   `remove_pattern` - Optional regex pattern to clean gene IDs (leave
    empty if not needed)

### Sample Metadata (`config/sample_metadata_full.txt`)

Tab-delimited format:

```         
SampleID    CellType    Dataset    Replicate
Sample1     hPSC        hESC_ours  1
Sample2     hPSC        hESC_ours  2
Sample3     pPSC        pESC_ours  1
```

Requirements: - First column must be sample IDs (matching count matrix
columns) - Include any columns you want to use for grouping or
annotation - Tab-delimited with header row

### Pathway Gene Lists (`config/WntSignalingPathway.txt`)

Simple text file with one gene symbol per line:

```         
WNT1
WNT3A
CTNNB1
AXIN2
```

------------------------------------------------------------------------

## Output Files

### Differential Expression Results

**allGenes.xlsx** - Complete results table with columns: - `GeneID` -
Gene symbol - `baseMean` - Average normalized count across all samples -
`log2FoldChange` - Log2 fold change (comparison vs baseline) - `lfcSE` -
Standard error of log2FC - `pvalue` - Raw p-value - `padj` - Adjusted
p-value (FDR) - `symbol` - Gene symbol (annotated) - `entrezid` - Entrez
gene ID (for enrichment) - Sample columns - Normalized counts for each
sample

**sigGenes.xlsx** - Filtered to genes with padj \< 0.05

**normalized_counts.xlsx** - DESeq2 normalized count matrix

### Enrichment Results

**GO_combined.xlsx** / **KEGG_combined.xlsx** - Columns: - `Cluster` -
Gene set (Upregulated/Downregulated) - `ID` - Term ID - `Description` -
Term description - `GeneRatio` - Ratio of genes in term - `BgRatio` -
Background ratio - `pvalue` - Enrichment p-value - `p.adjust` - Adjusted
p-value - `geneID` - Gene symbols in term

------------------------------------------------------------------------

## Troubleshooting {#troubleshooting}

### Common Errors

**Error: Cannot find file** - Check that all file paths are correct
relative to the RNAseq directory - Verify data files exist in the
specified locations - Ensure config files are properly formatted

**Error: No samples found matching criteria** - Check that metadata
column names match exactly (case-sensitive) - Verify sample IDs in
metadata match count matrix column names - Print metadata to check
available values: `unique(metadata$CellType)`

**Error: Gene not found in annotation database** - Some genes may not
have Entrez IDs in the annotation database - This is normal and doesn't
affect the analysis - Check that you're using the correct organism
("human" vs "pig")

**Error: Factor levels not set correctly** - Ensure BASELINE_LEVEL and
COMPARISON_LEVEL match values in your metadata - Check factor column
name matches FACTOR_NAME - Values are case-sensitive

**Error: Insufficient memory** - For large datasets, increase R memory:
`R --max-mem-size=16G` - Consider filtering to pathway genes for
heatmaps - Use fewer samples or genes if necessary

### Getting Help

1.  Check that all required packages are installed
2.  Verify your configuration files are properly formatted
3.  Check the R console output for specific error messages
4.  Review the script parameters to ensure they match your data

------------------------------------------------------------------------

## Tips and Best Practices

### Differential Expression Analysis

1.  **Always check PCA plots first** - Ensure samples cluster as
    expected
2.  **Start with PERFORM_ENRICHMENT = FALSE** for faster initial runs
3.  **Use meaningful ANALYSIS_NAME** - Makes it easier to track multiple
    comparisons
4.  **Check sample size** - Need at least 2 replicates per group

### UpSet Plots

1.  **Method 1 (multiple comparisons)** - Best for comparing different
    contrasts
2.  **Method 2 (single result categories)** - Best for examining
    magnitude of changes
3.  **Pathway filtering** - Use when interested in specific biological
    processes
4.  **Venn diagrams** - Only work well with 2-4 lists; use UpSet for
    more

### Heatmaps

1.  **Gene selection** - Start with pathway genes for cleaner
    visualization
2.  **Scaling method**:
    -   `"0to1"` - Best for comparing expression patterns
    -   `"row"` - Best for seeing gene-level changes
    -   `"none"` - Best for absolute expression levels
3.  **Clustering** - Set to FALSE if you want to preserve specific
    ordering
4.  **Sample filtering** - Use to focus on specific comparisons

------------------------------------------------------------------------

## Example Workflows

### Workflow 1: Standard Two-Group Comparison

``` bash
# 1. Run differential expression
Rscript run_differential_expression.R

# 2. Create heatmap of WNT pathway genes
# Edit run_heatmap.R:
#   - COUNTS_FILE = "figures/pPSC_vs_hPSC/normalized_counts.xlsx"
#   - USE_PATHWAY_GENES = TRUE
Rscript run_heatmap.R
```

### Workflow 2: Multi-Comparison Analysis

``` bash
# 1. Run DE analysis for each comparison separately
# Edit run_differential_expression.R for each comparison:
#   - ANALYSIS_NAME = "comparison1"
#   - COMPARISON_LEVEL = "GroupA"
Rscript run_differential_expression.R

# Repeat for other comparisons...

# 2. Compare results with UpSet plot
# Edit run_upset_venn.R:
#   - USE_MULTIPLE_COMPARISONS = TRUE
#   - List all comparison result files
Rscript run_upset_venn.R
```

### Workflow 3: Pathway-Focused Analysis

``` bash
# 1. Run full DE analysis with pathway filtering
# Edit run_differential_expression.R:
#   - RUN_PATHWAY_ANALYSIS = TRUE
Rscript run_differential_expression.R

# 2. Create heatmap of pathway genes
# Edit run_heatmap.R:
#   - USE_PATHWAY_GENES = TRUE
#   - PATHWAY_FILE = "config/WntSignalingPathway.txt"
Rscript run_heatmap.R

# 3. Compare up vs down regulated pathway genes
# Edit run_upset_venn.R:
#   - USE_SINGLE_RESULT = TRUE
#   - FILTER_PATHWAY = TRUE
Rscript run_upset_venn.R
```

------------------------------------------------------------------------

## Quick Reference

| Task | Script | Key Parameter |
|-----------------|-------------------|------------------------------------|
| Run DE analysis | `run_differential_expression.R` | `ANALYSIS_NAME`, `ORGANISM` |
| Compare multiple analyses | `run_upset_venn.R` | `USE_MULTIPLE_COMPARISONS = TRUE` |
| Compare up vs down regulation | `run_upset_venn.R` | `USE_SINGLE_RESULT = TRUE` |
| Create pathway heatmap | `run_heatmap.R` | `USE_PATHWAY_GENES = TRUE` |
| Create custom gene heatmap | `run_heatmap.R` | `USE_CUSTOM_GENES = TRUE` |
| Filter samples in heatmap | `run_heatmap.R` | `SAMPLE_CRITERIA` |

------------------------------------------------------------------------

**Last Updated:** December 2025
