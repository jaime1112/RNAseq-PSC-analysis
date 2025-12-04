# Dual PCA Plot Grouping Feature

**Date:** December 4, 2025
**Feature:** Multiple PCA plot groupings

---

## Overview

The pipeline now supports generating PCA plots with different sample groupings. This allows you to visualize your data at different levels of granularity:

1. **By CellType** - Broad comparison groups (e.g., hPSC vs pPSC)
2. **By Dataset** - Individual dataset sources (e.g., hESC_ours, hESC_Selmi, piPSC_ours, pESC_Choi2024)

---

## Why This Matters

### Broad Grouping (CellType)
Shows the main biological comparison you're testing. This is useful for:
- Seeing if your experimental groups separate clearly
- Understanding the overall experimental effect
- Publication figures showing the primary comparison

### Fine Grouping (Dataset)
Shows how individual datasets cluster. This is useful for:
- Quality control - checking if technical/batch effects exist
- Understanding dataset-specific variation
- Verifying that datasets from the same condition cluster together
- Identifying outlier datasets

---

## How to Use

### Method 1: Single PCA Plot

Generate one PCA plot with your chosen grouping:

```r
# Load utilities
source("scripts/utils/plotting_themes.R")

# By CellType
generate_pca_plot(
  rld = rld,
  output_file = "PCA_by_CellType.pdf",
  intgroup = "CellType"
)

# By Dataset
generate_pca_plot(
  rld = rld,
  output_file = "PCA_by_Dataset.pdf",
  intgroup = "Dataset"
)
```

### Method 2: Generate All PCA Plots Automatically

The `generate_pca_plots_multi()` function will create PCA plots for all available metadata columns:

```r
# Automatically generate PCA for all grouping factors
pca_files <- generate_pca_plots_multi(
  rld = rld,
  output_dir = "figures/my_analysis",
  base_filename = "PCA_plot",
  width = 8,
  height = 6
)
```

This will create:
- `PCA_plot_by_CellType.pdf`
- `PCA_plot_by_Dataset.pdf`
- Any other metadata columns you've added

---

## Configuration: Metadata File

Your metadata file must include the grouping columns. Example:

**config/sample_metadata_full.txt:**
```
Sample	CellType	Dataset
H1_rep1	hPSC	hESC_ours
H1_rep2	hPSC	hESC_ours
Selmi_H9.WT.1.1_RNA-seq	hPSC	hESC_Selmi
PiPSC8_rep1	pPSC	piPSC_ours
Choi_PA_1_pESC_FIW_rep1	pPSC	pESC_Choi2024
```

**Required columns:**
- `CellType` - Broad experimental groups (used in DESeq2 design)
- `Dataset` - Source of each sample

**Optional columns:**
You can add any additional metadata columns (e.g., `Batch`, `Replicate`, `Treatment`) and the multi-function will generate PCA plots for all of them.

---

## Test Script

A complete test script is provided: `test_dual_pca.R`

Run it to verify the functionality:

```bash
Rscript test_dual_pca.R
```

**Expected outputs:**
- `PCA_plot_by_CellType.pdf` - Shows 2 groups (hPSC, pPSC)
- `PCA_plot_by_Dataset.pdf` - Shows 4 groups (hESC_ours, hESC_Selmi, piPSC_ours, pESC_Choi2024)

---

## Updated Analysis Scripts

### full_test_analysis.R

Now automatically generates both PCA plots:

```r
# Step 4: Generate PCA Plots
cat("=== Step 4: Generating PCA Plots ===\n")

# By CellType (broad groups)
generate_pca_plot(
  rld = rld,
  output_file = file.path(output_dir, "PCA_plot_by_CellType.pdf"),
  intgroup = "CellType"
)

# By Dataset (individual datasets)
generate_pca_plot(
  rld = rld,
  output_file = file.path(output_dir, "PCA_plot_by_Dataset.pdf"),
  intgroup = "Dataset"
)

# Backwards compatibility
generate_pca_plot(
  rld = rld,
  output_file = file.path(output_dir, "PCA_plot.pdf"),
  intgroup = "CellType"
)
```

---

## Function Reference

### `generate_pca_plot()`

Generate a single PCA plot with specified grouping.

**Parameters:**
- `rld` - Variance-stabilized DESeq2 object
- `output_file` - Path to output PDF
- `intgroup` - Metadata column to color by (default: "CellType")
- `width` - Plot width in inches (default: 8)
- `height` - Plot height in inches (default: 6)

**Example:**
```r
generate_pca_plot(
  rld = rld,
  output_file = "my_pca.pdf",
  intgroup = "Dataset",
  width = 10,
  height = 8
)
```

### `generate_pca_plots_multi()`

Generate PCA plots for all available metadata columns.

**Parameters:**
- `rld` - Variance-stabilized DESeq2 object
- `output_dir` - Directory for output files
- `base_filename` - Base name for files (default: "PCA_plot")
- `width` - Plot width in inches (default: 8)
- `height` - Plot height in inches (default: 6)

**Returns:**
- Named list of output file paths

**Example:**
```r
pca_files <- generate_pca_plots_multi(
  rld = rld,
  output_dir = "figures/my_analysis",
  base_filename = "PCA"
)

# Access specific file paths
print(pca_files$CellType)
print(pca_files$Dataset)
```

---

## Tips

### Adding Custom Metadata Columns

You can add any metadata columns to your sample metadata file:

```
Sample	CellType	Dataset	Batch	Passage
H1_rep1	hPSC	hESC_ours	Batch1	P20
H1_rep2	hPSC	hESC_ours	Batch1	P20
H1_rep3	hPSC	hESC_ours	Batch2	P22
```

Then generate PCA plots for each:

```r
# By Batch
generate_pca_plot(rld, "PCA_by_Batch.pdf", intgroup = "Batch")

# By Passage
generate_pca_plot(rld, "PCA_by_Passage.pdf", intgroup = "Passage")

# Or use multi-function to get all
pca_files <- generate_pca_plots_multi(rld, output_dir)
```

### Interpreting PCA Plots

**CellType PCA:**
- Should show clear separation between experimental groups
- Indicates if your biological condition has a strong effect
- Used for publication main figures

**Dataset PCA:**
- Samples from the same dataset should cluster together
- Datasets from the same CellType should be near each other
- Large separation between datasets from same CellType suggests batch effects
- Useful for quality control and troubleshooting

### Batch Effect Detection

If datasets from the same CellType don't cluster near each other in the Dataset PCA, you may have batch effects. Consider:

1. Adding batch correction (e.g., using `limma::removeBatchEffect()`)
2. Including batch in your DESeq2 design formula: `~ Batch + CellType`
3. Using ComBat or similar batch correction methods

---

## File Locations

**Utility functions:**
- `scripts/utils/plotting_themes.R` - Contains PCA plotting functions

**Test scripts:**
- `test_dual_pca.R` - Standalone test of dual PCA functionality
- `full_test_analysis.R` - Full pipeline including dual PCA

**Configuration:**
- `config/sample_metadata_full.txt` - Sample metadata with CellType and Dataset

**Documentation:**
- `docs/PCA_DUAL_GROUPING.md` - This file

---

## Backwards Compatibility

The original `generate_pca_plot()` function still works exactly as before. The default grouping is by "CellType".

Existing code will continue to work without changes:

```r
# This still works
generate_pca_plot(rld, "PCA_plot.pdf")
```

The `full_test_analysis.R` script also maintains backwards compatibility by creating the original `PCA_plot.pdf` file.

---

## Examples

### Example 1: Basic Analysis

```r
# Simple two-group comparison
dds <- DESeqDataSetFromMatrix(counts, metadata, ~CellType)
dds <- DESeq(dds)
rld <- vst(dds)

# Generate both PCA plots
generate_pca_plot(rld, "PCA_CellType.pdf", intgroup = "CellType")
generate_pca_plot(rld, "PCA_Dataset.pdf", intgroup = "Dataset")
```

### Example 2: Batch Experiment

```r
# Multi-batch experiment
metadata <- read.table("metadata.txt", header = TRUE, row.names = 1)
# Columns: CellType, Dataset, Batch

dds <- DESeqDataSetFromMatrix(counts, metadata, ~ Batch + CellType)
dds <- DESeq(dds)
rld <- vst(dds)

# Generate PCA for all groupings
pca_files <- generate_pca_plots_multi(rld, "figures/batch_analysis")

# Will create:
# - PCA_plot_by_CellType.pdf
# - PCA_plot_by_Dataset.pdf
# - PCA_plot_by_Batch.pdf
```

### Example 3: Custom Analysis

```r
# Only generate Dataset PCA for QC
generate_pca_plot(
  rld = rld,
  output_file = "QC_dataset_clustering.pdf",
  intgroup = "Dataset",
  width = 12,
  height = 10
)
```

---

## Version History

**Version 1.0** (December 4, 2025)
- Initial implementation of dual PCA grouping
- Added `generate_pca_plots_multi()` function
- Updated metadata files with Dataset column
- Created test scripts and documentation
- Maintained backwards compatibility

---

## Related Documentation

- `README.md` - Main project documentation
- `METHODS.md` - Computational methods
- `DEPENDENCIES.md` - Software requirements
- `COMPLETE_TEST_RESULTS.md` - Testing documentation

---

**Status:** ✅ Tested and Production-Ready
**Last Updated:** December 4, 2025
