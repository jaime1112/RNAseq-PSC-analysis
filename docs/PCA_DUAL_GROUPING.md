# PCA Plots and Sample Correlation Heatmap

## Overview

During the QC step of the DE pipeline (Section 6 of `run_differential_expression.Rmd`), the script generates:

1. **One PCA plot per metadata column** (`PCA_by_[COLUMN].pdf`) — the pipeline loops over every column in `colData` (e.g. `PCA_by_CellType.pdf`, `PCA_by_Dataset.pdf`), so you can inspect clustering by biological group, data source, sex, batch, etc.
2. **Sample correlation heatmap** (`sample_correlation_heatmap.pdf`) — Pearson correlation matrix with hierarchical clustering, annotated by all metadata columns

## Sample Correlation Heatmap Details

- **Color scale:** white-to-blue (darker blue = higher correlation)
- **Clustering:** Ward's linkage on `1 - correlation` distance
- **Annotations:** All metadata columns from `colData`, ordered by number of unique values (most groups on top, fewest closest to heatmap)
- **Method:** Pearson by default, configurable via `method` parameter

## Function Reference

### `generate_pca_plot(rld, output_file, intgroup, width, height)`

Generate a single PCA plot. `intgroup` specifies the metadata column to color by.

### `generate_correlation_heatmap(rld, output_file, method, width, height)`

Generate sample-to-sample correlation heatmap. Automatically includes all metadata columns as annotations.

## Interpreting Results

- **PCA by analysis factor:** Groups should separate clearly if the biological effect is strong
- **PCA by Dataset:** Datasets from the same condition should cluster nearby; large separation suggests batch effects
- **Correlation heatmap:** High within-group correlation = consistent replicates; unexpected low correlation = potential outlier

---

**Last Updated:** April 2026
