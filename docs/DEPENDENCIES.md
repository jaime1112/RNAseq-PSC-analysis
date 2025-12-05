# Software Dependencies

This document lists all software dependencies required to run the RNAseq analysis pipeline.

## System Requirements

- **Operating System:** macOS, Linux, or Windows
- **R Version:** ≥ 4.2.0
- **RAM:** ≥ 8 GB recommended
- **Disk Space:** ≥ 10 GB for software and intermediate files

## R Installation

Download and install R from: https://cran.r-project.org/

### RStudio (Recommended)

Download and install RStudio from: https://posit.co/download/rstudio-desktop/

## R Package Dependencies

### Quick Installation

The easiest way to install all dependencies is to run the provided installation script:

```r
# From the project root directory
source("install_dependencies.R")
```

### Manual Installation

Alternatively, install packages manually in R console:

```r
# Install CRAN packages
install.packages(c(
  "tidyverse",
  "readxl",
  "writexl",
  "ggplot2",
  "ggrepel",
  "ggpubr",
  "pheatmap",
  "RColorBrewer",
  "VennDiagram",
  "UpSetR",
  "viridis",
  "rstudioapi"
))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "org.Hs.eg.db",
  "org.Ss.eg.db",
  "AnnotationDbi",
  "clusterProfiler",
  "enrichplot"
))
```

## Required Packages

### CRAN Packages

| Package | Purpose | Where Used |
|---------|---------|-----------|
| **tidyverse** | Data manipulation and visualization meta-package (includes dplyr, purrr, ggplot2) | All scripts |
| **readxl** | Read Excel files | run_heatmap.R, run_upset_venn.R |
| **writexl** | Export to Excel format | All main scripts |
| **ggplot2** | Grammar of graphics plotting (included in tidyverse) | Plotting utilities |
| **ggrepel** | Non-overlapping text labels in plots | MA plots |
| **ggpubr** | Publication-ready plots | MA plots, PCA |
| **pheatmap** | Heatmap generation | Heatmaps |
| **RColorBrewer** | Color palettes | Heatmaps |
| **viridis** | Perceptually uniform color scales | Heatmaps |
| **VennDiagram** | Venn diagram visualization | Gene set overlaps |
| **UpSetR** | UpSet plot visualization for multi-way comparisons | Gene set overlaps |
| **rstudioapi** | RStudio API for automatic path detection | All main scripts |

### Bioconductor Packages

| Package | Purpose | Where Used |
|---------|---------|-----------|
| **DESeq2** | RNA-seq differential expression analysis | Core analysis |
| **org.Hs.eg.db** | Human genome annotation database | Gene annotation, enrichment |
| **org.Ss.eg.db** | Pig genome annotation database | Gene annotation, enrichment |
| **AnnotationDbi** | Annotation database interface | Gene ID mapping |
| **clusterProfiler** | GO and KEGG enrichment analysis | Pathway enrichment |
| **enrichplot** | Enrichment result visualization | Enrichment plots |

### Package Dependencies (Auto-installed)

The following packages are automatically installed as dependencies of the packages above:

- **SummarizedExperiment** - Dependency of DESeq2
- **dplyr, purrr, readr, tidyr, stringr** - Included in tidyverse
- **grid** - Base R graphics (used for Venn diagrams)
- **stats** - Base R statistics

## Minimum Version Requirements

| Package | Minimum Version |
|---------|-----------------|
| R | 4.2.0 |
| DESeq2 | 1.40.0 |
| tidyverse | 2.0.0 |
| ggplot2 | 3.4.0 |
| clusterProfiler | 4.8.0 |
| org.Hs.eg.db | 3.17.0 |
| org.Ss.eg.db | 3.17.0 |

## Installation Notes

### rstudioapi

The `rstudioapi` package is required for automatic working directory detection when running scripts in RStudio. If you only run scripts from the command line, this package is optional but recommended.

### Organism Databases

- **org.Hs.eg.db** - Required for human samples or cross-species comparisons (using human as reference)
- **org.Ss.eg.db** - Required for pig-only analyses

If you only work with one organism, you can skip installing the other organism database.

## Troubleshooting

### Common Issues

**Issue:** Package installation fails with compilation errors

**Solution:** Install system dependencies (Ubuntu/Debian):
```bash
sudo apt-get install -y \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libpng-dev \
  libudunits2-dev
```

**Solution:** macOS (using Homebrew):
```bash
brew install curl openssl libxml2 libpng udunits
```

**Issue:** Bioconductor package version conflicts

**Solution:** Update all Bioconductor packages:
```r
BiocManager::install(ask = FALSE, update = TRUE)
```

**Issue:** `org.Ss.eg.db` not available for current R version

**Solution:** Install from source or use BiocManager:
```r
BiocManager::install("org.Ss.eg.db", force = TRUE)
```

**Issue:** tidyverse fails to install

**Solution:** Install individual tidyverse packages:
```r
install.packages(c("dplyr", "tidyr", "purrr", "readr", "ggplot2", "stringr"))
```

### Getting Help

- **Bioconductor Support:** https://support.bioconductor.org/
- **Stack Overflow:** Tag questions with `[r]` and `[deseq2]`
- **Package Documentation:** Use `?function_name` or `help(package="package_name")` in R console

## Verifying Installation

After installation, verify that all packages loaded successfully:

```r
# Test loading all packages
packages <- c(
  "tidyverse", "readxl", "writexl", "ggplot2", "ggrepel", "ggpubr",
  "pheatmap", "RColorBrewer", "viridis", "VennDiagram", "UpSetR",
  "rstudioapi", "DESeq2", "org.Hs.eg.db", "org.Ss.eg.db",
  "AnnotationDbi", "clusterProfiler", "enrichplot"
)

all_loaded <- sapply(packages, function(pkg) {
  result <- require(pkg, character.only = TRUE, quietly = TRUE)
  cat(sprintf("%-25s %s\n", pkg, ifelse(result, "✓", "✗")))
  return(result)
})

if (all(all_loaded)) {
  cat("\n✓ All packages loaded successfully!\n")
} else {
  cat("\n✗ Some packages failed to load. Please reinstall them.\n")
}
```

## Documenting Your Environment

To ensure reproducibility, save your session information:

```r
# Save session info to file
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

# Or view in console
sessionInfo()
```

Include this file with your analysis results.

## Updating Packages

To update all packages to the latest versions:

```r
# CRAN packages
update.packages(ask = FALSE)

# Bioconductor packages
BiocManager::install(ask = FALSE, update = TRUE)
```

**Note:** Major package updates may introduce breaking changes. Test your pipeline after updating.

## Alternative Installation Methods

### Using renv for Reproducibility

For exact version reproducibility:

```r
# Install renv
install.packages("renv")

# Initialize renv in project
renv::init()

# Install packages (using methods above)
# ...

# Create lockfile
renv::snapshot()

# On another machine, restore exact versions
renv::restore()
```

### Docker Container (Advanced)

For maximum reproducibility, a Docker container can be created. See the main README for details.

---

**Last Updated:** December 2025
