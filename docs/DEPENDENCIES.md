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

### Installation Script

Save and run the following script to install all dependencies:

```r
# install_dependencies.R

# Install CRAN packages
cran_packages <- c(
  # Data manipulation
  "tidyverse",      # Meta-package: dplyr, ggplot2, tidyr, readr, etc.
  "writexl",        # Export to Excel format

  # Visualization
  "ggplot2",        # Grammar of graphics plotting
  "ggrepel",        # Non-overlapping text labels
  "gplots",         # Additional plotting tools
  "vegan",          # Ecological statistics and visualization
  "pheatmap",       # Pretty heatmaps
  "RColorBrewer",   # Color palettes
  "FactoMineR",     # Multivariate analysis
  "ggpubr",         # Publication-ready plots
  "ggsci",          # Scientific journal color palettes
  "gghighlight",    # Highlight specific data in plots
  "VennDiagram",    # Venn diagrams
  "UpSetR",         # UpSet plots for set intersections
  "viridis",        # Perceptually uniform color scales
  "scico",          # Scientific color maps

  # R Markdown
  "knitr",          # Dynamic report generation
  "rmarkdown"       # R Markdown document rendering
)

# Install CRAN packages
install.packages(cran_packages)

# Install Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
bioc_packages <- c(
  # Differential expression
  "DESeq2",         # RNA-seq differential expression analysis
  "tximport",       # Import transcript-level counts
  "apeglm",         # Approximate posterior estimation for GLM

  # Annotation
  "org.Hs.eg.db",   # Human genome annotation
  "org.Ss.eg.db",   # Pig genome annotation
  "AnnotationDbi",  # Annotation database interface

  # Enrichment analysis
  "clusterProfiler", # GO/KEGG enrichment analysis
  "enrichplot",      # Enrichment visualization
  "DOSE",            # Disease ontology semantic and enrichment

  # General Bioconductor utilities
  "BiocGenerics",    # S4 generic functions
  "S4Vectors"        # S4 vector classes
)

BiocManager::install(bioc_packages)

# Verify installation
cat("\\n=== Verifying Installation ===\\n")
all_packages <- c(cran_packages, bioc_packages)
installed <- sapply(all_packages, function(pkg) {
  result <- require(pkg, character.only = TRUE, quietly = TRUE)
  cat(sprintf("%-20s %s\\n", pkg, ifelse(result, "✓", "✗")))
  return(result)
})

if (all(installed)) {
  cat("\\n✓ All packages installed successfully!\\n")
} else {
  failed <- names(installed)[!installed]
  cat("\\n✗ Failed to install:", paste(failed, collapse = ", "), "\\n")
}
```

### Core Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| DESeq2 | ≥ 1.40.0 | Differential expression analysis |
| apeglm | ≥ 1.22.0 | Log fold change shrinkage |
| tidyverse | ≥ 2.0.0 | Data manipulation and visualization |
| ggplot2 | ≥ 3.4.0 | Graphics and plotting |
| pheatmap | ≥ 1.0.12 | Heatmap generation |
| clusterProfiler | ≥ 4.8.0 | Enrichment analysis |
| org.Hs.eg.db | ≥ 3.17.0 | Human gene annotations |
| org.Ss.eg.db | ≥ 3.17.0 | Pig gene annotations |

### Complete Dependency List

#### CRAN Packages

```r
# Data Manipulation
- tidyverse (2.0.0)
- dplyr (1.1.0)
- tidyr (1.3.0)
- purrr (1.0.1)
- readr (2.1.4)

# Export
- writexl (1.4.2)

# Visualization - Core
- ggplot2 (3.4.3)
- ggrepel (0.9.3)
- scales (1.2.1)

# Visualization - Specialized
- gplots (3.1.3)
- vegan (2.6-4)
- pheatmap (1.0.12)
- RColorBrewer (1.1-3)
- viridis (0.6.3)
- scico (1.5.0)

# Visualization - Publication
- ggpubr (0.6.0)
- ggsci (3.0.0)
- gghighlight (0.4.1)

# Set Visualization
- VennDiagram (1.7.3)
- UpSetR (1.4.0)

# Statistical Analysis
- FactoMineR (2.8)

# Reporting
- knitr (1.43)
- rmarkdown (2.23)
```

#### Bioconductor Packages

```r
# Core Bioconductor
- BiocManager (1.30.21)
- BiocGenerics (0.46.0)
- S4Vectors (0.38.0)

# Differential Expression
- DESeq2 (1.40.0)
- tximport (1.28.0)
- apeglm (1.22.0)
- SummarizedExperiment (1.30.0)

# Annotation
- AnnotationDbi (1.62.0)
- org.Hs.eg.db (3.17.0)
- org.Ss.eg.db (3.17.0)

# Enrichment
- clusterProfiler (4.8.0)
- enrichplot (1.20.0)
- DOSE (3.26.0)
- GO.db (3.17.0)
- KEGGREST (1.40.0)
```

## Alternative: Using renv

For reproducibility, consider using `renv` for package management:

```r
# Install renv
install.packages("renv")

# Initialize renv in project
renv::init()

# Install packages (from the script above)
# ...

# Create lockfile
renv::snapshot()

# On another machine, restore packages
renv::restore()
```

The `renv.lock` file can be committed to version control for exact version reproducibility.

## Docker Container (Optional)

For maximum reproducibility, a Docker container can be created:

```dockerfile
# Dockerfile
FROM rocker/tidyverse:4.3.1

# Install system dependencies
RUN apt-get update && apt-get install -y \\
    libcurl4-openssl-dev \\
    libssl-dev \\
    libxml2-dev \\
    libpng-dev \\
    libudunits2-dev \\
    libgdal-dev

# Install Bioconductor
RUN R -e "install.packages('BiocManager')"

# Install packages
COPY install_dependencies.R /tmp/
RUN Rscript /tmp/install_dependencies.R

# Set working directory
WORKDIR /workspace

CMD ["R"]
```

Build and run:
```bash
docker build -t rnaseq-analysis .
docker run -v $(pwd):/workspace -it rnaseq-analysis
```

## Troubleshooting

### Common Issues

**Issue:** Package installation fails with compilation errors

**Solution:** Install system dependencies (Ubuntu/Debian):
```bash
sudo apt-get install -y \\
  libcurl4-openssl-dev \\
  libssl-dev \\
  libxml2-dev \\
  libpng-dev \\
  libudunits2-dev
```

**Issue:** Bioconductor package version conflicts

**Solution:** Update all Bioconductor packages:
```r
BiocManager::install(ask = FALSE, update = TRUE)
```

**Issue:** `org.Ss.eg.db` not available for current R version

**Solution:** Install from source or use older Bioconductor release:
```r
BiocManager::install("org.Ss.eg.db", force = TRUE)
```

### Getting Help

- **Bioconductor Support:** https://support.bioconductor.org/
- **Stack Overflow:** Tag questions with `[r]` and `[deseq2]`
- **Package Documentation:** Use `?function_name` in R console

## Version Information

To document the exact versions used in your analysis:

```r
# Save session info
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

# Or in R Markdown
sessionInfo()
```

This information should be included in all analysis reports.

## Updates

Packages are regularly updated. To update all packages:

```r
# CRAN packages
update.packages(ask = FALSE)

# Bioconductor packages
BiocManager::install(ask = FALSE, update = TRUE)
```

**Note:** Major updates may introduce breaking changes. Test thoroughly after updating.

---

**Last Updated:** December 2025
