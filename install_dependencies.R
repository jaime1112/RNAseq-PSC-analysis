# Installation Script for RNAseq Analysis Pipeline
# Run this script to install all required R packages

cat("===========================================\n")
cat("RNAseq Analysis Pipeline - Dependency Installation\n")
cat("===========================================\n\n")

# Install CRAN packages
cat("Installing CRAN packages...\n")
cran_packages <- c(
  # Data manipulation and file I/O
  "tidyverse",      # Meta-package: dplyr, ggplot2, tidyr, readr, purrr, etc.
  "readxl",         # Read Excel files
  "writexl",        # Export to Excel format

  # Visualization - Core
  "ggplot2",        # Grammar of graphics plotting (also in tidyverse)
  "ggrepel",        # Non-overlapping text labels
  "ggpubr",         # Publication-ready plots

  # Visualization - Heatmaps and color schemes
  "pheatmap",       # Pretty heatmaps
  "RColorBrewer",   # Color palettes
  "viridis",        # Perceptually uniform color scales

  # Set visualization
  "VennDiagram",    # Venn diagrams
  "UpSetR",         # UpSet plots for set intersections

  # RStudio integration
  "rstudioapi"      # RStudio API for path detection
)

install.packages(cran_packages, dependencies = TRUE)

# Install Bioconductor
cat("\nInstalling Bioconductor packages...\n")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

bioc_packages <- c(
  # Differential expression
  "DESeq2",           # RNA-seq differential expression analysis

  # Annotation databases
  "org.Hs.eg.db",     # Human genome annotation
  "org.Ss.eg.db",     # Pig genome annotation
  "AnnotationDbi",    # Annotation database interface

  # Enrichment analysis
  "clusterProfiler",  # GO/KEGG enrichment analysis
  "enrichplot"        # Enrichment visualization
)

BiocManager::install(bioc_packages, ask = FALSE)

# Verify installation
cat("\n===========================================\n")
cat("Verifying Installation\n")
cat("===========================================\n\n")

all_packages <- c(cran_packages, bioc_packages)
installed <- sapply(all_packages, function(pkg) {
  result <- require(pkg, character.only = TRUE, quietly = TRUE)
  status <- ifelse(result, "✓ INSTALLED", "✗ FAILED")
  cat(sprintf("%-25s %s\n", pkg, status))
  return(result)
})

cat("\n===========================================\n")
if (all(installed)) {
  cat("✓ All packages installed successfully!\n")
  cat("You can now run the analysis scripts.\n")
} else {
  failed <- names(installed)[!installed]
  cat("✗ Some packages failed to install:\n")
  cat(paste("  -", failed, collapse = "\n"), "\n")
  cat("\nPlease check the error messages above and try installing failed packages individually.\n")
}
cat("===========================================\n")
