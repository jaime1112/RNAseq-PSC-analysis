# Installation Script for RNAseq Analysis Pipeline
# Run this script to install all required R packages

cat("===========================================\n")
cat("RNAseq Analysis Pipeline - Dependency Installation\n")
cat("===========================================\n\n")

# Install CRAN packages
cat("Installing CRAN packages...\n")
cran_packages <- c(
  # Data manipulation
  "tidyverse",
  "writexl",

  # Visualization
  "ggplot2",
  "ggrepel",
  "gplots",
  "vegan",
  "pheatmap",
  "RColorBrewer",
  "FactoMineR",
  "ggpubr",
  "ggsci",
  "gghighlight",
  "VennDiagram",
  "UpSetR",
  "viridis",
  "scico",

  # R Markdown
  "knitr",
  "rmarkdown"
)

install.packages(cran_packages, dependencies = TRUE)

# Install Bioconductor
cat("\nInstalling Bioconductor...\n")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

bioc_packages <- c(
  # Differential expression
  "DESeq2",
  "tximport",
  "apeglm",

  # Annotation
  "org.Hs.eg.db",
  "org.Ss.eg.db",
  "AnnotationDbi",

  # Enrichment analysis
  "clusterProfiler",
  "enrichplot",
  "DOSE"
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
