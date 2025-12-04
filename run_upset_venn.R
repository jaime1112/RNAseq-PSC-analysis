# ============================================
# UPSET PLOT AND VENN DIAGRAM GENERATOR
# ============================================
# Compare gene lists across multiple comparisons
#
# USAGE (RStudio):
#   1. Open this script in RStudio
#   2. Edit the parameters below
#   3. Click "Source" or press Ctrl/Cmd+Shift+S
#
# USAGE (Terminal):
#   Rscript run_upset_venn.R
#
# CONFIGURATION:
#   Edit the parameters below before running

# ============================================
# PARAMETERS - EDIT THESE
# ============================================

# Analysis name (used for output directory)
ANALYSIS_NAME <- "WNT_multi_comparison"

# Output directory
OUTPUT_DIR <- "figures"

# Method 1: Load from multiple comparison results
USE_MULTIPLE_COMPARISONS <- FALSE  # Set to TRUE to use Method 1

# Paths to DE results files (for Method 1)
COMPARISON_RESULTS <- c(
  "Comparison1" = "figures/comp1/allGenes.xlsx",
  "Comparison2" = "figures/comp2/allGenes.xlsx",
  "Comparison3" = "figures/comp3/allGenes.xlsx"
)

# Method 2: Load from single result file (create categories)
USE_SINGLE_RESULT <- TRUE  # Set to TRUE to use Method 2

# Path to single DE results file (for Method 2)
SINGLE_RESULT_FILE <- "figures/full_test_pPSC_vs_hPSC/allGenes.xlsx"

# Gene filtering
DIRECTION <- "both"  # "up", "down", or "both"
LFC_THRESHOLD <- 1
PADJ_THRESHOLD <- 0.05

# Pathway filtering (optional)
FILTER_PATHWAY <- TRUE
PATHWAY_FILE <- "config/WntSignalingPathway.txt"

# ============================================
# SETUP
# ============================================

cat("=============================================\n")
cat("UpSet Plot and Venn Diagram Generator\n")
cat("=============================================\n\n")

# Set working directory to script location
# This works in both RStudio and command line
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  # Running in RStudio
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
} else {
  # Running from command line
  setwd(dirname(sys.frame(1)$ofile))
}

# Load libraries
cat("Loading libraries...\n")
suppressPackageStartupMessages({
  library("tidyverse")
  library("readxl")
  library("writexl")
  library("VennDiagram")
  library("UpSetR")
  library("grid")
})

# Source utility functions
source("scripts/utils/gene_set_operations.R")

# Create output directory
output_dir <- file.path(OUTPUT_DIR, ANALYSIS_NAME)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================
# 1. LOAD PATHWAY GENES (OPTIONAL)
# ============================================
if (FILTER_PATHWAY && file.exists(PATHWAY_FILE)) {
  cat("\n=== Loading Pathway Genes ===\n")
  pathway_genes <- load_pathway_genes(PATHWAY_FILE)
} else {
  pathway_genes <- NULL
}

# ============================================
# 2. LOAD GENE LISTS
# ============================================
cat("\n=== Loading Gene Lists ===\n")

if (USE_MULTIPLE_COMPARISONS) {
  # Method 1: Load from multiple comparison results
  gene_lists <- load_de_gene_lists(
    results_files = COMPARISON_RESULTS,
    lfc_threshold = LFC_THRESHOLD,
    padj_threshold = PADJ_THRESHOLD,
    direction = DIRECTION,
    gene_filter = pathway_genes
  )

} else if (USE_SINGLE_RESULT) {
  # Method 2: Create categories from single result
  res <- read_excel(SINGLE_RESULT_FILE)
  res <- as.data.frame(res)

  gene_lists <- list()

  if (DIRECTION %in% c("up", "both")) {
    gene_lists$High_Up <- res$Row.names[
      !is.na(res$log2FoldChange) & !is.na(res$padj) &
      res$log2FoldChange > 2 & res$padj < PADJ_THRESHOLD
    ]
    gene_lists$Moderate_Up <- res$Row.names[
      !is.na(res$log2FoldChange) & !is.na(res$padj) &
      res$log2FoldChange > 1 & res$log2FoldChange <= 2 & res$padj < PADJ_THRESHOLD
    ]
  }

  if (DIRECTION %in% c("down", "both")) {
    gene_lists$High_Down <- res$Row.names[
      !is.na(res$log2FoldChange) & !is.na(res$padj) &
      res$log2FoldChange < -2 & res$padj < PADJ_THRESHOLD
    ]
    gene_lists$Moderate_Down <- res$Row.names[
      !is.na(res$log2FoldChange) & !is.na(res$padj) &
      res$log2FoldChange < -1 & res$log2FoldChange >= -2 & res$padj < PADJ_THRESHOLD
    ]
  }

  # Filter for pathway if specified
  if (!is.null(pathway_genes)) {
    gene_lists <- lapply(gene_lists, function(x) intersect(x, pathway_genes))
  }

  # Print summary
  for (name in names(gene_lists)) {
    cat(sprintf("  %-20s: %d genes\n", name, length(gene_lists[[name]])))
  }
}

# ============================================
# 3. GENERATE UPSET PLOT
# ============================================
cat("\n=== Generating UpSet Plot ===\n")

generate_upset_plot(
  gene_lists = gene_lists,
  output_file = file.path(output_dir, "UpSet_gene_overlaps.pdf"),
  title = "Gene Set Overlaps",
  n_intersects = 20,
  width = 12,
  height = 7
)

# ============================================
# 4. GENERATE VENN DIAGRAM (if 2-4 lists)
# ============================================
if (length(gene_lists) >= 2 && length(gene_lists) <= 4) {
  cat("\n=== Generating Venn Diagram ===\n")

  generate_venn_diagram(
    gene_lists = gene_lists,
    output_file = file.path(output_dir, "Venn_gene_overlaps.pdf"),
    title = "Gene Set Overlaps",
    width = 10,
    height = 10
  )
} else {
  cat("\nSkipping Venn diagram (need 2-4 gene lists)\n")
}

# ============================================
# 5. FIND AND EXPORT OVERLAPS
# ============================================
cat("\n=== Finding Overlaps ===\n")

# Universal genes (in all lists)
universal_genes <- find_overlap_genes(gene_lists, pattern = "all")

# Export gene lists with overlaps
export_gene_lists(
  gene_lists = gene_lists,
  output_file = file.path(output_dir, "gene_lists_with_overlaps.xlsx"),
  include_overlaps = TRUE
)

# Create summary table
summary_table <- summarize_overlaps(gene_lists)
write_xlsx(summary_table, file.path(output_dir, "gene_list_summary.xlsx"))
print(summary_table)

# ============================================
# SUMMARY
# ============================================
cat("\n=============================================\n")
cat("✓ Complete!\n")
cat("=============================================\n\n")

cat("Output directory:", output_dir, "\n\n")
cat("Files created:\n")
cat("  - UpSet_gene_overlaps.pdf\n")
if (length(gene_lists) >= 2 && length(gene_lists) <= 4) {
  cat("  - Venn_gene_overlaps.pdf\n")
}
cat("  - gene_lists_with_overlaps.xlsx\n")
cat("  - gene_list_summary.xlsx\n\n")

if (length(universal_genes) > 0) {
  cat("Universal genes (in ALL lists):", length(universal_genes), "\n")
  cat("Genes:", paste(head(universal_genes, 20), collapse = ", "), "\n")
  if (length(universal_genes) > 20) cat("... and", length(universal_genes) - 20, "more\n")
}
