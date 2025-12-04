# ============================================
# HEATMAP GENERATOR
# ============================================
# Create heatmaps with flexible gene and sample selection
#
# USAGE (RStudio):
#   1. Open this script in RStudio
#   2. Edit the parameters below
#   3. Click "Source" or press Ctrl/Cmd+Shift+S
#
# USAGE (Terminal):
#   Rscript run_heatmap.R
#
# CONFIGURATION:
#   Edit the parameters below before running

# ============================================
# PARAMETERS - EDIT THESE
# ============================================

# Analysis name (used for output filename)
ANALYSIS_NAME <- "WNT_heatmap"

# Output directory
OUTPUT_DIR <- "figures"

# Input files
COUNTS_FILE <- "figures/full_test_pPSC_vs_hPSC/normalized_counts.xlsx"
METADATA_FILE <- "config/sample_metadata_full.txt"

# Gene selection
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

# Sample selection
# Leave NULL to use all samples
SAMPLE_CRITERIA <- NULL  # Or: list(CellType = "pPSC")
# Examples:
# list(CellType = "hPSC")  # Only human samples
# list(Dataset = c("hESC_ours", "hESC_Selmi"))  # Specific datasets

# Heatmap parameters
CLUSTER_ROWS <- TRUE
CLUSTER_COLS <- TRUE
COLOR_SCHEME <- "viridis"  # "viridis", "blue_red", or "blue_yellow"
SCALE_METHOD <- "0to1"  # "0to1", "row", "column", or "none"
ANNOTATION_COLS <- c("CellType", "Dataset")  # Metadata columns to display

# Plot dimensions
PLOT_WIDTH <- 12
PLOT_HEIGHT <- 10

# ============================================
# SETUP
# ============================================

cat("=============================================\n")
cat("Heatmap Generator\n")
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
  library("pheatmap")
  library("viridis")
  library("RColorBrewer")
})

# Source utility functions
source("scripts/utils/heatmap_helpers.R")
source("scripts/utils/gene_set_operations.R")

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================
# 1. LOAD DATA
# ============================================
cat("\n=== Loading Data ===\n")

# Load normalized counts
counts <- load_normalized_counts(COUNTS_FILE)

# Load metadata
metadata <- read.table(METADATA_FILE, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
cat("Loaded metadata for", nrow(metadata), "samples\n")

# ============================================
# 2. SELECT GENES
# ============================================
cat("\n=== Selecting Genes ===\n")

if (USE_ALL_GENES) {
  gene_list <- NULL
  cat("Using all genes\n")

} else if (USE_PATHWAY_GENES && file.exists(PATHWAY_FILE)) {
  gene_list <- load_pathway_genes(PATHWAY_FILE)
  cat("Using pathway genes:", length(gene_list), "\n")

} else if (USE_GENE_LIST_FILE && file.exists(GENE_LIST_FILE)) {
  gene_list <- read.table(GENE_LIST_FILE, stringsAsFactors = FALSE)[[1]]
  cat("Loaded gene list from file:", length(gene_list), "genes\n")

} else if (USE_CUSTOM_GENES) {
  gene_list <- CUSTOM_GENES
  cat("Using custom gene list:", length(gene_list), "genes\n")

} else {
  stop("No valid gene selection method specified!")
}

# ============================================
# 3. SELECT SAMPLES
# ============================================
cat("\n=== Selecting Samples ===\n")

if (!is.null(SAMPLE_CRITERIA)) {
  sample_list <- get_samples_by_criteria(metadata, SAMPLE_CRITERIA)
} else {
  sample_list <- NULL
  cat("Using all samples\n")
}

# ============================================
# 4. GENERATE HEATMAP
# ============================================
cat("\n=== Generating Heatmap ===\n")

output_file <- file.path(OUTPUT_DIR, paste0(ANALYSIS_NAME, ".pdf"))

generate_flexible_heatmap(
  count_data = counts,
  gene_list = gene_list,
  sample_list = sample_list,
  metadata = metadata,
  output_file = output_file,
  title = ANALYSIS_NAME,
  cluster_rows = CLUSTER_ROWS,
  cluster_cols = CLUSTER_COLS,
  color_scheme = COLOR_SCHEME,
  scale_method = SCALE_METHOD,
  annotation_col = ANNOTATION_COLS,
  width = PLOT_WIDTH,
  height = PLOT_HEIGHT
)

# ============================================
# SUMMARY
# ============================================
cat("\n=============================================\n")
cat("✓ Complete!\n")
cat("=============================================\n\n")

cat("Output file:", output_file, "\n")
