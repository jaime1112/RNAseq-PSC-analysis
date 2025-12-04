# ============================================
# DIFFERENTIAL EXPRESSION ANALYSIS PIPELINE
# ============================================
# Complete analysis including DE, enrichment, and visualizations
#
# USAGE (RStudio):
#   1. Open this script in RStudio
#   2. Edit the parameters below
#   3. Click "Source" or press Ctrl/Cmd+Shift+S
#
# USAGE (Terminal):
#   Rscript run_differential_expression.R
#
# CONFIGURATION:
#   Edit the parameters below before running

# ============================================
# PARAMETERS - EDIT THESE
# ============================================

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

# Output directory
OUTPUT_DIR <- "figures"

# ============================================
# SETUP
# ============================================

cat("=============================================\n")
cat("Differential Expression Analysis Pipeline\n")
cat("=============================================\n\n")

cat("Analysis:", ANALYSIS_NAME, "\n")
cat("Comparison:", COMPARISON_LEVEL, "vs", BASELINE_LEVEL, "\n")
cat("Organism:", ORGANISM, "\n\n")

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
  library("DESeq2")
  library("writexl")
  library("ggplot2")
  library("ggrepel")
  library("ggpubr")
  library("org.Hs.eg.db")
  library("org.Ss.eg.db")
  if (PERFORM_ENRICHMENT) {
    library("clusterProfiler")
    library("enrichplot")
  }
})

# Source utility functions
source("scripts/utils/data_processing.R")
source("scripts/utils/deseq_helpers.R")
source("scripts/utils/plotting_themes.R")
if (PERFORM_ENRICHMENT) {
  source("scripts/utils/enrichment_analysis.R")
}

# Create output directory
analysis_dir <- file.path(OUTPUT_DIR, ANALYSIS_NAME)
dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================
# 1. LOAD AND PROCESS DATA
# ============================================
cat("\n=== Step 1: Loading Data ===\n")

data_config <- read.csv(DATA_CONFIG, stringsAsFactors = FALSE)
cat("Loading", nrow(data_config), "datasets\n")

dataset_list <- list()
for (i in 1:nrow(data_config)) {
  dataset_name <- data_config$dataset_name[i]
  filepath <- data_config$filepath[i]
  remove_pat <- data_config$remove_pattern[i]
  if (is.na(remove_pat) || remove_pat == "") remove_pat <- NULL

  cat("  -", dataset_name, "\n")
  dataset_list[[dataset_name]] <- process_dataset(filepath, remove_pattern = remove_pat)
}

merged_data <- merge_datasets(dataset_list)
count_matrix <- prepare_count_matrix(merged_data)
cat("Final matrix:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")

# ============================================
# 2. CREATE DESEQ DATASET
# ============================================
cat("\n=== Step 2: Creating DESeq Dataset ===\n")

dds <- create_deseq_dataset(count_matrix, METADATA_FILE, design_formula = DESIGN_FORMULA)
dds <- set_factor_levels(dds, FACTOR_NAME, c(BASELINE_LEVEL, COMPARISON_LEVEL))

# ============================================
# 3. RUN DESEQ2 ANALYSIS
# ============================================
cat("\n=== Step 3: Running DESeq2 Analysis ===\n")

dds <- run_deseq_analysis(dds)
rld <- variance_stabilize(dds)

# ============================================
# 4. GENERATE PCA PLOTS
# ============================================
cat("\n=== Step 4: Generating PCA Plots ===\n")

generate_pca_plot(rld, file.path(analysis_dir, "PCA_by_CellType.pdf"), intgroup = "CellType")
generate_pca_plot(rld, file.path(analysis_dir, "PCA_by_Dataset.pdf"), intgroup = "Dataset")

# ============================================
# 5. EXTRACT AND ANNOTATE RESULTS
# ============================================
cat("\n=== Step 5: Extracting Results ===\n")

contrast_name <- paste0(FACTOR_NAME, "_", COMPARISON_LEVEL, "_vs_", BASELINE_LEVEL)
res <- results(dds, contrast = c(FACTOR_NAME, COMPARISON_LEVEL, BASELINE_LEVEL), alpha = 0.05)
res <- DESeq2::lfcShrink(dds, coef = contrast_name, type = "apeglm")

gene_symbols <- rownames(count_matrix)
res <- annotate_results(res, gene_symbols, organism = ORGANISM)
res <- merge_counts_with_results(res, dds)
res_df <- as.data.frame(res)

cat("Significant genes (padj < 0.05):", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")

# Export results
write_xlsx(res_df, file.path(analysis_dir, "allGenes.xlsx"))
write_xlsx(filter_significant_genes(res_df), file.path(analysis_dir, "sigGenes.xlsx"))

nm_counts_df <- as.data.frame(counts(dds, normalized = TRUE))
nm_counts_df <- cbind("GeneID" = rownames(nm_counts_df), nm_counts_df)
rownames(nm_counts_df) <- NULL
write_xlsx(nm_counts_df, file.path(analysis_dir, "normalized_counts.xlsx"))

# ============================================
# 6. GENERATE MA PLOTS
# ============================================
cat("\n=== Step 6: Generating MA Plots ===\n")

gene_sets <- get_gene_sets()

generate_ma_plot(res_df, file.path(analysis_dir, "MAplot.pdf"),
                title = paste(COMPARISON_LEVEL, "vs", BASELINE_LEVEL), top = 20)

generate_ma_plot(res_df, file.path(analysis_dir, "MAplot_WNTligands.pdf"),
                title = "WNT Ligands", highlight_genes = gene_sets$wnt_ligands)

generate_ma_plot(res_df, file.path(analysis_dir, "MAplot_WNTreceptors.pdf"),
                title = "WNT Receptors", highlight_genes = gene_sets$wnt_receptors)

generate_ma_plot(res_df, file.path(analysis_dir, "MAplot_Pluripotency.pdf"),
                title = "Pluripotency Genes", highlight_genes = gene_sets$pluripotency)

generate_ma_plot_labeled(res_df, file.path(analysis_dir, "MAplot_WNTlinked.pdf"),
                         title = "WNT Pathway Genes", highlight_genes = gene_sets$wnt_linked)

# ============================================
# 7. ENRICHMENT ANALYSIS (OPTIONAL)
# ============================================
if (PERFORM_ENRICHMENT) {
  cat("\n=== Step 7: Enrichment Analysis ===\n")

  gene_lists <- get_regulated_genes(res_df, lfc_cutoff = 1, padj_cutoff = 0.05)
  cat("Upregulated:", length(gene_lists$Upregulated), "\n")
  cat("Downregulated:", length(gene_lists$Downregulated), "\n")

  enrichment_results <- run_enrichment_analysis(
    gene_lists = gene_lists,
    universe = res_df$entrezid,
    organism = ORGANISM,
    ontology = "BP"
  )

  write_xlsx(enrichment_results$GO@compareClusterResult,
            file.path(analysis_dir, "GO_combined.xlsx"))
  write_xlsx(enrichment_results$KEGG@compareClusterResult,
            file.path(analysis_dir, "KEGG_combined.xlsx"))

  generate_enrichment_dotplot(enrichment_results$GO,
                              file.path(analysis_dir, "GO_dotplot.pdf"),
                              title = "GO Enrichment")
  generate_enrichment_dotplot(enrichment_results$KEGG,
                              file.path(analysis_dir, "KEGG_dotplot.pdf"),
                              title = "KEGG Enrichment")
}

# ============================================
# 8. PATHWAY ANALYSIS (OPTIONAL)
# ============================================
if (RUN_PATHWAY_ANALYSIS) {
  cat("\n=== Step 8: Pathway Analysis ===\n")

  wnt_file <- "config/WntSignalingPathway.txt"
  if (file.exists(wnt_file)) {
    wnt_results <- filter_pathway_genes(wnt_file, res_df, gene_col = "GeneID")
    write_xlsx(wnt_results, file.path(analysis_dir, "WntGenes.xlsx"))
    cat("WNT pathway genes exported\n")
  }
}

# ============================================
# SUMMARY
# ============================================
cat("\n=============================================\n")
cat("✓ Analysis Complete!\n")
cat("=============================================\n\n")

cat("Output directory:", analysis_dir, "\n\n")
cat("Files created:\n")
cat("  - allGenes.xlsx (all genes with statistics)\n")
cat("  - sigGenes.xlsx (significant genes only)\n")
cat("  - normalized_counts.xlsx (DESeq2 normalized counts)\n")
cat("  - PCA_by_CellType.pdf\n")
cat("  - PCA_by_Dataset.pdf\n")
cat("  - MAplot.pdf (main)\n")
cat("  - MAplot_*.pdf (pathway-specific, 4 plots)\n")
if (PERFORM_ENRICHMENT) {
  cat("  - GO_combined.xlsx\n")
  cat("  - KEGG_combined.xlsx\n")
  cat("  - GO_dotplot.pdf\n")
  cat("  - KEGG_dotplot.pdf\n")
}
if (RUN_PATHWAY_ANALYSIS) {
  cat("  - WntGenes.xlsx\n")
}
cat("\n")
