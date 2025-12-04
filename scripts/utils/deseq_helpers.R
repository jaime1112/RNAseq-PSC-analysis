# DESeq2 Analysis Helper Functions
# Functions for running differential expression analysis and annotation

#' Run DESeq2 analysis pipeline
#'
#' @param dds DESeqDataSet object
#' @return DESeqDataSet with results
run_deseq_analysis <- function(dds) {
  dds <- DESeq2::DESeq(dds)
  return(dds)
}


#' Perform variance stabilizing transformation for visualization
#'
#' @param dds DESeqDataSet object
#' @return Transformed data object
variance_stabilize <- function(dds) {
  rld <- DESeq2::vst(dds)
  return(rld)
}


#' Extract and shrink differential expression results
#'
#' @param dds DESeqDataSet object
#' @param contrast Vector specifying contrast (e.g., c("CellType", "pPSC", "hPSC"))
#' @param coef Coefficient name for shrinkage (e.g., "CellType_pPSC_vs_hPSC")
#' @param alpha Significance threshold (default: 0.05)
#' @param shrink_method Shrinkage method (default: "apeglm")
#' @return Results object with shrunken log2 fold changes
extract_results <- function(dds, contrast, coef, alpha = 0.05, shrink_method = "apeglm") {

  # Get raw results
  res <- DESeq2::results(dds, contrast = contrast, alpha = alpha)

  # Apply shrinkage to log2 fold changes
  res <- DESeq2::lfcShrink(dds, coef = coef, type = shrink_method)

  return(res)
}


#' Annotate results with gene information
#'
#' @param res DESeq2 results object
#' @param gene_symbols Gene symbols from rownames
#' @param organism Organism database ("Hs" for human, "Ss" for pig)
#' @return Annotated results data frame
annotate_results <- function(res, gene_symbols, organism = "human") {

  # Get organism-specific database
  # Normalize organism identifier
  org_lower <- tolower(organism)

  if (org_lower %in% c("hs", "human", "hsa", "homo sapiens")) {
    org_db <- org.Hs.eg.db::org.Hs.eg.db
    org_name <- "Human"
  } else if (org_lower %in% c("ss", "pig", "ssc", "sus scrofa", "porcine")) {
    org_db <- org.Ss.eg.db::org.Ss.eg.db
    org_name <- "Pig"
  } else {
    stop(paste("Unknown organism:", organism,
               "\nSupported: 'Hs'/'human'/'hsa' (human) or 'Ss'/'pig'/'ssc' (pig)"))
  }

  message(paste("Annotating results with", org_name, "gene information"))

  # Define columns to retrieve
  columns <- c("GENENAME", "ENTREZID", "ENSEMBL")

  # Retrieve annotations
  ann <- AnnotationDbi::select(
    org_db,
    keys = gene_symbols,
    keytype = "SYMBOL",
    columns = columns
  )

  # Remove duplicates (keep first occurrence)
  ann <- ann[!duplicated(ann$SYMBOL), ]

  # Set row names
  rownames(ann) <- ann$SYMBOL

  # Add annotations to results
  res@listData$genename <- ann$GENENAME
  res@listData$ensembl <- ann$ENSEMBL
  res@listData$entrezid <- ann$ENTREZID

  return(res)
}


#' Merge results with normalized counts
#'
#' @param res Annotated DESeq2 results
#' @param dds DESeqDataSet object
#' @return Data frame with results and normalized counts
merge_counts_with_results <- function(res, dds) {

  # Extract normalized counts
  normalized_counts <- DESeq2::counts(dds, normalized = TRUE)

  # Merge with results
  res <- merge(res, normalized_counts, by = "row.names")

  return(res)
}


#' Filter and export significant genes
#'
#' @param res_df Results data frame
#' @param padj_cutoff Adjusted p-value cutoff (default: 0.05)
#' @return Data frame with significant genes only
filter_significant_genes <- function(res_df, padj_cutoff = 0.05) {
  sig_genes <- res_df[which(res_df$padj < padj_cutoff), ]
  return(sig_genes)
}


#' Extract gene lists by regulation direction
#'
#' @param res_df Results data frame
#' @param lfc_cutoff Log2 fold change cutoff (default: 1)
#' @param padj_cutoff Adjusted p-value cutoff (default: 0.05)
#' @return List with upregulated and downregulated genes
get_regulated_genes <- function(res_df, lfc_cutoff = 1, padj_cutoff = 0.05) {

  # Upregulated genes
  up_genes <- res_df$entrezid[
    res_df$log2FoldChange > lfc_cutoff &
    res_df$padj < padj_cutoff &
    !is.na(res_df$padj)
  ]

  # Downregulated genes
  down_genes <- res_df$entrezid[
    res_df$log2FoldChange < -lfc_cutoff &
    res_df$padj < padj_cutoff &
    !is.na(res_df$padj)
  ]

  return(list(
    Upregulated = up_genes,
    Downregulated = down_genes
  ))
}


#' Filter pathway genes by differential expression
#'
#' @param pathway_file Path to file with pathway gene list
#' @param res_df Results data frame
#' @param gene_col Column name with gene symbols in pathway file (default: "GeneID")
#' @param padj_cutoff Adjusted p-value cutoff (default: 0.05)
#' @return Data frame with pathway genes sorted by regulation
filter_pathway_genes <- function(pathway_file, res_df, gene_col = "GeneID", padj_cutoff = 0.05) {

  # Read pathway genes
  pathway <- read.table(pathway_file, header = TRUE, sep = "\t")
  pathway_symbols <- pathway[[gene_col]]

  # Filter results by regulation
  up_symbols <- res_df$Row.names[res_df$log2FoldChange > 0 & res_df$padj < padj_cutoff]
  nonsig_symbols <- res_df$Row.names[res_df$padj > padj_cutoff]
  down_symbols <- res_df$Row.names[res_df$log2FoldChange < 0 & res_df$padj < padj_cutoff]

  # Find pathway genes in each category
  common_up <- intersect(pathway_symbols, up_symbols)
  common_nonsig <- intersect(pathway_symbols, nonsig_symbols)
  common_down <- intersect(pathway_symbols, down_symbols)

  # Pad lists to equal length
  max_length <- max(length(common_up), length(common_nonsig), length(common_down))
  common_up <- c(common_up, rep(NA, max_length - length(common_up)))
  common_nonsig <- c(common_nonsig, rep(NA, max_length - length(common_nonsig)))
  common_down <- c(common_down, rep(NA, max_length - length(common_down)))

  # Create combined data frame
  combined_df <- data.frame(
    Common_Symbols_Up = common_up,
    Common_Symbols_NonSig = common_nonsig,
    Common_Symbols_Down = common_down
  )

  return(combined_df)
}
