# Data Processing Utilities for RNAseq Analysis
# This file contains shared functions for data import, processing, and merging

#' Process a single dataset by reading, filtering, and aggregating
#'
#' @param filepath Path to the count data file (relative to data directory)
#' @param remove_pattern Optional regex pattern to remove rows (e.g., "^ENSSSCG" for pig genes)
#' @return Data frame with processed gene counts
process_dataset <- function(filepath, remove_pattern = NULL) {

  # Read in the data
  df <- read.table(filepath, header = TRUE, sep = "\t")

  # Optionally remove rows whose GeneID matches the remove_pattern
  if (!is.null(remove_pattern)) {
    df <- subset(df, !grepl(remove_pattern, df$GeneID))
  }

  # Remove underscore and everything following it from GeneID
  # This standardizes gene identifiers
  df$GeneID <- sub("_.*", "", df$GeneID)

  # Aggregate duplicates by summing numeric columns
  df <- aggregate(. ~ GeneID, data = df, FUN = sum)

  return(df)
}


#' Merge multiple datasets on common gene identifiers
#'
#' @param dataset_list List of data frames to merge
#' @param join_by Column name to join on (default: "GeneID")
#' @param join_type Type of join: "inner" or "full" (default: "inner")
#' @return Merged data frame
merge_datasets <- function(dataset_list, join_by = "GeneID", join_type = "inner") {

  if (join_type == "inner") {
    merged_data <- purrr::reduce(dataset_list, dplyr::inner_join, by = join_by)
  } else if (join_type == "full") {
    merged_data <- purrr::reduce(dataset_list, dplyr::full_join, by = join_by)
  } else {
    stop("join_type must be 'inner' or 'full'")
  }

  return(merged_data)
}


#' Prepare count matrix for DESeq2 analysis
#'
#' @param merged_data Data frame with gene counts
#' @return Count matrix with genes as rownames
prepare_count_matrix <- function(merged_data) {

  # Set gene IDs as rownames
  rownames(merged_data) <- merged_data[, 1]
  merged_data <- merged_data[, -1]

  # Convert to integer counts (required for DESeq2)
  merged_data[] <- as.data.frame(lapply(merged_data, function(x) as.integer(round(x))))

  return(merged_data)
}


#' Load sample metadata and create DESeqDataSet
#'
#' @param count_matrix Count matrix with genes as rownames
#' @param metadata_file Path to sample metadata file
#' @param design_formula Design formula for DESeq2 (as string, e.g., "~CellType")
#' @return DESeqDataSet object
create_deseq_dataset <- function(count_matrix, metadata_file, design_formula = "~CellType") {

  # Load sample metadata
  colData <- read.table(metadata_file, header = TRUE, sep = "\t")

  # Convert design formula string to formula object
  design <- stats::as.formula(design_formula)

  # Create DESeqDataSet
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = colData,
    design = design
  )

  return(dds)
}


#' Set factor levels for DESeq2 comparison
#'
#' @param dds DESeqDataSet object
#' @param factor_name Name of the factor column (e.g., "CellType")
#' @param level_order Vector of factor levels in desired order (baseline first)
#' @return DESeqDataSet with reordered factor levels
set_factor_levels <- function(dds, factor_name, level_order) {

  SummarizedExperiment::colData(dds)[[factor_name]] <- factor(
    SummarizedExperiment::colData(dds)[[factor_name]],
    levels = level_order
  )

  return(dds)
}
