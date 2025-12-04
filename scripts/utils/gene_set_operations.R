# Gene Set Operations and Comparison Functions
# Functions for comparing gene lists across multiple analyses

#' Load differentially expressed gene lists from analysis results
#'
#' @param results_files Named vector of paths to DE results files (Excel or RDS)
#' @param lfc_threshold Log2 fold change threshold (default: 1)
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @param direction "up", "down", or "both" (default: "both")
#' @param gene_filter Optional vector of genes to filter for (e.g., WNT pathway)
#' @return Named list of gene vectors
load_de_gene_lists <- function(results_files, lfc_threshold = 1, padj_threshold = 0.05,
                                direction = "both", gene_filter = NULL) {

  gene_lists <- list()

  for (comparison_name in names(results_files)) {
    file_path <- results_files[comparison_name]

    # Load results
    if (grepl("\\.xlsx$", file_path)) {
      res <- readxl::read_excel(file_path)
      res <- as.data.frame(res)
    } else if (grepl("\\.rds$", file_path)) {
      res <- readRDS(file_path)
    } else {
      stop(paste("Unsupported file format:", file_path))
    }

    # Filter by significance and fold change
    res_sig <- res[!is.na(res$padj) & !is.na(res$log2FoldChange) &
                   res$padj < padj_threshold &
                   abs(res$log2FoldChange) > lfc_threshold, ]

    # Filter by direction
    if (direction == "up") {
      res_sig <- res_sig[res_sig$log2FoldChange > 0, ]
    } else if (direction == "down") {
      res_sig <- res_sig[res_sig$log2FoldChange < 0, ]
    }

    # Get gene symbols
    if ("Row.names" %in% colnames(res_sig)) {
      genes <- res_sig$Row.names
    } else if ("GeneSymbol" %in% colnames(res_sig)) {
      genes <- res_sig$GeneSymbol
    } else {
      genes <- rownames(res_sig)
    }

    # Apply gene filter if provided
    if (!is.null(gene_filter)) {
      genes <- intersect(genes, gene_filter)
    }

    gene_lists[[comparison_name]] <- genes
    message(paste(comparison_name, ":", length(genes), "genes"))
  }

  return(gene_lists)
}


#' Load pathway genes from file
#'
#' @param pathway_file Path to pathway gene file
#' @param gene_column Name of column containing gene symbols (default: "GeneID")
#' @return Vector of gene symbols
load_pathway_genes <- function(pathway_file, gene_column = "GeneID") {

  if (!file.exists(pathway_file)) {
    stop(paste("Pathway file not found:", pathway_file))
  }

  pathway_data <- read.table(pathway_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  if (!gene_column %in% colnames(pathway_data)) {
    stop(paste("Column", gene_column, "not found in pathway file"))
  }

  genes <- pathway_data[[gene_column]]
  message(paste("Loaded", length(genes), "genes from pathway file"))

  return(genes)
}


#' Find genes in specific overlap patterns
#'
#' @param gene_lists Named list of gene vectors
#' @param pattern "all" (in all lists), "any" (in any list), "unique_to" (in only one list),
#'                or custom pattern like c("Comp1", "Comp2") for specific overlap
#' @param comparison_name For "unique_to" pattern, name of the comparison
#' @return Vector of genes matching the pattern
find_overlap_genes <- function(gene_lists, pattern = "all", comparison_name = NULL) {

  # Check if pattern is a vector of comparison names
  if (is.character(pattern) && length(pattern) > 1) {
    # Custom pattern: genes in specified comparisons
    if (!all(pattern %in% names(gene_lists))) {
      stop("All pattern names must be in gene_lists")
    }
    overlap_genes <- Reduce(intersect, gene_lists[pattern])
    message(paste("Genes in", paste(pattern, collapse = " AND "), ":", length(overlap_genes)))

  } else if (length(pattern) == 1 && pattern == "all") {
    # Genes in all comparisons
    overlap_genes <- Reduce(intersect, gene_lists)
    message(paste("Genes in ALL comparisons:", length(overlap_genes)))

  } else if (length(pattern) == 1 && pattern == "any") {
    # Genes in any comparison
    overlap_genes <- Reduce(union, gene_lists)
    message(paste("Genes in ANY comparison:", length(overlap_genes)))

  } else if (length(pattern) == 1 && pattern == "unique_to") {
    # Genes unique to one comparison
    if (is.null(comparison_name) || !comparison_name %in% names(gene_lists)) {
      stop("comparison_name must be specified and valid for 'unique_to' pattern")
    }

    target_genes <- gene_lists[[comparison_name]]
    other_genes <- Reduce(union, gene_lists[names(gene_lists) != comparison_name])
    overlap_genes <- setdiff(target_genes, other_genes)
    message(paste("Genes unique to", comparison_name, ":", length(overlap_genes)))

  } else {
    stop("Invalid pattern. Use 'all', 'any', 'unique_to', or vector of comparison names")
  }

  return(overlap_genes)
}


#' Generate UpSet plot from gene lists
#'
#' @param gene_lists Named list of gene vectors
#' @param output_file Path to output PDF
#' @param title Plot title
#' @param n_intersects Number of intersections to show (default: 20)
#' @param width Plot width (default: 10)
#' @param height Plot height (default: 6)
generate_upset_plot <- function(gene_lists, output_file, title = "Gene Set Overlaps",
                                n_intersects = 20, width = 10, height = 6) {

  # Convert to binary matrix
  upset_data <- UpSetR::fromList(gene_lists)

  # Generate plot
  pdf(file = output_file, width = width, height = height)

  upset_plot <- UpSetR::upset(
    upset_data,
    sets = names(gene_lists),
    nintersects = n_intersects,
    order.by = "freq",
    mainbar.y.label = "Gene Intersection Size",
    sets.x.label = "Total Genes Per Set",
    text.scale = 1.5,
    main.bar.color = "#2166ac",
    sets.bar.color = "#b2182b"
  )

  print(upset_plot)
  dev.off()

  message(paste("UpSet plot saved to:", output_file))
  return(invisible(output_file))
}


#' Generate Venn diagram from gene lists (2-4 way)
#'
#' @param gene_lists Named list of gene vectors (2-4 lists)
#' @param output_file Path to output PDF
#' @param title Plot title
#' @param width Plot width (default: 8)
#' @param height Plot height (default: 8)
generate_venn_diagram <- function(gene_lists, output_file, title = NULL,
                                  width = 8, height = 8) {

  # Check number of lists
  n_lists <- length(gene_lists)
  if (n_lists < 2 || n_lists > 4) {
    stop("Venn diagrams support 2-4 gene lists. Use UpSet plots for more.")
  }

  # Define colors
  colors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")[1:n_lists]

  # Create Venn diagram
  venn_plot <- VennDiagram::venn.diagram(
    x = gene_lists,
    category.names = names(gene_lists),
    filename = NULL,
    output = TRUE,
    imagetype = "png",
    height = 480,
    width = 480,
    resolution = 300,
    compression = "lzw",
    lwd = 2,
    lty = "solid",
    fill = colors,
    cex = 1.2,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.dist = c(0.055, 0.055, 0.055, 0.055)[1:n_lists],
    margin = 0.1
  )

  # Save to PDF
  pdf(file = output_file, width = width, height = height)
  if (!is.null(title)) {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 2, ncol = 1, heights = grid::unit(c(1, 10), "null"))))
    grid::grid.text(title, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1), gp = grid::gpar(fontsize = 16, fontface = "bold"))
    grid::grid.draw(venn_plot)
  } else {
    grid::grid.draw(venn_plot)
  }
  dev.off()

  message(paste("Venn diagram saved to:", output_file))
  return(invisible(output_file))
}


#' Export gene lists to Excel with overlap information
#'
#' @param gene_lists Named list of gene vectors
#' @param output_file Path to output Excel file
#' @param include_overlaps Include sheets for overlaps (default: TRUE)
export_gene_lists <- function(gene_lists, output_file, include_overlaps = TRUE) {

  # Create list of data frames
  excel_list <- list()

  # Individual gene lists
  max_length <- max(sapply(gene_lists, length))
  for (name in names(gene_lists)) {
    genes <- gene_lists[[name]]
    # Pad with NA to equal length
    padded_genes <- c(genes, rep(NA, max_length - length(genes)))
    excel_list[[name]] <- data.frame(Gene = padded_genes, stringsAsFactors = FALSE)
  }

  # Add overlap sheets
  if (include_overlaps && length(gene_lists) >= 2) {
    # All overlaps
    all_overlap <- find_overlap_genes(gene_lists, pattern = "all")
    if (length(all_overlap) > 0) {
      excel_list[["All_Overlaps"]] <- data.frame(Gene = all_overlap, stringsAsFactors = FALSE)
    }

    # Pairwise overlaps (for 2-4 lists)
    if (length(gene_lists) <= 4) {
      list_names <- names(gene_lists)
      for (i in 1:(length(list_names) - 1)) {
        for (j in (i + 1):length(list_names)) {
          pair_name <- paste(list_names[i], list_names[j], sep = "_AND_")
          pair_overlap <- intersect(gene_lists[[i]], gene_lists[[j]])
          if (length(pair_overlap) > 0) {
            excel_list[[pair_name]] <- data.frame(Gene = pair_overlap, stringsAsFactors = FALSE)
          }
        }
      }
    }
  }

  # Write to Excel
  writexl::write_xlsx(excel_list, output_file)
  message(paste("Gene lists exported to:", output_file))
  message(paste("Total sheets:", length(excel_list)))

  return(invisible(output_file))
}


#' Create summary table of gene list overlaps
#'
#' @param gene_lists Named list of gene vectors
#' @return Data frame with overlap statistics
summarize_overlaps <- function(gene_lists) {

  n_lists <- length(gene_lists)
  list_names <- names(gene_lists)

  # Create summary data frame
  summary_df <- data.frame(
    Comparison = list_names,
    Total_Genes = sapply(gene_lists, length),
    stringsAsFactors = FALSE
  )

  # Add overlap columns
  if (n_lists >= 2) {
    # Common to all
    all_overlap <- find_overlap_genes(gene_lists, pattern = "all")
    summary_df$In_All <- sapply(gene_lists, function(x) sum(x %in% all_overlap))

    # Unique genes
    for (i in 1:n_lists) {
      comparison <- list_names[i]
      unique_genes <- find_overlap_genes(gene_lists, pattern = "unique_to", comparison_name = comparison)
      summary_df$Unique[i] <- length(unique_genes)
    }

    # Percentage in all
    summary_df$Percent_In_All <- round(100 * summary_df$In_All / summary_df$Total_Genes, 2)
  }

  return(summary_df)
}
