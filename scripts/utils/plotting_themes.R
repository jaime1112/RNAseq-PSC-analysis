# Plotting Functions and Themes
# Standardized visualization functions for RNAseq analysis

#' Generate PCA plot and save to PDF
#'
#' @param rld Variance-stabilized data object
#' @param output_file Path to output PDF file
#' @param intgroup Factor to color by (default: "CellType")
#' @param width Plot width in inches (default: 8)
#' @param height Plot height in inches (default: 6)
generate_pca_plot <- function(rld, output_file, intgroup = "CellType", width = 8, height = 6) {

  # Open PDF device
  pdf(file = output_file, width = width, height = height)

  # Generate and print PCA plot
  pca_plot <- DESeq2::plotPCA(rld, intgroup = intgroup)
  print(pca_plot)

  # Close device
  dev.off()

  message(paste("PCA plot saved to:", output_file))
}


#' Generate MA plot with highlighted genes
#'
#' @param res_df Results data frame
#' @param output_file Path to output PDF file
#' @param title Plot title
#' @param highlight_genes Vector of gene symbols to highlight (optional)
#' @param fdr FDR threshold for significance (default: 0.05)
#' @param fc Fold change threshold for significance (default: 4)
#' @param top Number of top genes to label (default: 20)
#' @param width Plot width (default: 8)
#' @param height Plot height (default: 8)
generate_ma_plot <- function(res_df, output_file, title, highlight_genes = NULL,
                             fdr = 0.05, fc = 4, top = 20, width = 8, height = 8) {

  pdf(file = output_file, width = width, height = height)

  if (is.null(highlight_genes)) {
    # Standard MA plot with top genes
    maplot <- ggpubr::ggmaplot(
      res_df,
      main = title,
      fdr = fdr,
      fc = fc,
      size = 2,
      genenames = as.vector(res_df$Row.names),
      legend = "top",
      top = top,
      font.label = c("bold", 12),
      label.rectangle = FALSE,
      alpha = 0.5,
      font.legend = "bold",
      font.main = "bold",
      ggtheme = ggpubr::theme_pubr()
    )
    print(maplot)

  } else {
    # MA plot with specific highlighted genes
    highlighted_data <- res_df[
      res_df$Row.names %in% highlight_genes &
      abs(res_df$log2FoldChange) > log2(fc) &
      res_df$padj < fdr,
    ]

    maplot <- ggpubr::ggmaplot(
      res_df,
      main = title,
      fdr = fdr,
      fc = fc,
      size = 2,
      genenames = as.vector(res_df$Row.names),
      legend = "top",
      top = 0,
      label.select = highlighted_data$Row.names,
      font.label = c("bold", 12),
      label.rectangle = FALSE,
      alpha = 0.5,
      font.legend = "bold",
      font.main = "bold",
      ggtheme = ggpubr::theme_pubr()
    )

    # Add highlighted points
    maplot_with_points <- maplot +
      ggplot2::geom_point(
        data = highlighted_data,
        ggplot2::aes(x = log2(highlighted_data$baseMean), y = log2FoldChange),
        shape = 21,
        colour = "black",
        fill = "red",
        size = 3,
        stroke = 1
      )

    print(maplot_with_points)
  }

  dev.off()

  message(paste("MA plot saved to:", output_file))
}


#' Generate MA plot with gene labels using ggrepel
#'
#' @param res_df Results data frame (must have complete cases)
#' @param output_file Path to output PDF file
#' @param title Plot title
#' @param highlight_genes Vector of gene symbols to highlight
#' @param fdr FDR threshold (default: 0.05)
#' @param fc Fold change threshold (default: 4)
#' @param width Plot width (default: 8)
#' @param height Plot height (default: 8)
generate_ma_plot_labeled <- function(res_df, output_file, title, highlight_genes,
                                     fdr = 0.05, fc = 4, width = 8, height = 8) {

  # Filter out missing values
  res_filtered <- res_df[complete.cases(res_df$log2FoldChange, res_df$baseMean, res_df$padj), ]

  # Get highlighted gene data
  highlighted_data <- res_filtered[
    res_filtered$Row.names %in% highlight_genes &
    abs(res_filtered$log2FoldChange) > log2(fc) &
    res_filtered$padj < fdr,
  ]

  pdf(file = output_file, width = width, height = height)

  maplot <- ggpubr::ggmaplot(
    res_filtered,
    main = title,
    fdr = fdr,
    fc = fc,
    size = 2,
    genenames = as.vector(res_filtered$Row.names),
    legend = "top",
    top = 0,
    font.label = c("bold", 12),
    label.rectangle = FALSE,
    alpha = 0.5,
    font.legend = "bold",
    font.main = "bold",
    ggtheme = ggpubr::theme_pubr()
  )

  labeled_plot <- maplot +
    ggplot2::geom_point(
      data = highlighted_data,
      ggplot2::aes(x = log2(baseMean), y = log2FoldChange),
      shape = 21,
      colour = "black",
      fill = "red",
      size = 3,
      stroke = 1
    ) +
    ggrepel::geom_text_repel(
      data = highlighted_data,
      ggplot2::aes(x = log2(baseMean), y = log2FoldChange, label = Row.names),
      size = 4,
      fontface = "bold",
      box.padding = 1,
      point.padding = 1,
      segment.color = "black",
      segment.size = 1,
      max.overlaps = Inf
    )

  print(labeled_plot)
  dev.off()

  message(paste("Labeled MA plot saved to:", output_file))
}


#' Generate sample correlation heatmap and save to PDF
#'
#' @param rld Variance-stabilized data object (from vst or rlog)
#' @param output_file Path to output PDF file
#' @param method Correlation method (default: "pearson")
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
generate_correlation_heatmap <- function(rld, output_file,
                                         method = "pearson",
                                         width = NULL, height = NULL) {

  # Compute sample-sample correlation matrix
  sample_cor <- cor(SummarizedExperiment::assay(rld), method = method)
  n_samples <- ncol(sample_cor)

  # Scale dimensions and font size to sample count
  cellsize <- 12
  if (is.null(width))  width  <- max(10, n_samples * cellsize / 72 + 4)
  if (is.null(height)) height <- max(8,  n_samples * cellsize / 72 + 4)
  fontsize_labels <- max(5, min(8, 200 / n_samples))

  # Prepare metadata annotations — use all available metadata columns
  metadata <- as.data.frame(SummarizedExperiment::colData(rld))
  # Exclude non-informative columns (e.g., sizeFactor added by DESeq2)
  exclude_cols <- c("sizeFactor", "replaceable")
  annotation_cols <- setdiff(colnames(metadata), exclude_cols)
  # Order by number of unique values (fewest on top)
  annotation_cols <- annotation_cols[order(sapply(annotation_cols, function(col) length(unique(metadata[[col]]))), decreasing = TRUE)]

  if (length(annotation_cols) > 0) {
    annotation_df <- metadata[, annotation_cols, drop = FALSE]

    # Create annotation colors
    annotation_colors <- list()
    for (col in annotation_cols) {
      unique_vals <- unique(annotation_df[[col]])
      n_vals <- length(unique_vals)
      if (n_vals <= 12) {
        palette <- RColorBrewer::brewer.pal(min(max(n_vals, 3), 8), "Set2")
        if (n_vals > 8) {
          palette <- c(palette, RColorBrewer::brewer.pal(n_vals - 8, "Set3"))
        }
        annotation_colors[[col]] <- setNames(palette[1:n_vals], unique_vals)
      }
    }
  } else {
    annotation_df <- NULL
    annotation_colors <- NULL
  }

  # White-to-blue color scale (darker blue = higher correlation)
  colors <- grDevices::colorRampPalette(c("white", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C"))(100)

  pdf(file = output_file, width = width, height = height)

  pheatmap::pheatmap(
    mat = sample_cor,
    color = colors,
    clustering_distance_rows = as.dist(1 - sample_cor),
    clustering_distance_cols = as.dist(1 - sample_cor),
    clustering_method = "ward.D2",
    show_rownames = TRUE,
    show_colnames = TRUE,
    cellwidth = cellsize,
    cellheight = cellsize,
    fontsize_row = fontsize_labels,
    fontsize_col = fontsize_labels,
    border_color = NA,
    main = paste("Sample Correlation (", method, ")", sep = ""),
    annotation_col = annotation_df,
    annotation_row = annotation_df,
    annotation_colors = annotation_colors,
    treeheight_row = 20,
    treeheight_col = 20
  )

  dev.off()

  # Export correlation matrix to Excel
  cor_df <- as.data.frame(sample_cor)
  cor_df <- cbind(Sample = rownames(cor_df), cor_df)
  xlsx_file <- sub("\\.pdf$", ".xlsx", output_file)
  writexl::write_xlsx(cor_df, xlsx_file)

  message(paste("Correlation heatmap saved to:", output_file))
  message(paste("Correlation matrix saved to:", xlsx_file))
}


#' Define gene sets for pathway analysis
#'
#' @return List of gene vectors for different pathways
get_gene_sets <- function() {
  gene_sets <- list(
    wnt_ligands = c("WNT1", "WNT2", "WNT2B", "WNT3", "WNT3A", "WNT4", "WNT5A",
                    "WNT5B", "WNT6", "WNT7A", "WNT7B", "WNT8A", "WNT8B", "WNT9A",
                    "WNT9B", "WNT10A", "WNT10B", "WNT11", "WNT16"),

    wnt_receptors = c("LRP5", "LRP6", "FZD1", "FZD2", "FZD3", "FZD4", "FZD5",
                     "FZD6", "FZD7", "FZD8", "FZD9", "FZD10", "FZR1", "SMO",
                     "PORCN", "WLS"),

    pluripotency = c("POU5F1", "NANOG", "SOX2", "LIN28A", "OTX2", "DNMT3A",
                    "DNMT3B", "PRDM14", "DPPA2", "DPPA4", "DPPA5", "ESRRB",
                    "ZIC2", "ZIC3", "ZFP42"),

    wnt_targets = c("CCND1", "MYC", "PDK1", "FN1", "PTGS2", "AXIN2", "LEF1",
                   "CLDN1", "VEGFA", "FGF18", "ID2", "TERT", "LGR5", "FZD7",
                   "FST", "EN2", "GJA1", "MITF", "T", "TBX3", "BGLAP",
                   "NEUROG1", "SP5", "NEUROD1", "NKX2.2", "GBX2", "CDX1",
                   "CDC25A", "PITX2", "CDH1", "FGF4", "VCAN", "TNFRSF19"),

    wnt_linked = c("ALPK2", "AMFR", "ANKRD6", "CAV1", "CPE", "CSNK1G1",
                  "CSNK2A2", "CTHRC1", "CTNND1", "DAAM1", "DACT1", "DACT3",
                  "DKKL1", "DRD2", "EDA", "EGFR", "FGF10", "FZD1", "FZD10",
                  "GRK5", "GSC", "GSK3B", "HHEX", "KANK1", "KLF4", "KREMEN1",
                  "LATS1", "LRRFIP2", "LRRK1", "LRRK2", "MITF", "MLLT3",
                  "NFKB1", "NID1", "NOTUM", "NR4A2", "PLCG2", "PRICKLE2",
                  "RNF146", "RNF43", "ROR2", "RSPO1", "SCYL2", "SFRP4",
                  "STK3", "SULF1", "TMEM67", "TMEM88", "TNFAIP3", "TTC21B",
                  "UBE2B", "USP47", "VPS35", "WNT2", "WNT2B", "WNT3",
                  "WNT3A", "WNT7A", "WNT8A", "WNT9A", "ZNF703")
  )

  return(gene_sets)
}
