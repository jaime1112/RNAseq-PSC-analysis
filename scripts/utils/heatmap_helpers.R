# Heatmap Generation Functions
# Flexible heatmap creation with gene and sample selection

#' Generate flexible heatmap with gene and sample selection
#'
#' @param count_data Data frame or matrix of normalized counts (genes x samples)
#' @param gene_list Vector of genes to include (optional, uses all if NULL)
#' @param sample_list Vector of samples to include (optional, uses all if NULL)
#' @param metadata Data frame with sample annotations (optional)
#' @param output_file Path to output PDF
#' @param title Plot title
#' @param cluster_rows Cluster rows (default: TRUE)
#' @param cluster_cols Cluster columns (default: TRUE)
#' @param show_rownames Show gene names (default: TRUE)
#' @param show_colnames Show sample names (default: FALSE)
#' @param color_scheme Color scheme: "viridis", "blue_red", or custom palette
#' @param scale_method Scaling method: "row" (per gene), "column" (per sample), "none", or "0to1"
#' @param annotation_col Sample annotations to display
#' @param width Plot width (default: 10)
#' @param height Plot height (default: 8)
#' @return Path to output file
generate_flexible_heatmap <- function(count_data, gene_list = NULL, sample_list = NULL,
                                     metadata = NULL, output_file, title = NULL,
                                     cluster_rows = TRUE, cluster_cols = TRUE,
                                     show_rownames = TRUE, show_colnames = FALSE,
                                     color_scheme = "viridis", scale_method = "0to1",
                                     annotation_col = NULL, width = 10, height = 8) {

  # Convert to matrix if data frame
  if (is.data.frame(count_data)) {
    # Check if first column is gene names
    if (!is.numeric(count_data[, 1])) {
      rownames(count_data) <- count_data[, 1]
      count_data <- count_data[, -1]
    }
    count_data <- as.matrix(count_data)
  }

  # Subset genes
  if (!is.null(gene_list)) {
    available_genes <- intersect(gene_list, rownames(count_data))
    if (length(available_genes) == 0) {
      stop("No genes from gene_list found in count_data")
    }
    if (length(available_genes) < length(gene_list)) {
      missing <- length(gene_list) - length(available_genes)
      message(paste("Warning:", missing, "genes from gene_list not found in data"))
    }
    count_data <- count_data[available_genes, , drop = FALSE]
    message(paste("Using", length(available_genes), "genes from provided gene_list"))
  }

  # Subset samples
  if (!is.null(sample_list)) {
    available_samples <- intersect(sample_list, colnames(count_data))
    if (length(available_samples) == 0) {
      stop("No samples from sample_list found in count_data")
    }
    if (length(available_samples) < length(sample_list)) {
      missing <- length(sample_list) - length(available_samples)
      message(paste("Warning:", missing, "samples from sample_list not found in data"))
    }
    count_data <- count_data[, available_samples, drop = FALSE]
    message(paste("Using", length(available_samples), "samples from provided sample_list"))
  }

  # Scale data
  if (scale_method == "0to1") {
    # 0-1 normalization per gene
    scaled_data <- t(apply(count_data, 1, function(x) {
      if (max(x, na.rm = TRUE) > 0) {
        return(x / max(x, na.rm = TRUE))
      } else {
        return(x)
      }
    }))
  } else if (scale_method == "row") {
    # Z-score per gene
    scaled_data <- t(scale(t(count_data)))
  } else if (scale_method == "column") {
    # Z-score per sample
    scaled_data <- scale(count_data)
  } else {
    # No scaling
    scaled_data <- count_data
  }

  # Select color scheme
  if (color_scheme == "viridis") {
    colors <- viridis::viridis(100)
  } else if (color_scheme == "blue_red") {
    colors <- grDevices::colorRampPalette(c("blue", "white", "red"))(100)
  } else if (color_scheme == "blue_yellow") {
    colors <- grDevices::colorRampPalette(c("#313695", "#4575b4", "#abd9e9", "#ffffbf", "#fee090", "#f46d43", "#a50026"))(100)
  } else if (is.vector(color_scheme)) {
    colors <- color_scheme
  } else {
    stop("Invalid color_scheme. Use 'viridis', 'blue_red', 'blue_yellow', or provide color vector")
  }

  # Prepare metadata annotations
  if (!is.null(metadata) && !is.null(annotation_col)) {
    # Subset metadata to samples in heatmap
    metadata_subset <- metadata[colnames(scaled_data), annotation_col, drop = FALSE]

    # Create annotation colors automatically
    annotation_colors <- list()
    for (col in annotation_col) {
      unique_vals <- unique(metadata_subset[[col]])
      n_vals <- length(unique_vals)

      if (n_vals <= 12) {
        # Use predefined palette for categorical data
        palette <- RColorBrewer::brewer.pal(min(n_vals, 8), "Set2")
        if (n_vals > 8) {
          palette <- c(palette, RColorBrewer::brewer.pal(n_vals - 8, "Set3"))
        }
        names(palette) <- unique_vals
        annotation_colors[[col]] <- palette
      }
    }
  } else {
    metadata_subset <- NULL
    annotation_colors <- NULL
  }

  # Determine font size based on number of genes
  n_genes <- nrow(scaled_data)
  if (n_genes > 100) {
    fontsize_row <- 4
  } else if (n_genes > 50) {
    fontsize_row <- 6
  } else {
    fontsize_row <- 8
  }

  # Generate heatmap
  pdf(file = output_file, width = width, height = height)

  heatmap_plot <- pheatmap::pheatmap(
    mat = scaled_data,
    color = colors,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "ward.D2",
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    fontsize_row = fontsize_row,
    fontsize_col = 8,
    border_color = NA,
    main = title,
    annotation_col = metadata_subset,
    annotation_colors = annotation_colors,
    treeheight_row = 20,
    treeheight_col = 20
  )

  print(heatmap_plot)
  dev.off()

  message(paste("Heatmap saved to:", output_file))
  message(paste("Dimensions:", nrow(scaled_data), "genes x", ncol(scaled_data), "samples"))

  return(invisible(output_file))
}


#' Load normalized counts from file
#'
#' @param count_file Path to normalized counts file (Excel or RDS)
#' @return Data frame of normalized counts
load_normalized_counts <- function(count_file) {

  if (!file.exists(count_file)) {
    stop(paste("Count file not found:", count_file))
  }

  if (grepl("\\.xlsx$", count_file)) {
    counts <- readxl::read_excel(count_file)
    counts <- as.data.frame(counts)

    # Set gene names as rownames
    if (!is.numeric(counts[, 1])) {
      rownames(counts) <- counts[, 1]
      counts <- counts[, -1]
    }
  } else if (grepl("\\.rds$", count_file)) {
    counts <- readRDS(count_file)
  } else {
    stop(paste("Unsupported file format:", count_file))
  }

  message(paste("Loaded counts:", nrow(counts), "genes x", ncol(counts), "samples"))

  return(counts)
}


#' Get samples by metadata criteria
#'
#' @param metadata Data frame with sample metadata
#' @param criteria Named list of column-value pairs to filter by
#' @return Vector of sample names matching criteria
get_samples_by_criteria <- function(metadata, criteria) {

  # Start with all samples
  selected_samples <- rownames(metadata)

  # Apply each criterion
  for (col in names(criteria)) {
    if (!col %in% colnames(metadata)) {
      stop(paste("Column", col, "not found in metadata"))
    }

    values <- criteria[[col]]
    selected_samples <- selected_samples[metadata[selected_samples, col] %in% values]
  }

  message(paste("Selected", length(selected_samples), "samples matching criteria"))

  return(selected_samples)
}


#' Create multiple heatmaps from gene list collection
#'
#' @param count_data Normalized count matrix
#' @param gene_list_collection Named list of gene vectors
#' @param output_dir Directory for output files
#' @param base_filename Base name for output files
#' @param metadata Sample metadata (optional)
#' @param sample_list Samples to include (optional)
#' @param ... Additional arguments passed to generate_flexible_heatmap
#' @return Vector of output file paths
generate_multiple_heatmaps <- function(count_data, gene_list_collection, output_dir,
                                      base_filename = "heatmap", metadata = NULL,
                                      sample_list = NULL, ...) {

  output_files <- c()

  for (list_name in names(gene_list_collection)) {
    gene_list <- gene_list_collection[[list_name]]

    if (length(gene_list) == 0) {
      message(paste("Skipping", list_name, "- empty gene list"))
      next
    }

    output_file <- file.path(output_dir, paste0(base_filename, "_", list_name, ".pdf"))

    generate_flexible_heatmap(
      count_data = count_data,
      gene_list = gene_list,
      sample_list = sample_list,
      metadata = metadata,
      output_file = output_file,
      title = paste(list_name, "Expression"),
      ...
    )

    output_files <- c(output_files, output_file)
  }

  message(paste("\nGenerated", length(output_files), "heatmaps in", output_dir))

  return(invisible(output_files))
}
