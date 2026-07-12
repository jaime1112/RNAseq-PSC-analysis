# Gene Set Enrichment Analysis Functions
# Functions for GO term and KEGG pathway enrichment

#' Get organism database and code based on organism identifier
#'
#' @param organism Organism identifier: "Hs"/"human"/"hsa" for human, "Ss"/"pig"/"ssc" for pig
#' @return List with org_db (annotation database) and kegg_code (for KEGG)
get_organism_db <- function(organism) {

  # Normalize organism identifier to lowercase
  org_lower <- tolower(organism)

  # Determine organism type
  if (org_lower %in% c("hs", "human", "hsa", "homo sapiens")) {
    org_db <- org.Hs.eg.db::org.Hs.eg.db
    kegg_code <- "hsa"
    species_name <- "Human"
    annot_code <- "Hs"
  } else if (org_lower %in% c("ss", "pig", "ssc", "sus scrofa", "porcine")) {
    org_db <- org.Ss.eg.db::org.Ss.eg.db
    kegg_code <- "ssc"
    species_name <- "Pig"
    annot_code <- "Ss"
  } else {
    stop(paste("Unknown organism:", organism,
               "\nSupported: 'Hs'/'human'/'hsa' (human) or 'Ss'/'pig'/'ssc' (pig)"))
  }

  message(paste("Using organism:", species_name, "(KEGG:", kegg_code, ", Annotation:", annot_code, ")"))

  return(list(
    org_db = org_db,
    kegg_code = kegg_code,
    species_name = species_name,
    annot_code = annot_code
  ))
}


#' Perform GO and KEGG enrichment for up/downregulated genes
#'
#' @param gene_lists Named list with "Upregulated" and "Downregulated" gene vectors (Entrez IDs)
#' @param universe Vector of all genes tested (Entrez IDs)
#' @param organism Organism identifier ("Hs"/"human"/"hsa" for human, "Ss"/"pig"/"ssc" for pig)
#' @param ontology GO ontology ("BP", "MF", or "CC")
#' @param pval_cutoff P-value cutoff (default: 0.05)
#' @return List with GO and KEGG enrichment results
run_enrichment_analysis <- function(gene_lists, universe, organism = "human",
                                   ontology = "BP", pval_cutoff = 0.05) {

  # Get organism-specific databases
  org_info <- get_organism_db(organism)

  # GO enrichment
  ora_combined <- clusterProfiler::compareCluster(
    geneClusters = gene_lists,
    fun = "enrichGO",
    universe = universe,
    OrgDb = org_info$org_db,
    ont = ontology,
    keyType = "ENTREZID",
    pAdjustMethod = "BH",
    pvalueCutoff = pval_cutoff
  )

  if (is.null(ora_combined)) {
    message("No GO enrichment found for any gene cluster — skipping GO downstream steps.")
  }

  # Add gene symbols to GO results
  ora_combined <- add_gene_symbols(ora_combined, org_info$org_db)

  # KEGG enrichment
  # enrichKEGG fetches pathway lists live from https://rest.kegg.jp/, which can
  # time out or be unreachable. Give it a longer download window and trap network
  # errors so a KEGG outage doesn't kill the whole pipeline.
  prev_timeout <- getOption("timeout")
  options(timeout = max(300, prev_timeout))
  on.exit(options(timeout = prev_timeout), add = TRUE)

  kegg_combined <- tryCatch(
    clusterProfiler::compareCluster(
      geneClusters = gene_lists,
      fun = "enrichKEGG",
      universe = universe,
      organism = org_info$kegg_code,
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff = pval_cutoff
    ),
    error = function(e) {
      message("KEGG enrichment failed (", conditionMessage(e),
              ") — skipping KEGG; GO results are unaffected.")
      NULL
    }
  )

  if (is.null(kegg_combined)) {
    message("No KEGG enrichment found for any gene cluster — skipping KEGG downstream steps.")
  }

  # Add gene symbols to KEGG results
  kegg_combined <- add_gene_symbols(kegg_combined, org_info$org_db)

  return(list(
    GO = ora_combined,
    KEGG = kegg_combined,
    organism_info = org_info
  ))
}


#' Add gene symbols to enrichment results
#'
#' @param enrichment_result compareCluster result object
#' @param org_db Organism annotation database
#' @return enrichment_result with added GeneSymbol column
add_gene_symbols <- function(enrichment_result, org_db) {

  # compareCluster returns NULL when no terms enrich in any cluster
  if (is.null(enrichment_result) ||
      is.null(enrichment_result@compareClusterResult) ||
      nrow(enrichment_result@compareClusterResult) == 0) {
    return(enrichment_result)
  }

  # Split GeneID strings into lists of Entrez IDs
  entrez_lists <- strsplit(enrichment_result@compareClusterResult$geneID, "/")

  # Map Entrez IDs to symbols
  symbol_lists <- lapply(entrez_lists, function(entrez_ids) {
    symbols <- AnnotationDbi::mapIds(
      org_db,
      keys = entrez_ids,
      keytype = "ENTREZID",
      column = "SYMBOL",
      multiVals = "first"
    )
    # Replace NA symbols with original Entrez ID
    symbols[is.na(symbols)] <- entrez_ids[is.na(symbols)]
    return(symbols)
  })

  # Collapse symbols back into slash-separated strings
  enrichment_result@compareClusterResult$GeneSymbol <- sapply(
    symbol_lists,
    function(x) paste(x, collapse = "/")
  )

  return(enrichment_result)
}


#' Generate enrichment dotplot
#'
#' @param enrichment_result Enrichment result object
#' @param output_file Path to output PDF
#' @param title Plot title
#' @param n_categories Number of categories to show (default: 10)
#' @param width Plot width (default: 8)
#' @param height Plot height (default: 8)
generate_enrichment_dotplot <- function(enrichment_result, output_file, title,
                                       n_categories = 10, width = 8, height = 8) {

  if (is.null(enrichment_result) ||
      is.null(enrichment_result@compareClusterResult) ||
      nrow(enrichment_result@compareClusterResult) == 0) {
    message(paste("Skipping dotplot — no enrichment results for:", title))
    return(invisible(NULL))
  }

  pdf(file = output_file, width = width, height = height)

  dotplot <- enrichplot::dotplot(
    enrichment_result,
    showCategory = n_categories,
    title = title
  ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  print(dotplot)
  dev.off()

  if (interactive()) try(print(dotplot), silent = TRUE)

  message(paste("Enrichment dotplot saved to:", output_file))
}
