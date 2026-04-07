# Organism and Enrichment Database Selection

## Overview

The pipeline uses an `ENRICHMENT_ORGANISM` parameter to select the annotation database for gene annotation (Entrez ID mapping) and pathway enrichment (GO/KEGG).

## Recommended Setting

```r
ENRICHMENT_ORGANISM <- "human"
```

**Use `"human"` in most cases**, even for pig samples. Reasons:
- Gene symbols in the pipeline are standardized to human nomenclature
- Human databases (`org.Hs.eg.db`) have far better gene symbol-to-Entrez ID coverage
- Human GO/KEGG databases have much richer pathway annotations (including stem cell pluripotency, WNT signaling, etc.)

Use `"pig"` only if your gene symbols are pig-specific identifiers.

## Accepted Identifiers

| Human | Pig |
|-------|-----|
| `"human"` (recommended) | `"pig"` (recommended) |
| `"Hs"`, `"hsa"` | `"Ss"`, `"ssc"` |

## What It Controls

1. **Gene annotation** (`annotate_results()`) — selects `org.Hs.eg.db` vs `org.Ss.eg.db` for mapping gene symbols to Entrez IDs, gene names, and ENSEMBL IDs
2. **GO enrichment** — selects the organism database for GO term annotations
3. **KEGG enrichment** — selects the KEGG organism code (`"hsa"` vs `"ssc"`)

## Required Packages

```r
BiocManager::install("org.Hs.eg.db")  # For human
BiocManager::install("org.Ss.eg.db")  # For pig
```

---

**Last Updated:** April 2026
