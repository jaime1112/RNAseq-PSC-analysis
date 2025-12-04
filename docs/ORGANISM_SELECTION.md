# Organism Selection Guide

**Date:** December 4, 2025
**Feature:** Flexible organism database selection for annotation and enrichment

---

## Overview

The pipeline now supports flexible organism selection for gene annotation and enrichment analysis. This allows you to automatically use the correct genome database based on your comparison type:

- **Human** (Homo sapiens) - uses `org.Hs.eg.db` and KEGG code `hsa`
- **Pig** (Sus scrofa) - uses `org.Ss.eg.db` and KEGG code `ssc`

---

## When to Use Each Organism

### Use `organism = "human"` when:

1. **Human-only comparisons**
   - hESC line A vs hESC line B
   - hIPSC vs hESC
   - Any comparison involving only human samples

2. **Cross-species comparisons (human vs pig)**
   - hPSC vs pPSC
   - hESC vs piPSC
   - Human cells vs pig cells
   - **Note:** Use human as the reference genome when comparing across species

### Use `organism = "pig"` when:

1. **Pig-only comparisons**
   - piPSC line A vs piPSC line B
   - pESC vs piPSC
   - Pig embryo vs pig fibroblast
   - Any comparison involving only pig samples

---

## Accepted Organism Identifiers

The system accepts multiple synonyms for each organism (case-insensitive):

### Human
- `"human"` (recommended)
- `"Hs"`
- `"hsa"`
- `"Homo sapiens"`

### Pig
- `"pig"` (recommended)
- `"Ss"`
- `"ssc"`
- `"Sus scrofa"`
- `"porcine"`

---

## How to Specify Organism

### Method 1: In Analysis Scripts

Update the `ORGANISM` parameter at the top of your analysis script:

```r
# In full_test_analysis.R or custom scripts
ORGANISM <- "human"  # or "pig"
```

### Method 2: In Configuration Files

Use the analysis configuration template:

```
# config/analysis_config_template.txt
organism: human  # or pig
```

### Method 3: In Function Calls

Pass directly to functions:

```r
# Annotation
res <- annotate_results(res, gene_symbols, organism = "human")

# Enrichment
enrichment <- run_enrichment_analysis(
  gene_lists = gene_lists,
  universe = all_genes,
  organism = "pig"
)
```

---

## What Changes Based on Organism

### 1. Gene Annotation (`annotate_results()`)

The function automatically selects the correct annotation database:

**Human (org.Hs.eg.db):**
```r
res <- annotate_results(res, gene_symbols, organism = "human")
# Uses: org.Hs.eg.db
# Retrieves: Human gene names, Entrez IDs, ENSEMBL IDs
```

**Pig (org.Ss.eg.db):**
```r
res <- annotate_results(res, gene_symbols, organism = "pig")
# Uses: org.Ss.eg.db
# Retrieves: Pig gene names, Entrez IDs, ENSEMBL IDs
```

### 2. Enrichment Analysis (`run_enrichment_analysis()`)

The function automatically configures both GO and KEGG enrichment:

**Human:**
```r
enrichment <- run_enrichment_analysis(
  gene_lists = gene_lists,
  universe = all_genes,
  organism = "human"
)
# GO: Uses org.Hs.eg.db for human GO terms
# KEGG: Uses organism code "hsa" for human pathways
```

**Pig:**
```r
enrichment <- run_enrichment_analysis(
  gene_lists = gene_lists,
  universe = all_genes,
  organism = "pig"
)
# GO: Uses org.Ss.eg.db for pig GO terms
# KEGG: Uses organism code "ssc" for pig pathways
```

---

## Example Workflows

### Example 1: Human-only Comparison

Comparing two human ESC lines:

**Data Configuration** (`config/data_sources_human_only.csv`):
```csv
dataset_name,filepath,remove_pattern
hESC_ours,data/hESC_ours/hESC_ours.txt,
hESC_Selmi,data/hESC_Selmi_2021/hESC_Selmi_2021.txt,
```

**Analysis Script:**
```r
# Set organism
ORGANISM <- "human"

# Load and process data
# ... (data loading code)

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)

# Annotate with human gene info
res <- annotate_results(res, gene_symbols, organism = ORGANISM)

# Enrichment with human pathways
enrichment <- run_enrichment_analysis(
  gene_lists = gene_lists,
  universe = all_genes,
  organism = ORGANISM
)
```

**Output:**
- Gene names from human genome
- GO terms from human biology
- KEGG pathways for human (hsa)

### Example 2: Pig-only Comparison

Comparing pig iPSC vs pig ESC:

**Data Configuration** (`config/data_sources_pig_only.csv`):
```csv
dataset_name,filepath,remove_pattern
piPSC_ours,data/piPSC_ours/piPSC_ours.txt,
pESC_Choi2024,data/pESC_Choi_2024/pESC_Choi_2024.txt,
```

**Analysis Script:**
```r
# Set organism
ORGANISM <- "pig"

# Load and process data
# ... (data loading code)

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)

# Annotate with pig gene info
res <- annotate_results(res, gene_symbols, organism = ORGANISM)

# Enrichment with pig pathways
enrichment <- run_enrichment_analysis(
  gene_lists = gene_lists,
  universe = all_genes,
  organism = ORGANISM
)
```

**Output:**
- Gene names from pig genome
- GO terms from pig biology
- KEGG pathways for pig (ssc)

### Example 3: Cross-species Comparison

Comparing human PSC vs pig PSC:

**Data Configuration** (`config/data_sources_full.csv`):
```csv
dataset_name,filepath,remove_pattern
hESC_ours,data/hESC_ours/hESC_ours.txt,
piPSC_ours,data/piPSC_ours/piPSC_ours.txt,
```

**Analysis Script:**
```r
# Use human as reference for cross-species comparison
ORGANISM <- "human"

# Load and process data
# Gene symbols are already standardized across species in the data files

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)

# Annotate with human gene info (reference organism)
res <- annotate_results(res, gene_symbols, organism = ORGANISM)

# Enrichment with human pathways (reference organism)
enrichment <- run_enrichment_analysis(
  gene_lists = gene_lists,
  universe = all_genes,
  organism = ORGANISM
)
```

**Output:**
- Gene names from human genome (reference)
- GO terms from human biology (reference)
- KEGG pathways for human (reference)

**Important:** For cross-species comparisons, use the species that represents your reference genome. In most cases, this is human.

---

## Function Reference

### `get_organism_db(organism)`

Helper function that returns organism-specific database information.

**Parameters:**
- `organism` - Organism identifier (human/pig and synonyms)

**Returns:**
- List with:
  - `org_db` - Annotation database object (org.Hs.eg.db or org.Ss.eg.db)
  - `kegg_code` - KEGG organism code ("hsa" or "ssc")
  - `species_name` - Full species name ("Human" or "Pig")
  - `annot_code` - Annotation code ("Hs" or "Ss")

**Example:**
```r
org_info <- get_organism_db("human")
# $org_db: org.Hs.eg.db object
# $kegg_code: "hsa"
# $species_name: "Human"
# $annot_code: "Hs"
```

### `annotate_results(res, gene_symbols, organism = "human")`

Annotates DESeq2 results with gene metadata.

**Parameters:**
- `res` - DESeq2 results object
- `gene_symbols` - Gene symbols from rownames
- `organism` - Organism identifier (default: "human")

**Returns:**
- Annotated results object with added columns:
  - `genename` - Full gene name
  - `entrezid` - Entrez gene ID
  - `ensembl` - ENSEMBL gene ID

**Example:**
```r
res <- annotate_results(res, rownames(dds), organism = "pig")
```

### `run_enrichment_analysis(..., organism = "human")`

Performs GO and KEGG enrichment analysis.

**Parameters:**
- `gene_lists` - Named list with "Upregulated" and "Downregulated" genes (Entrez IDs)
- `universe` - Vector of all genes tested (Entrez IDs)
- `organism` - Organism identifier (default: "human")
- `ontology` - GO ontology ("BP", "MF", or "CC", default: "BP")
- `pval_cutoff` - P-value cutoff (default: 0.05)

**Returns:**
- List with:
  - `GO` - GO enrichment results
  - `KEGG` - KEGG enrichment results
  - `organism_info` - Organism database information

**Example:**
```r
enrichment <- run_enrichment_analysis(
  gene_lists = list(
    Upregulated = up_genes,
    Downregulated = down_genes
  ),
  universe = all_genes,
  organism = "pig",
  ontology = "BP"
)
```

---

## Configuration Files

### Analysis Configuration Template

`config/analysis_config_template.txt`:

```
# Organism Selection
# Options: "human"/"Hs"/"hsa" OR "pig"/"Ss"/"ssc"
organism: human

# For human-only comparisons: organism: human
# For pig-only comparisons: organism: pig
# For cross-species comparisons: organism: human (use as reference)
```

### Data Source Configurations

**Human-only** (`config/data_sources_human_only.csv`):
```csv
dataset_name,filepath,remove_pattern
hESC_ours,data/hESC_ours/hESC_ours.txt,
hESC_Selmi,data/hESC_Selmi_2021/hESC_Selmi_2021.txt,
```

**Pig-only** (`config/data_sources_pig_only.csv`):
```csv
dataset_name,filepath,remove_pattern
piPSC_ours,data/piPSC_ours/piPSC_ours.txt,
pESC_Choi2024,data/pESC_Choi_2024/pESC_Choi_2024.txt,
```

**Cross-species** (`config/data_sources_full.csv`):
```csv
dataset_name,filepath,remove_pattern
hESC_ours,data/hESC_ours/hESC_ours.txt,
hESC_Selmi,data/hESC_Selmi_2021/hESC_Selmi_2021.txt,
piPSC_ours,data/piPSC_ours/piPSC_ours.txt,
pESC_Choi2024,data/pESC_Choi_2024/pESC_Choi_2024.txt,
```

---

## Error Handling

The system provides clear error messages for invalid organism identifiers:

```r
res <- annotate_results(res, genes, organism = "mouse")
# Error: Unknown organism: mouse
# Supported: 'Hs'/'human'/'hsa' (human) or 'Ss'/'pig'/'ssc' (pig)
```

---

## Confirmation Messages

The functions print confirmation messages showing which organism is being used:

```r
annotate_results(res, genes, organism = "human")
# Message: Annotating results with Human gene information

run_enrichment_analysis(..., organism = "pig")
# Message: Using organism: Pig (KEGG: ssc , Annotation: Ss )
```

---

## Tips and Best Practices

### 1. Cross-species Comparisons
Always use **human as the reference** when comparing human vs pig samples. This ensures:
- Consistent gene nomenclature
- Access to more comprehensive pathway databases
- Better annotation coverage

### 2. Organism Code Flexibility
You can use any accepted synonym. These are all equivalent:
```r
organism = "human"
organism = "Human"  # Case-insensitive
organism = "Hs"
organism = "hsa"
```

### 3. Configuration Consistency
Set the organism parameter in ONE place at the top of your script:
```r
ORGANISM <- "human"
# Then use ORGANISM variable throughout
annotate_results(..., organism = ORGANISM)
run_enrichment_analysis(..., organism = ORGANISM)
```

### 4. Check Your Data
Before running analysis, verify:
- Are all samples from one species? → Use that species
- Are samples from multiple species? → Use human as reference
- Are gene symbols already standardized? → Yes (done during data processing)

---

## Database Requirements

Ensure you have the required Bioconductor packages installed:

**For human analyses:**
```r
BiocManager::install("org.Hs.eg.db")
```

**For pig analyses:**
```r
BiocManager::install("org.Ss.eg.db")
```

**For cross-species:**
```r
# Install both
BiocManager::install(c("org.Hs.eg.db", "org.Ss.eg.db"))
```

---

## File Locations

**Utility functions:**
- `scripts/utils/enrichment_analysis.R` - Contains `get_organism_db()` and enrichment functions
- `scripts/utils/deseq_helpers.R` - Contains `annotate_results()` function

**Configuration:**
- `config/analysis_config_template.txt` - Analysis configuration template
- `config/data_sources_human_only.csv` - Human-only data config
- `config/data_sources_pig_only.csv` - Pig-only data config
- `config/data_sources_full.csv` - Cross-species data config

**Documentation:**
- `docs/ORGANISM_SELECTION.md` - This file

**Test scripts:**
- `full_test_analysis.R` - Demonstrates organism selection

---

## Backwards Compatibility

The default organism is "human", maintaining compatibility with existing code:

```r
# Old code still works (uses human by default)
annotate_results(res, genes)
run_enrichment_analysis(gene_lists, universe)

# New code can specify organism
annotate_results(res, genes, organism = "pig")
run_enrichment_analysis(gene_lists, universe, organism = "pig")
```

---

## Related Documentation

- `README.md` - Main project documentation
- `METHODS.md` - Computational methods
- `DEPENDENCIES.md` - Software requirements
- `PCA_DUAL_GROUPING.md` - PCA visualization options

---

**Status:** ✅ Tested and Production-Ready
**Last Updated:** December 4, 2025
