# RNAseq Analysis: Comparative Transcriptomics of Porcine and Human Pluripotent Stem Cells

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the bioinformatics analysis pipeline and code for comparative transcriptomic analysis of porcine and human pluripotent stem cells (PSCs), with a focus on WNT signaling pathway regulation.

## Publication [Details Pending]

**Title:** [Your Publication Title]

**Authors:** [Author List]

**Journal:** [Journal Name], [Year]

**DOI:** [Publication DOI]

## Project Structure

```
RNAseq/
├── README.md                           # Project overview
├── USER_GUIDE.md                       # Detailed usage instructions
├── LICENSE                             # MIT License
├── .gitignore                          # Git ignore rules
│
├── run_differential_expression.Rmd        # Main DE analysis (annotated notebook)
├── run_differential_expression_local.Rmd  # Working copy (edit freely)
├── run_upset_venn.Rmd                     # UpSet plot and Venn diagram generator
├── run_heatmap.Rmd                        # Heatmap generator
├── install_dependencies.R                 # One-shot package installer
│
├── scripts/
│   └── utils/                         # Shared utility functions
│       ├── data_processing.R          # Data loading and processing
│       ├── deseq_helpers.R            # DESeq2 analysis helpers
│       ├── plotting_themes.R          # Visualization functions
│       ├── enrichment_analysis.R      # GO/KEGG enrichment
│       ├── gene_set_operations.R      # Gene set operations and overlaps
│       └── heatmap_helpers.R          # Heatmap generation
│
├── config/                            # Configuration files
│   ├── data_sources_TEMPLATE.csv      # Dataset sources template
│   ├── sample_metadata_TEMPLATE.txt   # Sample metadata template
│   ├── WntSignalingPathway.txt        # WNT pathway gene list
│   ├── PositiveRegulationOfWntSignalingPathway.txt  # Positive WNT regulators
│   └── NegativeRegulationOfWntSignalingPathway.txt  # Negative WNT regulators
│
├── data/                              # Data directory (not tracked)
│   ├── README.md                      # Data availability information
│   └── [dataset subdirectories]
│
├── figures/                           # Generated outputs (not tracked)
│
├── docs/                              # Documentation
│   ├── METHODS.md                     # Computational methods
│   ├── DEPENDENCIES.md                # Software requirements
│   ├── MULTI_COMPARISON_WORKFLOW.md   # Advanced workflow guide
│   ├── ORGANISM_SELECTION.md          # Enrichment organism selection
│   └── PCA_DUAL_GROUPING.md          # PCA and sample correlation details
```

## Research Context

This project investigates the transcriptional differences between porcine and human pluripotent stem cells, with emphasis on:

- Differential gene expression patterns between species
- WNT signaling pathway activity and regulation
- Pluripotency marker expression
- Comparative analysis across multiple published PSC datasets

## Key Analyses

1. **Differential Expression Analysis**
   - DESeq2-based pairwise comparisons
   - Log2 fold change shrinkage with apeglm
   - Multiple testing correction (Benjamini-Hochberg)

2. **Principal Component Analysis (PCA) and Sample QC**
   - Sample clustering and quality control
   - Variance-stabilized transformation
   - Visualization of experimental groups
   - Sample-to-sample correlation heatmap

3. **Pathway Enrichment Analysis**
   - Gene Ontology (GO) enrichment
   - KEGG pathway analysis
   - WNT signaling pathway-specific analysis

4. **Visualization**
   - MA plots with pathway-specific gene highlighting
   - Heatmaps of differentially expressed genes
   - Sample correlation heatmaps
   - Venn diagrams for gene set overlaps
   - UpSet plots for multi-way comparisons


## Installation

### Prerequisites

- R version ≥ 4.2.0
- **RStudio** (highly recommended - download from https://posit.co/download/rstudio-desktop/)
- Git (optional)

### R Package Dependencies

**Recommended Setup:** Create an RStudio Project first:
1. In RStudio: `File` → `New Project` → `Existing Directory`
2. Browse to your RNAseq folder
3. Click "Create Project"

Install required packages in RStudio Console:

```r
# CRAN packages
install.packages(c(
  "tidyverse",
  "readxl",
  "writexl",
  "ggplot2",
  "ggrepel",
  "ggpubr",
  "pheatmap",
  "RColorBrewer",
  "VennDiagram",
  "UpSetR",
  "viridis",
  "rstudioapi",
  "rmarkdown",
  "knitr"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "org.Hs.eg.db",
  "org.Ss.eg.db",
  "AnnotationDbi",
  "clusterProfiler",
  "enrichplot"
))
```

## Usage

### Quick Start

**For detailed instructions, see [USER_GUIDE.md](USER_GUIDE.md)**

1. **Prepare your data:**
   - Place count data files in `data/` directory
   - Copy and configure `config/data_sources_TEMPLATE.csv` with dataset paths
   - Copy and configure `config/sample_metadata_TEMPLATE.txt` with sample information

2. **Run differential expression analysis (in RStudio):**

   - Open `run_differential_expression.Rmd` (or duplicate it to a `_local.Rmd`
     copy and edit there — the `_local.Rmd` file is git-tracked but treated as
     your working copy; keep the template at its defaults).
   - Edit the parameter chunk near the top.
   - Step through chunks with **`Ctrl/Cmd+Shift+Enter`** (run current chunk) or
     click the green ▶ arrow on each chunk. Plots render inline and in the
     RStudio **Plots** pane in real time, so you can supervise QC before the
     pipeline finishes.
   - To run end-to-end: **Run → Run All**, or **Knit** to render an HTML report.

3. **Generate visualizations:**

   - Open `run_upset_venn.Rmd` or `run_heatmap.Rmd`.
   - Edit the parameter chunk, then run chunks interactively as above.

> **Terminal use:** The `.Rmd` files can be batch-rendered with
> `Rscript -e 'rmarkdown::render("run_differential_expression.Rmd")'`, but the
> intended workflow is interactive use inside RStudio so you can monitor each
> step.

### Three Main Notebooks

This pipeline provides three annotated R Markdown notebooks with configurable parameters near the top of each file:

1. **`run_differential_expression.Rmd`** - Complete DE analysis pipeline
   - Loads and merges datasets
   - Runs DESeq2 analysis
   - Generates PCA plots and sample correlation heatmap
   - Generates MA plots (overall and pathway-specific)
   - Optional GO/KEGG enrichment
   - Exports results to Excel

2. **`run_upset_venn.Rmd`** - Gene set overlap analysis
   - Compare multiple DE results or create categories from single result
   - Generate UpSet plots and Venn diagrams
   - Export overlap gene lists

3. **`run_heatmap.Rmd`** - Flexible heatmap generation
   - Multiple gene selection methods (pathway, file, custom, all)
   - Sample filtering by metadata
   - Customizable clustering and colors

### Customization

All notebooks expose their parameters in a single chunk near the top — edit
the values, then run downstream chunks:

````
```{r parameters}
ANALYSIS_NAME       <- "GroupA_vs_GroupB"
BASELINE_LEVEL      <- "GroupB"
COMPARISON_LEVEL    <- "GroupA"
DESIGN_FORMULA      <- "~CellType"
FACTOR_NAME         <- "CellType"
ENRICHMENT_ORGANISM <- "human"
PERFORM_ENRICHMENT  <- TRUE
```
````

**See [USER_GUIDE.md](USER_GUIDE.md) for complete documentation of all parameters.**

## Analysis Details

### Differential Expression Workflow

1. **Data Processing**
   - Import gene count matrices
   - Remove species-specific genes (optional)
   - Aggregate duplicate gene entries
   - Merge datasets on common gene identifiers

2. **DESeq2 Analysis**
   - Create DESeqDataSet with experimental design
   - Estimate size factors and dispersions
   - Fit negative binomial GLM
   - Perform Wald tests for comparisons

3. **Result Annotation**
   - Shrink log2 fold changes (apeglm method)
   - Annotate with gene names and IDs
   - Filter by significance (padj < 0.05)

4. **Quality Control and Visualization**
   - PCA plots for sample clustering
   - Sample correlation heatmap for QC
   - MA plots with pathway gene highlighting
   - Export results to Excel format

### Gene Set Definitions

The pipeline includes predefined gene sets for WNT pathway analysis:

- **WNT Ligands:** WNT1-16 family members
- **WNT Receptors:** Frizzled receptors, LRP5/6, etc.
- **Pluripotency Markers:** POU5F1, NANOG, SOX2, etc.
- **WNT Target Genes:** MYC, CCND1, AXIN2, etc.

Custom gene sets can be added in `config/` directory as tab-delimited text files.

## Citation

If you use this code or data, please cite: [PUBLICATION PENDING]

```bibtex
@article{yourname2025,
  title={Your Publication Title},
  author={Author List},
  journal={Journal Name},
  year={2025},
  doi={10.xxxx/xxxxx}
}
```

## Contributing

This repository is primarily for publication reproducibility. For questions or issues:

1. Open an issue on GitHub
2. Contact the corresponding author: [lifangjack.chu@ucalgary.ca]

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Chu Laboratory members for experimental work
- Data contributors (see Data Availability section)
- Funding sources: [Funding agencies and grant numbers]
- Claude Code used for repository preparation

## Contact

**Corresponding Author:** [Name]
**Email:** [email@institution.edu]
**Lab Website:** [URL]
**Institution:** [Institution Name]

---

**Last Updated:** December 2025
