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
├── run_differential_expression.R      # Main DE analysis script
├── run_upset_venn.R                   # UpSet plot and Venn diagram generator
├── run_heatmap.R                      # Heatmap generator
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
│   ├── data_sources_full.csv          # Dataset sources configuration
│   ├── sample_metadata_full.txt       # Sample metadata
│   └── WntSignalingPathway.txt        # WNT pathway gene list
│
├── data/                              # Data directory (not tracked)
│   ├── README.md                      # Data availability information
│   └── [dataset subdirectories]
│
├── figures/                           # Generated outputs (not tracked)
│
├── docs/                              # Documentation
│   └── MULTI_COMPARISON_WORKFLOW.md   # Advanced workflow guide
│
└── old_analysis_files/                # Archived analysis files (not tracked)
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

2. **Principal Component Analysis (PCA)**
   - Sample clustering and quality control
   - Variance-stabilized transformation
   - Visualization of experimental groups

3. **Pathway Enrichment Analysis**
   - Gene Ontology (GO) enrichment
   - KEGG pathway analysis
   - WNT signaling pathway-specific analysis

4. **Visualization**
   - MA plots with pathway-specific gene highlighting
   - Heatmaps of differentially expressed genes
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
# Install rstudioapi for RStudio integration
install.packages("rstudioapi")

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
  "viridis"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "org.Hs.eg.db",
  "org.Ss.eg.db",
  "clusterProfiler",
  "enrichplot"
))
```

## Usage

### Quick Start

**For detailed instructions, see [USER_GUIDE.md](USER_GUIDE.md)**

1. **Prepare your data:**
   - Place count data files in `data/` directory
   - Configure `config/data_sources_full.csv` with dataset paths
   - Configure `config/sample_metadata_full.txt` with sample information

2. **Run differential expression analysis:**

   **In RStudio:**
   - Open `run_differential_expression.R`
   - Edit the parameters at the top of the script
   - Click **"Source"** button or press `Ctrl/Cmd+Shift+S`

   **Or from Terminal:**
   ```bash
   Rscript run_differential_expression.R
   ```

3. **Generate visualizations:**

   **In RStudio:**
   - Open `run_upset_venn.R` or `run_heatmap.R`
   - Edit parameters as needed
   - Click **"Source"** to run

   **Or from Terminal:**
   ```bash
   Rscript run_upset_venn.R    # For UpSet plots
   Rscript run_heatmap.R       # For heatmaps
   ```

### Three Main Scripts

This pipeline provides three user-friendly scripts with configurable parameters at the top of each file:

1. **`run_differential_expression.R`** - Complete DE analysis pipeline
   - Loads and merges datasets
   - Runs DESeq2 analysis
   - Generates PCA and MA plots
   - Optional GO/KEGG enrichment
   - Exports results to Excel

2. **`run_upset_venn.R`** - Gene set overlap analysis
   - Compare multiple DE results or create categories from single result
   - Generate UpSet plots and Venn diagrams
   - Export overlap gene lists

3. **`run_heatmap.R`** - Flexible heatmap generation
   - Multiple gene selection methods (pathway, file, custom, all)
   - Sample filtering by metadata
   - Customizable clustering and colors

### Customization

All scripts use a parameter section at the top for easy configuration. Simply edit the values and run:

```r
# Example from run_differential_expression.R
ANALYSIS_NAME <- "pPSC_vs_hPSC"
BASELINE_LEVEL <- "hPSC"
COMPARISON_LEVEL <- "pPSC"
ORGANISM <- "human"  # or "pig"
PERFORM_ENRICHMENT <- TRUE
```

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

4. **Visualization and Export**
   - PCA plots for QC
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

## Contact

**Corresponding Author:** [Name]
**Email:** [email@institution.edu]
**Lab Website:** [URL]
**Institution:** [Institution Name]

---

**Last Updated:** December 2025
