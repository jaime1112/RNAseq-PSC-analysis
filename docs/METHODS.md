# Computational Methods

This document provides detailed descriptions of the computational methods used in this RNAseq analysis project.

## Table of Contents

1. [Data Processing](#data-processing)
2. [Differential Expression Analysis](#differential-expression-analysis)
3. [Normalization and Transformation](#normalization-and-transformation)
4. [Statistical Testing](#statistical-testing)
5. [Enrichment Analysis](#enrichment-analysis)
6. [Visualization Methods](#visualization-methods)

---

## Data Processing

### Gene Count Matrices

**Input Format:**
- Tab-delimited text files
- Column 1: Gene identifiers (SYMBOL or ENSEMBL)
- Remaining columns: Raw counts per sample

**Processing Steps:**

1. **Gene ID Standardization**
   - Remove version numbers from ENSEMBL IDs (text after underscore)
   - Standardize gene symbols to primary identifiers

2. **Species-Specific Filtering** (Optional)
   - Remove genes matching species-specific patterns
   - Example: Remove porcine genes (`^ENSSSCG`) from human comparisons

3. **Duplicate Gene Aggregation**
   - Sum counts for genes with multiple entries
   - Ensures unique gene identifiers per matrix

4. **Dataset Merging**
   - Inner join on gene identifiers
   - Retains only genes present in all datasets
   - Alternative: Full outer join with NA handling

**Implementation:** See `scripts/utils/data_processing.R`

---

## Differential Expression Analysis

### DESeq2 Method

We use DESeq2 (Love et al., Genome Biology 2014) for differential expression analysis.

**Model:**
```
K_ij ~ NB(μ_ij, α_i)
log2(μ_ij) = Σ_r x_jr β_ir
```

Where:
- K_ij: count for gene i, sample j
- μ_ij: expected mean count
- α_i: dispersion parameter for gene i
- x_jr: design matrix
- β_ir: coefficients

### Analysis Pipeline

1. **DESeqDataSet Creation**
   ```r
   dds <- DESeqDataSetFromMatrix(
     countData = counts,
     colData = metadata,
     design = ~CellType
   )
   ```

2. **Dispersion Estimation**
   - Per-gene dispersion estimates
   - Fit dispersion-mean relationship
   - Shrink gene-wise dispersions toward fitted values

3. **Negative Binomial GLM Fitting**
   - Maximum likelihood estimation
   - Beta-binomial distribution for counts

4. **Hypothesis Testing**
   - Wald test for pairwise comparisons
   - H0: log2(fold change) = 0
   - Alternative: Likelihood ratio test (LRT) for complex designs

5. **Log Fold Change Shrinkage**
   - Method: apeglm (Zhu et al., Bioinformatics 2019)
   - Reduces noise from low-count genes
   - Improves ranking and visualization

**Reference Level:**
- The first level of the factor is the reference (denominator)
- Positive log2FC: higher in comparison vs. reference
- Negative log2FC: lower in comparison vs. reference

---

## Normalization and Transformation

### Size Factor Normalization

DESeq2 uses median-of-ratios method:

1. Calculate geometric mean of each gene across samples
2. Compute ratio of each sample to geometric mean
3. Take median of ratios as size factor
4. Normalized count = raw count / size factor

**Purpose:** Account for sequencing depth differences

### Variance Stabilizing Transformation (VST)

For visualization and clustering:

```r
vsd <- vst(dds, blind = FALSE)
```

**Properties:**
- Stabilizes variance across mean count range
- Log2-like transformation for high counts
- Useful for PCA, heatmaps, sample QC

**Alternative:** rlog transformation (more conservative, slower)

---

## Statistical Testing

### Multiple Testing Correction

**Method:** Benjamini-Hochberg (FDR control)

Adjusted p-value (padj):
```
padj_i = min(p_i × n / rank(p_i), 1)
```

**Significance Threshold:**
- Adjusted p-value < 0.05 (default)
- Corresponds to 5% false discovery rate

### Independent Filtering

DESeq2 automatically filters genes to maximize power:
- Removes genes with very low counts
- Increases statistical power by reducing multiple testing burden
- Filter threshold determined to optimize number of rejections

### Effect Size Thresholds

In addition to statistical significance, we often filter by:
- |log2 fold change| > 1 (2-fold change)
- |log2 fold change| > 2 (4-fold change)

---

## Enrichment Analysis

### Gene Ontology (GO) Enrichment

**Method:** Over-representation analysis (ORA) via clusterProfiler

**Test:** Hypergeometric test
```
P(X ≥ k) = Σ [(M choose i)(N-M choose n-i)] / (N choose n)
```

Where:
- N: total genes in background
- M: genes in GO term
- n: genes in test set
- k: overlap between test set and GO term

**Background:** All genes with count data (universe)

**Ontologies:**
- BP: Biological Process
- MF: Molecular Function
- CC: Cellular Component

**Multiple Testing:** Benjamini-Hochberg correction

### KEGG Pathway Enrichment

Similar hypergeometric test for KEGG pathways:
- Organism: "hsa" (human) or "ssc" (pig)
- Database: KEGG pathway annotations
- Updates: Requires internet connection for current pathways

### Pathway Analysis Output

For each enriched term:
- **ID:** GO/KEGG identifier
- **Description:** Term name
- **GeneRatio:** k/n (genes in test set in term)
- **BgRatio:** M/N (genes in background in term)
- **pvalue:** Raw p-value
- **p.adjust:** Adjusted p-value (BH method)
- **qvalue:** Q-value (positive FDR)
- **GeneID:** Entrez IDs of overlapping genes
- **GeneSymbol:** Gene symbols of overlapping genes

---

## Visualization Methods

### Principal Component Analysis (PCA)

**Input:** Variance-stabilized counts

**Method:**
1. Center data (mean = 0 for each gene)
2. Compute covariance matrix
3. Calculate eigenvectors and eigenvalues
4. Project samples onto top principal components

**DESeq2 Implementation:**
- Uses top 500 most variable genes by default
- PC1 vs PC2 typically shows largest variation

**Interpretation:**
- Distance between points: transcriptional similarity
- Clustering: similar expression profiles
- Outliers: potential quality issues

### MA Plots

**Axes:**
- x-axis: log2(mean normalized counts) - "M" for "mean"
- y-axis: log2(fold change) - "A" for "abundance"

**Features:**
- Significant genes: colored points
- Horizontal line at y = 0: no change
- Fold change thresholds: horizontal lines

**Gene Highlighting:**
- Custom gene sets overlaid in red
- Labels for genes of interest
- ggrepel for non-overlapping labels

### Heatmaps

**Method:** pheatmap package

**Clustering:**
- Distance: correlation distance
- Method: Ward's linkage (ward.D2)

**Scaling:**
- Row scaling: z-score or 0-1 normalization
- Visualizes relative expression across samples

**Color Schemes:**
- viridis: perceptually uniform, colorblind-friendly
- RdBu: diverging palette for fold changes

### Venn Diagrams and UpSet Plots

**Venn Diagrams:**
- Up to 3-4 way comparisons
- Shows overlaps between gene lists

**UpSet Plots:**
- Better for >4 comparisons
- Shows intersection sizes as bars
- More scalable than Venn diagrams

---

## Quality Control

### Sample QC Metrics

1. **Total Read Counts**
   - Check for outliers
   - Assess library size distribution

2. **PCA Analysis**
   - Identify batch effects
   - Detect outlier samples
   - Verify experimental groupings

3. **Sample-to-Sample Distances**
   - Hierarchical clustering
   - Identify mismatched samples

4. **Gene Detection**
   - Number of genes with >0 counts
   - Number of genes passing filters

### Filtering Criteria

**Gene Filtering:**
- Minimum count threshold (e.g., ≥10 counts total)
- Present in minimum number of samples
- Implemented automatically by DESeq2

**Sample Filtering:**
- Remove low-quality samples if identified
- Document reasons for exclusion

---

## Software and Versions

All analyses use open-source software available via CRAN and Bioconductor.

**Core Packages:**
- DESeq2: Differential expression analysis
- apeglm: Log fold change shrinkage
- clusterProfiler: Enrichment analysis
- ggplot2: Visualization
- pheatmap: Heatmap generation

**Annotation:**
- org.Hs.eg.db: Human gene annotations
- org.Ss.eg.db: Pig gene annotations

See session info in analysis reports for exact versions used.

**Last Updated:** December 2025
