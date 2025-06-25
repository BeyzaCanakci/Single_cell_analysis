
# Seurat Single-Cell RNA-Seq Analysis Pipeline

This project provides an end-to-end workflow for single-cell RNA-seq data analysis using [Seurat](https://satijalab.org/seurat/), covering data preprocessing, quality control, clustering, cell type annotation, visualization, pseudotime analysis, and gene ontology enrichment.  
All steps are demonstrated in the script [`Seurat/seurat-single-cell-analysis.R`](Seurat/seurat-single-cell-analysis.R).

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Outputs](#outputs)
- [References](#references)

---

## Overview

This pipeline processes 10X Genomics matrix-formatted data, performing:

- Quality control (QC) and visualization
- Filtering based on QC metrics
- Normalization and variable gene selection
- Dimensionality reduction (PCA, UMAP)
- Clustering and marker gene identification
- Cell type annotation and visualization
- Cell cycle scoring
- Pseudotime analysis (Monocle 3)
- Gene Ontology (GO) enrichment analysis

## Requirements

- R (≥ 4.0.0)
- [Seurat](https://satijalab.org/seurat/) (≥ 4.0)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [dplyr](https://dplyr.tidyverse.org/)
- [monocle3](https://cole-trapnell-lab.github.io/monocle3/)
- [clusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/)
- [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)

Install missing R packages as needed:
```r
install.packages(c("Seurat", "ggplot2", "dplyr"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("monocle3", "clusterProfiler", "org.Hs.eg.db"))
```

## Usage

1. **Prepare your data:**  
   Place `matrix.mtx`, `genes.tsv`, and `barcodes.tsv` in the working directory or adjust file paths in the script.

2. **Run the analysis:**  
   Open and execute the script:  
   [`Seurat/seurat-single-cell-analysis.R`](Seurat/seurat-single-cell-analysis.R)

   You can run the script in RStudio or via R command line:
   ```r
   source("Seurat/seurat-single-cell-analysis.R")
   ```

3. **Review outputs:**  
   Visualizations and results will be displayed or saved as specified in the script.

## Pipeline Steps

1. **Load Data:**  
   - Reads 10X matrix-format data.

2. **Create Seurat Object:**  
   - Filtering genes/cells based on minimum expression.

3. **Quality Control:**  
   - Calculates mitochondrial gene percentage and log10 genes per UMI.
   - Generates violin and scatter plots for QC metrics.

4. **Filtering:**  
   - Removes cells/genes not passing QC thresholds.

5. **Normalization & Feature Selection:**  
   - Normalizes data and selects variable features.

6. **Dimensionality Reduction:**  
   - Performs scaling, PCA, and UMAP.

7. **Clustering:**  
   - Identifies cell clusters using graph-based methods.

8. **Cell Cycle Scoring:**  
   - Assigns cell cycle phase based on gene expression.

9. **Cluster Annotation:**  
   - Finds marker genes and annotates clusters using canonical markers.

10. **Pseudotime Analysis:**  
    - Performs trajectory inference with Monocle 3.

11. **GO Enrichment:**  
    - Runs GO analysis for marker genes.

## Outputs

- **Plots and Figures:**  
  - QC violin/scatter plots
  - UMAP and PCA plots
  - Marker heatmaps and feature plots
  - Barplots of cell type composition
  - Pseudotime plots
  - GO enrichment dotplots

- **Tables:**  
  - Cluster marker genes
  - Cell type fractions

## References

- [Seurat Documentation](https://satijalab.org/seurat/)
- [Monocle 3 Documentation](https://cole-trapnell-lab.github.io/monocle3/)
- [clusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/)

