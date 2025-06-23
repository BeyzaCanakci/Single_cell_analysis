# Single_cell_analysis

##scRNA-seq Data Analysis with Scanpy in Colab

This project contains Python code written to analyze single-cell RNA sequencing (scRNA-seq) data using the Scanpy library. The analysis pipeline covers data reading, preprocessing, dimensionality reduction (PCA, UMAP), clustering (Leiden), and differential gene expression analysis. Additionally, it includes pseudotime analysis to study cell developmental trajectories and pathway enrichment analysis.

## Installation

The following libraries are required to run the project:

- scanpy
- matplotlib
- seaborn
- numpy
- pandas
- leidenalg
- gseapy

You can install these libraries using pip:
## Usage

The code is designed to be run in a Jupyter Notebook or Google Colab environment. You can execute the code cells sequentially to perform all analysis steps.

1.  **Data Loading:** The project uses Scanpy's built-in `pbmc3k` dataset. If you wish to use your own data, you need to modify the data reading section accordingly.
2.  **Preprocessing:** Steps for calculating quality control metrics, filtering (removing low-quality cells and genes), normalization, and log transformation are applied.
3.  **Highly Variable Gene Selection:** The most variable genes are identified for downstream analysis.
4.  **Dimensionality Reduction:** PCA and UMAP techniques are used to reduce data dimensionality and visualize it in a 2D space.
5.  **Clustering:** Cells are grouped based on similar gene expression profiles using the Leiden algorithm.
6.  **Differential Gene Expression:** Marker genes are identified for each cluster.
7.  **Pseudotime Analysis:** Cell developmental or differentiation pathways are investigated using Diffusion Pseudotime (DPT).
8.  **Pathway Enrichment Analysis:** The Gseapy library is used to determine which biological pathways are active in genes that change along pseudotime.
9.  **Cell Type Specific Analysis:** Optionally, similar analysis steps can be repeated on a subset of a specific cell type.

## Analysis Outputs

Running the code will generate the following outputs:

-   Quality control plots (histograms of total RNA, number of genes, mitochondrial ratio)
-   Scatter plot of highly variable genes
-   UMAP plots (colored by clustering and cell types)
-   Heatmap of marker genes
-   Pseudotime UMAP plot
-   Heatmap of genes varying along pseudotime
-   Optional, cell type specific analysis output
-   CSV file containing the list of marker genes (`marker_genler_pbmc3k.csv`)
-   CSV file containing the list of genes ranked by correlation with pseudotime (`pseudotime_genes_ranked.csv`)
-   Output directories for GSEA (Gene Set Enrichment Analysis) results (`gsea_results`, optional `mono_pathways`)

## Contributing

If you would like to contribute to the project, please open a pull request.

