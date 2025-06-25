# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Read 10X matrix-formatted data
expr <- ReadMtx(
  mtx = "matrix.mtx",
  features = "genes.tsv",
  cells = "barcodes.tsv"
)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = expr, project = "ZenodoPBMC", min.cells = 3, min.features = 200)
# min.cells = 3: Keep genes expressed in at least 3 cells
# min.features = 200: Keep cells with at least 200 expressed genes

# Extract raw count matrix
counts_data <- GetAssayData(seurat_obj, slot = "counts")
head(counts_data[, 1:5])

# -------- QUALITY CONTROL --------

# Calculate the percentage of mitochondrial gene expression
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", assay = "RNA") / 100

# Compute log10 ratio of genes per UMI (complexity measure)
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)

# View QC metrics after addition
counts_data_qc <- GetAssayData(seurat_obj, slot = "counts")
head(counts_data_qc[, 1:5])

# QC Visualization
VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Bottom-left: low nCount and nFeature = likely empty droplets
FeatureScatter(seurat_obj, feature1 = "percent.mt", feature2 = "nFeature_RNA")
# Top-right: unusually high values = potential doublets

# -------- FILTERING BASED ON QC THRESHOLDS --------

seurat_obj <- subset(seurat_obj,
                     subset = nCount_RNA >= 500 &
                       nFeature_RNA >= 250 &
                       percent.mt <= 0.2 &
                       log10GenesPerUMI >= 0.8
)

# Re-visualize QC after filtering
VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# -------- NORMALIZATION --------

seurat_obj <- NormalizeData(seurat_obj)
normalized_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
head(normalized_data[, 1:5])

# -------- VARIABLE GENE SELECTION --------

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat_obj), 10)

# Visualize variable features
plot1 <- VariableFeaturePlot(seurat_obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# -------- SCALING & PCA --------

seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
ElbowPlot(seurat_obj)
# ElbowPlot: helps determine the number of PCs to use

# -------- CLUSTERING --------

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:6)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)
table(Idents(seurat_obj))

# -------- UMAP VISUALIZATION --------

seurat_obj <- RunUMAP(seurat_obj, dims = 1:6)
DimPlot(seurat_obj, group.by = "orig.ident")
DimPlot(seurat_obj, label = TRUE, group.by = "orig.ident", split.by = "orig.ident") + NoLegend()

# View metadata
head(seurat_obj@meta.data)

# -------- CELL CYCLE SCORING --------

# Load pre-defined cell cycle gene lists from Seurat
cc_genes <- Seurat::cc.genes.updated.2019
s_genes <- cc_genes$s.genes
g2m_genes <- cc_genes$g2m.genes

# Assign cell cycle phase scores to cells
seurat_obj <- CellCycleScoring(seurat_obj,
                               s.features = s_genes,
                               g2m.features = g2m_genes,
                               set.ident = TRUE
)
head(seurat_obj[[]], 10)  # Check metadata for 'S.Score', 'G2M.Score', 'Phase'

# Visualize cell cycle phases
DimPlot(seurat_obj, group.by = "Phase", label = TRUE)
DimPlot(seurat_obj, label = TRUE, split.by = "Phase") + NoLegend()
FeaturePlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score"), min.cutoff = 'q10', sort.cell = TRUE)

# Final UMAP check for technical biases (e.g., phase separation)
DimPlot(seurat_obj, reduction = "umap", group.by = "Phase", label = TRUE)

# -------- CLUSTER ANNOTATION USING MARKER GENES --------

# Identify marker genes for all clusters
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

# Visualize marker genes
DoHeatmap(seurat_obj, features = top10$gene)
FeaturePlot(seurat_obj, features = top10$gene)
subset(markers, gene %in% c("CD3D", "MS4A1", "CD14"))  # Marker inspection

# -------- PCA CONTRIBUTION VISUALIZATION --------

columns <- c(paste0("PC_", 1:16), "ident", "umap_1", "umap_2")
pc_data <- FetchData(seurat_obj, vars = columns)
FeaturePlot(seurat_obj, features = paste0("PC_", 1:6), reduction = "umap")

# -------- MARKER-BASED CLUSTER IDENTITY CHECK --------

FeaturePlot(seurat_obj, features = c("CD3D", "MS4A1", "CD14", "FCGR3A"))  # Key lineage markers
FeaturePlot(seurat_obj, features = c("CD3D", "CD3E"))       # T cells
FeaturePlot(seurat_obj, features = "MS4A1")                 # B cells
FeaturePlot(seurat_obj, features = c("CD14", "LYZ"))        # Monocytes
FeaturePlot(seurat_obj, features = c("FCGR3A", "MS4A7"))    # Macrophage-like
FeaturePlot(seurat_obj, features = c("GNLY", "NKG7"))       # NK cells
FeaturePlot(seurat_obj, features = "PPBP")                  # Platelets

# Violin plots for gene expression by cluster
VlnPlot(seurat_obj, features = c("CD3D", "CD3E", "LYZ", "TYMP", "S100A9", "NKG7"), pt.size = 0)
VlnPlot(seurat_obj, features = c("CD3D", "CD3E", "LYZ", "TYMP", "S100A9", "NKG7"), group.by = "seurat_clusters", pt.size = 0)

# -------- CLUSTER RENAMING --------

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- SetIdent(seurat_obj, value = "seurat_clusters")
levels(seurat_obj)

new.cluster.ids <- c(
  "T cells",
  "Monocytes (Inflammatory)",
  "CD8+ T/NK",
  "T subpopulation",
  "Myeloid (Low Activity)",
  "Monocytes/DC",
  "Pre-monocytes",
  "Neutrophil-like"
)
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
levels(seurat_obj)
Idents(seurat_obj) %>% table()

# Final UMAP with renamed clusters
DimPlot(seurat_obj, label = TRUE, label.size = 5) + NoLegend()

# -------- DIFFERENTIAL EXPRESSION BY CELL TYPE --------

markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DoHeatmap(seurat_obj, features = top10$gene)

# -------- CONDITION-BASED COMPOSITION ANALYSIS --------
# If two samples (e.g., healthy vs diseased) are present

prop.table(table(seurat_obj$orig.ident, Idents(seurat_obj)), 1) %>% barplot()

# -------- BARPLOT FOR CELL TYPE FRACTIONS --------

cell_counts <- table(Idents(seurat_obj)) %>% as.data.frame()
colnames(cell_counts) <- c("CellType", "Count")
cell_counts$Fraction <- cell_counts$Count / sum(cell_counts$Count)

ggplot(cell_counts, aes(x = reorder(CellType, -Fraction), y = Fraction, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(x = "Cell Type", y = "Fraction", title = "Cell Type Composition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -------- PSEUDOTIME ANALYSIS WITH SLINGSHOT --------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("monocle3")
install.packages("monocle3")
library(Seurat)
library(monocle3)

# Take Seurat object
seurat_obj <- seurat_obj

# Monocle3 transformation
cds <- as.cell_data_set(seurat_obj)

# Cluster ID and UMAP coordination 
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# Select root node seç ( T cells)
cds <- order_cells(cds, root_cells = colnames(seurat_obj)[which(Idents(seurat_obj) == "T cells")[1]])

# Pseudotime 
plot_cells(cds, color_cells_by = "pseudotime")


# -------- CELL CYCLE PHASE DISTRIBUTION ACROSS CLUSTERS --------

table(Phase = seurat_obj$Phase, Cluster = Idents(seurat_obj)) %>% prop.table(2)

# -------- GO ENRICHMENT ANALYSIS FOR A CELL TYPE --------

library(clusterProfiler)
library(org.Hs.eg.db)

# Get marker genes for T cells
tcell_markers <- FindMarkers(seurat_obj, ident.1 = "T cells", min.pct = 0.25)
genes <- rownames(tcell_markers)[tcell_markers$p_val_adj < 0.05]

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform GO enrichment
ego <- enrichGO(gene = entrez_ids$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

# Visualize results
dotplot(ego, showCategory = 15) + ggtitle("GO Enrichment – T Cells")
