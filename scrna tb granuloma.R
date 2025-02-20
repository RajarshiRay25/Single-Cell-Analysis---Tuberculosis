# Load Libraries

library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# Access data files

cts <- ReadMtx(mtx = "4Week_countsmatrix.mtx",
               features = "4Week_features.tsv",
               cells = "4Week_barcodes.tsv",
               feature.column = 1)

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = cts)

metadata <- read.delim("metadata.txt", header = TRUE, sep = "\t")
head(metadata)
View(seurat_obj@meta.data)

rownames(metadata) <- metadata[,1]  # Set first column as rownames
metadata <- metadata[,-1]  # Remove the first column (since it's now rownames)



# Add metadata to Seurat object
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)

seurat_obj$sample <- rownames(seurat_obj@meta.data)

View(seurat_obj)

# Perform mitochondrial quantification through pattern matching prior QC - non human data may not have

seurat_obj$mitoPercent <- PercentageFeatureSet(seurat_obj,pattern = '^MT-')

# Perform QC filtering 

# seurat_obj_filter <- subset(seurat_obj, subset = nCount_RNA > 800 & nFeature_RNA > 500)
seurat_obj_filter <- seurat_obj
View(seurat_obj@meta.data)

# Perform the standard preprocessing

seurat_obj_filter <- NormalizeData(object = seurat_obj_filter,normalization.method = "LogNormalize", scale.factor = 10000)

seurat_obj_filter <- FindVariableFeatures(object = seurat_obj_filter)

seurat_obj_filter <- ScaleData(object = seurat_obj_filter)

# See the variable features plot

# Find the top 10 variable genes

variable.genes.10 <- head(VariableFeatures(seurat_obj_filter),10)

# Plot the Variable genes 

plot.1 <- VariableFeaturePlot(seurat_obj_filter)
LabelPoints(plot = plot.1,points = variable.genes.10,repel = TRUE)


# Perform Dimensionality Test through PCA and Elbow 

seurat_obj_filter <- RunPCA(object = seurat_obj_filter)

ElbowPlot(seurat_obj_filter)
VizDimLoadings(seurat_obj_filter, dims = 1:15, reduction = "pca")

DimHeatmap(seurat_obj_filter,dims = 1:10, cells = 500, balanced = TRUE)

PCAPlot(object = seurat_obj_filter)

# Perform Clustering - Find neighbor cells, cluster together and UMAP

seurat_obj_filter <- FindNeighbors(object = seurat_obj_filter,dims = 1:15)
seurat_obj_filter <- FindClusters(object = seurat_obj_filter)

seurat_obj_filter <- RunUMAP(object = seurat_obj_filter, dims = 1:15)

plot1 <- DimPlot(seurat_obj_filter,reduction = "umap",label = TRUE)
DimPlot(seurat_obj_filter,reduction = "umap",label = TRUE)


# plot
p1 <- DimPlot(seurat_obj_filter, reduction = 'umap', group.by = 'donor_id')
p1

# Meta Data

View(seurat_obj_filter@meta.data)


# Feature Extraction and Annotation

# Find all the unique DEG biomarkers in the dataset across each cell types clusters

# seurat_obj_filter@assays$RNA$counts <- seurat_obj_filter@assays$RNA$counts + 1

marker.genes <- FindAllMarkers(seurat_obj_filter,
               logfc.threshold = 0.1,
               min.pct = 0.25,
               only.pos = TRUE,
               test.use = 'wilcox',
               slot = 'data')

# Save the results as a CSV file
write.csv(marker.genes, file = "seurat_markers.csv", row.names = FALSE)

# Create Graphical representations for marker specific maps across cell types

FeaturePlot(seurat_obj_filter,features = c('GPX1'),min.cutoff = 'q10')
FeaturePlot(seurat_obj_filter, features = c('GPX1'))


# Observe the identities of the object

Idents(seurat_obj_filter)

# Assign Identities to names

View(seurat_obj_filter@meta.data)

Idents(seurat_obj_filter) <- seurat_obj_filter@meta.data$CellTypeAnnotations

Idents(seurat_obj_filter)

DimPlot(seurat_obj_filter,reduction = 'umap', label = T)


# Read the marker list

# Load the CSV file with the list of marker genes
markers <- read.csv("seurat_markers.csv")

# Extract the unique gene names from the column containing gene names
gene_list <- unique(markers$gene)

# Extract the specific number of genes you want to see
gene_list <- gene_list[1:5]

# Create Graphical representations for marker specific maps across cell types - Annotated

FeaturePlot(seurat_obj_filter,features = gene_list, min.cutoff = 'q10',label = T)

# Create boxplot to assess the expression levels with percentage of marker genes across cell types

DotPlot(seurat_obj_filter, features = gene_list, cols = c("blue", "red"), dot.scale = 8)

# Create Graphical representations for marker specific maps across cell types by donor ID - Annotated

FeaturePlot(seurat_obj_filter,features = gene_list, min.cutoff = 'q10',label = T,split.by = 'donor_id')

# Create violin plots to assess the expression levels of the gene list across cell 

VlnPlot(seurat_obj_filter,features = gene_list, group.by = 'CellTypeAnnotations',ncol = 3)

# Single cell heatmap of feature expression
DoHeatmap(subset(seurat_obj_filter, downsample = 100), features = gene_list, size = 3)

# Find marker features and genes across 2 cell types as per choice

cell_to_cell.marker <- FindMarkers(seurat_obj_filter, ident.1 = 'Macrophage', ident.2 = 'Neutrophil')

# Observe the significant biomarkers between the 2 cell types

View(cell_to_cell.marker)

# Save the cell to cell biomarkers
# Convert row names (genes) into a separate column

cell_to_cell.marker <- cell_to_cell.marker %>%  rownames_to_column(var = "Gene")
write.csv(cell_to_cell.marker, file = "seurat_markers_celltocell.csv", row.names = T)

# Visualise the biomarkers

# Load the CSV file with the list of marker genes
markers_celltocell <- read.csv("seurat_markers_celltocell.csv")

# Extract the unique gene names from the column containing gene names
gene_list_celltocell <- unique(markers_celltocell$Gene)

# Extract the specific number of genes you want to see
gene_list_celltocell <- gene_list_celltocell[1:5]

FeaturePlot(seurat_obj_filter,features = gene_list_celltocell, min.cutoff = 'q10',label = T)





