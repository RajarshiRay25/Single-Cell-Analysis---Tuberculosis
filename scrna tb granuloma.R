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
VizDimLoadings(seurat_obj_filter, dims = 1:2, reduction = "pca")

DimHeatmap(seurat_obj_filter,dims = 1:10, cells = 500, balanced = TRUE)

PCAPlot(object = seurat_obj_filter)

# Perform Clustering - Find neighbor cells, cluster together and UMAP

seurat_obj_filter <- FindNeighbors(object = seurat_obj_filter,dims = 1:15)
seurat_obj_filter <- FindClusters(object = seurat_obj_filter)

seurat_obj_filter <- RunUMAP(object = seurat_obj_filter, dims = 1:15)

plot1 <- DimPlot(seurat_obj_filter,reduction = "umap",label = TRUE)
DimPlot(seurat_obj_filter,reduction = "umap",label = TRUE)


# plot
p1 <- DimPlot(seurat_obj_filter, reduction = 'umap', group.by = 'sex')
p1

# Meta Data

View(seurat_obj_filter@meta.data)


# Feature Extraction and Annotation

# Find all the DEG biomarkers in the dataset across each cell types

# seurat_obj_filter@assays$RNA$counts <- seurat_obj_filter@assays$RNA$counts + 1

FindAllMarkers(seurat_obj_filter,
               logfc.threshold = 0.1,
               min.pct = 0.25,
               only.pos = TRUE,
               test.use = 'wilcox',
               slot = 'data')
FeaturePlot(seurat_obj_filter,features = c('GPX1'),min.cutoff = 'q10')


# Observe the identities of the object

Idents(seurat_obj_filter)

# Assign Identities to names

View(seurat_obj_filter@meta.data)

Idents(seurat_obj_filter) <- seurat_obj_filter@meta.data$CellTypeAnnotations

Idents(seurat_obj_filter)

DimPlot(seurat_obj_filter,reduction = 'umap', label = T)

FeaturePlot(seurat_obj_filter,features = c('GPX1'),min.cutoff = 'q10',label = T)

DotPlot(seurat_obj_filter, features = c('GPX1'), cols = c("blue", "red"), dot.scale = 8)