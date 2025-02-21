# Load required libraries
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# Set file paths (Modify these based on actual filenames)
mtx_file <- "sparse_all_data_cd8.mtx.gz"
features_file <- "features_all_data_cd8.txt"
barcodes_file <- "barcodes_all_data_cd8.txt"
metadata_file <- "formatted_metadata.csv"

# Load the count matrix
cts <- ReadMtx(mtx = mtx_file,
               features = features_file,
               cells = barcodes_file,
               feature.column = 1)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = cts)

# Load metadata
metadata <- read.csv(metadata_file, header = TRUE, sep = ",")
rownames(metadata) <- metadata[,1]  # Set first column as row names
metadata <- metadata[,-1]  # Remove the first column (since it's now rownames)

# Add metadata to Seurat object
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)
seurat_obj$sample <- rownames(seurat_obj@meta.data)

# Quality control: Mitochondrial percentage calculation
seurat_obj$mitoPercent <- PercentageFeatureSet(seurat_obj, pattern = '^MT-')

# Filter low-quality cells (adjust thresholds if needed)
seurat_obj <- subset(seurat_obj, subset = nCount_RNA > 800 & nFeature_RNA > 500)

# Normalize data
seurat_obj <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify variable features
seurat_obj <- FindVariableFeatures(object = seurat_obj)

# Scale data
seurat_obj <- ScaleData(object = seurat_obj)

# PCA for dimensionality reduction
seurat_obj <- RunPCA(object = seurat_obj)

# Determine optimal number of PCs using an elbow plot
ElbowPlot(seurat_obj)

# Find neighbors and clusters (adjust resolution if needed)
seurat_obj <- FindNeighbors(object = seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(object = seurat_obj, resolution = 0.5)

# Run UMAP for visualization
seurat_obj <- RunUMAP(object = seurat_obj, dims = 1:15)

# Plot UMAP
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# Explore metadata
View(seurat_obj@meta.data)

# Identify marker genes for each cluster
marker.genes <- FindAllMarkers(seurat_obj,
                               logfc.threshold = 0.1,
                               min.pct = 0.25,
                               only.pos = TRUE,
                               test.use = 'wilcox')

# Save marker genes to CSV
write.csv(marker.genes, file = "seurat_markers_new.csv", row.names = FALSE)

# Visualize top marker genes
top_genes <- head(VariableFeatures(seurat_obj), 10)
FeaturePlot(seurat_obj, features = top_genes, min.cutoff = 'q10')

# Annotate clusters using metadata if available
Idents(seurat_obj) <- seurat_obj@meta.data$General_Celltypes

unique(Idents(seurat_obj))


# Compare Macrophages/Monocytes vs. Neutrophils
cell_marker_diff <- FindMarkers(seurat_obj, 
                                ident.1 = "Macrophages/Monocytes", 
                                ident.2 = "Neutrophils")

# Save the results
write.csv(cell_marker_diff, file = "seurat_markers_macrophage_vs_neutrophil.csv", row.names = FALSE)

# View results
View(cell_marker_diff)


# Save cell-to-cell comparison results
cell_marker_diff <- cell_marker_diff %>% rownames_to_column(var = "Gene")
write.csv(cell_marker_diff, file = "seurat_markers_celltocell_new.csv", row.names = FALSE)

# Visualize marker gene expression
gene_list <- unique(marker.genes$gene)[1:5]
FeaturePlot(seurat_obj, features = gene_list, min.cutoff = 'q10', label = TRUE)

# Dot plot for expression visualization
DotPlot(seurat_obj, features = gene_list, cols = c("blue", "red"), dot.scale = 8)

# Heatmap for marker genes
DoHeatmap(subset(seurat_obj, downsample = 100), features = gene_list, size = 3)

# Violin plot to check expression levels
VlnPlot(seurat_obj, features = gene_list, group.by = 'General_Celltypes', ncol = 3)

