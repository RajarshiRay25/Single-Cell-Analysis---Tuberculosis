# Load the history

seurat_obj_filter <- readRDS("data files for tb granuloma/seurat_tb_granuloma.rds")

# Load Libraries

library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(monocle3)

# Observe the dataset

View(seurat_obj_filter@meta.data$CellTypeAnnotations)
unique(seurat_obj_filter@meta.data$CellTypeAnnotation)

# Upload dataset

counts <- ReadMtx(mtx = "data files for tb granuloma/4Week_countsmatrix.mtx",
                  features = "data files for tb granuloma/4Week_features.tsv",
                  cells = "data files for tb granuloma/4Week_barcodes.tsv",
                  feature.column = 1)

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = counts)

metadata <- read.delim("data files for tb granuloma/metadata.txt", header = TRUE, sep = "\t", row.names = 1)

rownames(counts)

# metadata <- metadata[rownames(metadata) != "TYPE", , drop = FALSE]

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

# Monocle3

cds <- as.cell_data_set(seurat_obj_filter)
cds

# Cell metdata

colData(cds)

# Gene metadata

fData(cds)
rownames(fData(cds))[1:10]

# since it misses the gene_short_name column, let's add it

fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)


# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds@clusters$UMAP$partitions <- reacreate.partition


list_cluster <- seurat_obj_filter@active.ident
cds@clusters$UMAP$clusters <- list_cluster

list_cluster


# Assign UMAP coordinate - cell embeddings

cds@int_colData@listData$reducedDims$UMAP <- seurat_obj_filter@reductions$umap@cell.embeddings

# plot

cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'CellTypeAnnotations',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")



cluster.before.trajectory

cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by = 'CellTypeAnnotations',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)
# ...4. Order the cells in pseudotime -------------------

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 16]))

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)

# cells ordered by monocle3 pseudotime

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(CellTypeAnnotations, monocle3_pseudotime, median), fill = CellTypeAnnotations)) +
  geom_boxplot()

# ...5. Finding genes that change as a function of pseudotime --------------------
deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

deg_bcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()

FeaturePlot(seurat_obj_filter, features = c('ABCA1', 'ABCA3', 'ABCB11'))
