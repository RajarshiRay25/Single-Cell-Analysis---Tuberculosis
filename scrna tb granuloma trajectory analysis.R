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

View(seurat_obj@meta.data)

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




# Trajectory analysis using Monocle3 

# Convert Seurat object to Monocle object 

cds <- as.cell_data_set(seurat_obj_filter)
cds

# Observe the cell metadata of the scRNA object

View(colData(cds))

# Gene metadata - View the Gene candidates across the cell samples

fData(cds)
rownames(fData(cds))[1:10]

# Adding the gene short name column to cds to assign the gene names

fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)


# assign the data to a single partition in Monocle data to make a single trajectory

reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds@clusters$UMAP$partitions <- reacreate.partition

list_cluster <- seurat_obj_filter@active.ident
cds@clusters$UMAP$clusters <- list_cluster
list_cluster

# OR

# Let Monocle automatically assign cluster

cds <- cluster_cells(cds, reduction_method = "UMAP")  # Let Monocle3 find clusters


# Assign UMAP coordinate - cell embeddings

cds@int_colData@listData$reducedDims$UMAP <- seurat_obj_filter@reductions$umap@cell.embeddings


View(seurat_obj_filter@reductions$umap@cell.embeddings)


# plot

cds <- learn_graph(cds, use_partition = FALSE)


cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
                            color_cells_by = "CellTypeAnnotations",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan','orange','pink','coral','salmon')) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names


plot_cells(cds,
           color_cells_by = 'CellTypeAnnotations',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)


# Some code tests and debugs

table(clusters(cds))  # Verify presence of clusters

# Get the barcodes (cell names) for cluster 5
cluster5_cells <- colnames(cds[, clusters(cds) == 5])

# Extract the cell identities (Idents) for these cells
idents_cluster5 <- cds@colData[cluster5_cells, "CellTypeAnnotations"] 

# View the unique cell types in cluster 5
table(idents_cluster5)

# See cluster wise cell presence and number
table(clusters(cds), colData(cds)$CellTypeAnnotations)
 

# Order the cells by setting up the root through the cluster number we extracted in the last code

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 4]))

# Plot the trajectory

traj_1 <- plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)

traj_2 <-  plot_cells(cds,
           color_cells_by = 'CellTypeAnnotations',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

traj_1 | traj_2

# cells ordered by monocle3 pseudotime

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(CellTypeAnnotations, monocle3_pseudotime, median), fill = CellTypeAnnotations)) +
  geom_boxplot()

# Extracting the DEG from the pseudotime experiment to analyse which genes involve in the overall progression as per trajectory

deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

deg_bcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()

# Feature plot to observe the top 3 genes obtained from the above  code

FeaturePlot(seurat_obj_filter, features = c('ABCA1', 'ABCA3', 'ABCB11'))

