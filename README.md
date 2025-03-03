# 📌 Pipeline for Single-Cell RNA-Seq Analysis of *M. tuberculosis* Samples

## 🔬 Overview

This repository contains the R code for processing and analyzing single-cell RNA sequencing (scRNA-seq) data from the following datasets:

### **4-Week *M. tuberculosis* Granulomas Dataset**
The dataset has been sourced from the **Broadway Single Cell Portal**:

🔗 [Dataset Link](https://singlecell.broadinstitute.org/single_cell/study/SCP1749/cellular-ecology-of-m-tuberculosis-granulomas-4-week-dataset?cluster=4Week_ClusteringDF.csv\&spatialGroups=--\&annotation=donor_id--group--study\&subsample=all#study-visualize)

### **Functional Role of CD8+ Lymphocytes in Tuberculosis**
This study explores how the depletion of innate and/or adaptive CD8+ lymphocytes in macaques impacts *Mycobacterium tuberculosis* (Mtb) infection control, leading to increased granuloma numbers, lung inflammation, and bacterial burden.

🔗 [Dataset Link](https://singlecell.broadinstitute.org/single_cell/study/SCP642/cd8-lymphocytes-are-critical-for-early-control-of-tuberculosis-in-macaques#study-visualize)

This pipeline includes **data pre-processing, quality control, normalization, clustering, marker gene identification, visualization, and trajectory analysis** using the **Seurat** and **Monocle3** packages.

---

## 🛠️ Installation & Dependencies

Ensure you have R and the following packages installed:

```r
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

install.packages("devtools")
remotes::install_github('satijalab/seurat-wrappers')
devtools::install_github('cole-trapnell-lab/monocle3')

install.packages(c("Seurat", "SeuratData","patchwork", "ggplot2", "tidyverse", "gridExtra"))
```

Load the necessary libraries:

```r
library(Seurat)
library(SeuratWrappers)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(monocle3)
```

---

## 💀 Data Processing Workflow

### 🔹 1. Read the Data

```r
cts <- ReadMtx(mtx = "4Week_countsmatrix.mtx",
               features = "4Week_features.tsv",
               cells = "4Week_barcodes.tsv",
               feature.column = 1)
seurat_obj <- CreateSeuratObject(counts = cts)
```

### 🔹 2. Load Metadata

```r
metadata <- read.delim("metadata.txt", header = TRUE, sep = "\t")
rownames(metadata) <- metadata[,1]  # Set first column as rownames
metadata <- metadata[,-1]  # Remove the first column
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)
```

### 🔹 3. Quality Control (QC)

```r
seurat_obj$mitoPercent <- PercentageFeatureSet(seurat_obj, pattern = '^MT-')
seurat_obj_filter <- seurat_obj
```

### 🔹 4. Normalization & Feature Selection

```r
seurat_obj_filter <- NormalizeData(seurat_obj_filter, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj_filter <- FindVariableFeatures(seurat_obj_filter)
seurat_obj_filter <- ScaleData(seurat_obj_filter)
```

### 🔹 5. Clustering & UMAP Visualization

```r
seurat_obj_filter <- FindNeighbors(seurat_obj_filter, dims = 1:15)
seurat_obj_filter <- FindClusters(seurat_obj_filter)
seurat_obj_filter <- RunUMAP(seurat_obj_filter, dims = 1:15)
DimPlot(seurat_obj_filter, reduction = 'umap', label = TRUE)
```

### 🔹 6. Convert Seurat Object to Monocle3

```r
cds <- as.cell_data_set(seurat_obj_filter)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
```

### 🔹 7. Assign UMAP Coordinates

```r
cds@int_colData@listData$reducedDims$UMAP <- seurat_obj_filter@reductions$umap@cell.embeddings
```

### 🔹 8. Define Root Cell for Pseudotime Analysis

Selecting the root cell based on the biology of tuberculosis:

```r
cds <- order_cells(cds, root_cells = "Macrophage")
```

### 🔹 9. Visualizing Pseudotime Trajectory

```r
plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
```

---

## 📊 Results & Visualizations

### **🔸 UMAP Clustering of Cells**

This plot shows clusters of cells identified in the dataset.

```r
DimPlot(seurat_obj_filter, reduction = 'umap', label = TRUE)
```

### **🔸 Pseudotime Trajectory Analysis**

```r
plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
```

---

## 💡 Key Insights

✔️ **Dimensionality reduction & clustering** reveal distinct cell populations in the dataset.
✔️ **Differential gene expression (DGE) analysis** identifies key marker genes for different cell types.
✔️ **Monocle3-based pseudotime analysis** provides insights into the progression of cellular states in tuberculosis granulomas.

---

## 🤝 Acknowledgements

This analysis is based on data from the **Broadway Single Cell Portal**. Special thanks to the research community contributing to open-access single-cell transcriptomics datasets.

---

## 📝 License

This project is open-source.

📌 **Developer-1**: [Rajarshi Ray]\
📧 **Contact**: [rajarshi.ray@tuni.fi]

📌 **Developer-2**: [Ratul Bhowmik]\
📧 **Contact**: [ratul.bhowmik@tuni.fi]

📌 **Research Mentor**: [Dr. Ashok Aspatwar]\
📧 **Contact**: [ashok.aspatwar@tuni.fi]

