# 📌 Pipeline for Single-Cell RNA-Seq Analysis of *M. tuberculosis* samples.

## 🔬 Overview

This repository contains the R code for processing and analyzing single-cell RNA sequencing (scRNA-seq) data from the following datasets - 

**4-Week *****M. tuberculosis***** Granulomas Dataset**. The dataset has been sourced from the **Broadway Single Cell Portal**:

🔗 [Dataset Link](https://singlecell.broadinstitute.org/single_cell/study/SCP1749/cellular-ecology-of-m-tuberculosis-granulomas-4-week-dataset?cluster=4Week_ClusteringDF.csv\&spatialGroups=--\&annotation=donor_id--group--study\&subsample=all#study-visualize)

Functional role of CD8+ lymphocytes in tuberculosis. The study explores how the depletion of innate and/or adaptive CD8+ lymphocytes in macaques impacts Mycobacterium tuberculosis (Mtb) infection control, leading to increased granuloma numbers, lung inflammation, and bacterial burden.

🔗 [Dataset Link](https://singlecell.broadinstitute.org/single_cell/study/SCP642/cd8-lymphocytes-are-critical-for-early-control-of-tuberculosis-in-macaques#study-visualize)

This pipeline includes **data pre-processing, quality control, normalization, clustering, marker gene identification, and visualization** using the Seurat package.

---

## 🛠️ Installation & Dependencies

Ensure you have R and the following packages installed:

```r
install.packages(c("Seurat", "SeuratData", "patchwork", "ggplot2", "tidyverse", "gridExtra"))
```

Load the necessary libraries:

```r
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(gridExtra)
```

---

## 📂 Data Processing Workflow

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

### 🔹 5. Principal Component Analysis (PCA) & Visualization

```r
seurat_obj_filter <- RunPCA(seurat_obj_filter)
ElbowPlot(seurat_obj_filter)
VizDimLoadings(seurat_obj_filter, dims = 1:15, reduction = "pca")
DimHeatmap(seurat_obj_filter, dims = 1:10, cells = 500, balanced = TRUE)
PCAPlot(seurat_obj_filter)
```

### 🔹 6. Clustering & UMAP Visualization

```r
seurat_obj_filter <- FindNeighbors(seurat_obj_filter, dims = 1:15)
seurat_obj_filter <- FindClusters(seurat_obj_filter)
seurat_obj_filter <- RunUMAP(seurat_obj_filter, dims = 1:15)
DimPlot(seurat_obj_filter, reduction = 'umap', label = TRUE)
```

### 🔹 7. Differential Gene Expression (DGE) Analysis

```r
marker.genes <- FindAllMarkers(seurat_obj_filter,
               logfc.threshold = 0.1,
               min.pct = 0.25,
               only.pos = TRUE,
               test.use = 'wilcox',
               slot = 'data')

# Save results
write.csv(marker.genes, file = "seurat_markers.csv", row.names = FALSE)
```

### 🔹 8. Marker Gene Visualization

```r
gene_list <- unique(marker.genes$gene)[1:5]
FeaturePlot(seurat_obj_filter, features = gene_list, min.cutoff = 'q10', label = TRUE)
DotPlot(seurat_obj_filter, features = gene_list, cols = c("blue", "red"), dot.scale = 8)
VlnPlot(seurat_obj_filter, features = gene_list, group.by = 'CellTypeAnnotations', ncol = 3)
DoHeatmap(subset(seurat_obj_filter, downsample = 100), features = gene_list, size = 3)
```

### 🔹 9. Find Differentially Expressed Genes Between Two Cell Types

```r
cell_to_cell.marker <- FindMarkers(seurat_obj_filter, ident.1 = 'Macrophage', ident.2 = 'Neutrophil')
```

---

## 📊 Results & Visualizations

### **🔸 UMAP Clustering of Cells**

This plot shows clusters of cells identified in the dataset.

```r
DimPlot(seurat_obj_filter, reduction = 'umap', label = TRUE)
```

### **🔸 Expression of Selected Marker Genes**

**Heatmap of top marker genes across cell types:**

```r
DoHeatmap(subset(seurat_obj_filter, downsample = 100), features = gene_list, size = 3)
```

**Violin plots of gene expression across clusters:**

```r
VlnPlot(seurat_obj_filter, features = gene_list, group.by = 'CellTypeAnnotations', ncol = 3)
```

---

## 💡 Key Insights

✔️ **Dimensionality reduction & clustering** reveal distinct cell populations in the dataset. ✔️ **Differential gene expression (DGE) analysis** identifies key marker genes for different cell types. ✔️ **Visualization tools (UMAP, heatmaps, dot plots, violin plots)** provide clear insights into gene expression patterns.

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

