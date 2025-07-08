# SlimR: Marker-Based R Package for Single-Cell and Spatial Transcriptomics Annotation

## Overview
SlimR is an R package designed for annotating single-cell and spatial transcriptomics datasets. It supports the creation of a unified marker list ("Markers_list") through multiple sources: 
- Built-in curated species-specific cell type and marker reference databases (e.g., CellMarker2, PanglaoDB)
- Seurat objects containing cell label information
- User-provided Excel tables mapping cell types to markers

Based on "Markers_list", SlimR enables:
- Iterative visualization of cell type annotation plots (Markers_dotplot, Metric_heatmap, Mean_expression_bar_plot)
- One-click generation of annotation heatmaps visualizing relationships between input cell types and reference markers

---

## 1. Preparation

### 1.1 Installation
Install SlimR directly from GitHub using:
```r
devtools::install_github("Zhaoqing-wang/SlimR")
```

### 1.2 Package Loading
Load SlimR in your R environment:
```r
library(SlimR)
```

### 1.3 Dependencies
SlimR requires R (≥ 3.5) and depends on these packages: 
`cowplot`, `dplyr`, `ggplot2`, `patchwork`, `pheatmap`, `readxl`, `scales`, `Seurat`, `tidyr`, `tools`. 

If installation fails, install dependencies first:
```r
install.packages("Package_names")
BiocManager::install("Package_names")
```

---

## 2. Standardized "Marker_list" Import

SlimR requires a standardized list format for marker information, metrics, and cell type annotations (Format: List names = cell types, Column 1 = markers, Subsequent columns = metrics). This can be generated from:

### 2.1 Excel File Import
Input format requirements:
- Sheet name = Cell type annotation
- First row = Header
- First column = Marker genes
- Subsequent columns = Metric statistics

Import with:
```r
Markers_list_Excel <- read_excel_markers("D:/Laboratory/Marker_load.xlsx")
```
*Note:* This Marker_list can be used in Sections 3.1, 3.2, and 4.1

### 2.2 Seurat Object Import
Extract cluster markers first:
```r
seurat_markers <- FindAllMarkers(sce.all, group.by = "Cell_type", only.pos = TRUE)
```
Generate standardized list:
```r
Markers_list_Seurat <- read_seurat_markers(
  seurat_markers,
  sort_by = "avg_log2FC",
  gene_filter = 10
)
```
*Note:* This Marker_list supports Sections 3.1, 3.2, and 4.2

### 2.3 CellMarker2 Database
Load preprocessed database:
```r
Cellmarker2 <- SlimR::Cellmarker2
```
Optional metadata queries:
```r
Cellmarker2_table <- SlimR::Cellmarker2_table  # Tissue type information
Cellmarker2_raw <- SlimR::Cellmarker2_raw      # Original data
```
Filter markers:
```r
Markers_list_Cellmarker2 <- Markers_filter_Cellmarker2(
  Cellmarker2,
  species = "Human",
  tissue_class = "Intestine"
)
```
*Note:* This Marker_list supports Sections 3.1, 3.2, and 4.3

### 2.4 PanglaoDB Database
Load database:
```r
PanglaoDB <- SlimR::PanglaoDB
```
Optional metadata queries:
```r
PanglaoDB_table <- SlimR::PanglaoDB_table  # Organ information
PanglaoDB_raw <- SlimR::PanglaoDB_raw      # Original data
```
Filter markers:
```r
Markers_list_PanglaoDB <- Markers_filter_PanglaoDB(
  PanglaoDB,
  species_input = "Human",
  organ_input = "GI tract"
)
```
*Note:* This Marker_list supports Sections 3.1, 3.2, and 4.4

---

## 3. Automated Annotation Pipeline

### 3.1 Annotation Heatmap
Visualize marker expression patterns across clusters:
```r
Celltype_annotation_Heatmap(
  seurat_obj = sce.all,
  gene_list = Markers_list,
  species = "Human",
  cluster_col = "seurat_cluster"
)
```

### 3.2 Annotation Bar Plot
Generate mean expression bar plots:
```r
Celltype_annotation_Bar(
  seurat_obj = sce.all,
  gene_list = Markers_list,
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_annotation_Bar/"
)
```

---

## 4. Semi-Automated Annotation Pipeline

### 4.1 Excel-Based Annotation
Generate integrated dotplot and metric heatmap visualizations:
```r
Celltype_annotation_Excel(
  seurat_obj = sce.all,
  gene_list = Markers_list_Excel,
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_annotation_Excel/"
)
```

### 4.2 Seurat-Based Annotation
For single-cell to spatial transcriptomics mapping:
```r
Celltype_annotation_Seurat(
  seurat_obj = sce.all,
  gene_list = Markers_list_Seurat,
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_annotation_Seurat/"
)
```

### 4.3 CellMarker2 Annotation
Integrated visualization using CellMarker2 database:
```r
Celltype_annotation_Cellmarker2(
  seurat_obj = sce.all,
  gene_list = Markers_list_Cellmarker2,
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_annotation_Cellmarker2/"
)
```

### 4.4 PanglaoDB Annotation
Integrated visualization using PanglaoDB database:
```r
Celltype_annotation_PanglaoDB(
  seurat_obj = sce.all,
  gene_list = Markers_list_PanglaoDB,
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_annotation_PanglaoDB/"
)
```

---

## 5. Contact
Thank you for using SlimR. For technical support or feedback, please contact:

Zhaoqing Wang  
[851091628@qq.com](mailto:851091628@qq.com)
