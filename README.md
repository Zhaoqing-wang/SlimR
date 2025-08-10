# SlimR: Marker-Based R Package for Single-Cell and Spatial-Transcriptomic Annotation

[![CRAN Version](https://img.shields.io/cran/v/SlimR?label=CRAN)](https://cran.r-project.org/package=SlimR)
[![GitHub License](https://img.shields.io/github/license/Zhaoqing-wang/SlimR?label=License)](https://github.com/Zhaoqing-wang/SlimR/blob/main/LICENSE)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/SlimR?text=Downloads)](https://cran.r-project.org/package=SlimR)
[![GitHub R package version](https://img.shields.io/github/r-package/v/Zhaoqing-wang/SlimR?label=GitHub&color=green)](https://github.com/Zhaoqing-wang/SlimR/releases)
[![GitHub commit activity](https://img.shields.io/github/commit-activity/w/Zhaoqing-wang/SlimR?label=Commit%20activity)](https://github.com/Zhaoqing-wang/SlimR/commits/main/)

## Overview

<img width="233.28" height="270" alt="Sticker" src="inst/Sticker.png" align="right">

SlimR is an R package designed for annotating single-cell and spatial-transcriptomics (ST) datasets. It supports the creation of a unified marker list, Markers_list, using sources including: the package's built-in curated species-specific cell type and marker reference databases (e.g., 'Cellmarker2', 'PanglaoDB', 'scIBD', 'TCellSI'), Seurat objects containing cell label information, or user-provided Excel tables mapping cell types to markers.

Based on the Markers_list, SlimR can calculate gene expression of different cell types and predict annotation information and calculate AUC ('Celltype_Calculate') with one click, and annotate it ('Celltype_Annotation'). At the same time, it can calculate gene expression corresponding to the cell type to generate the corresponding annotation reference map for manual annotation (for example, 'Heatmap', 'Dot Plot', 'Box plot').

## Table of Contents
1. [Preparation](#1-preparation)  
   - [1.1 Installation](#11-installation)  
   - [1.2 Loading SlimR](#12-loading-slimr)  
   - [1.3 Dependencies (if installation fails)](#13-dependencies-if-installation-fails)  

2. [Standardized Marker_list Input](#2-standardized-marker_list-input)  
   - [2.1 From Cellmarker2 Database](#21-from-cellmarker2-database)  
   - [2.2 From PanglaoDB Database](#22-from-panglaodb-database)  
   - [2.3 From scIBD Database](#23-from-scibd-database)
   - [2.4 From TCellSI Database](#24-from-tcellsi-database)  
   - [2.5 From Seurat Objects](#25-from-seurat-objects)  
   - [2.6 From Excel Tables](#26-from-excel-tables)  
  
3. [Automated Annotation Workflow](#3-automated-annotation-workflow) 
   - [3.1 Calculate Celltype](#31-calculate-celltype)  
   - [3.2 Annotation Celltype](#32-annotation-celltype) 

4. [Semi-Automated Annotation Workflow](#4-semi-automated-annotation-workflow)  
   - [4.1 Annotation Heatmap](#41-annotation-heatmap)  
   - [4.2 Annotation Dot Plot](#42-annotation-dot-plot) 
   - [4.3 Annotation Box Plot](#43-annotation-box-plot)

5. [Other functions provided by SlimR](#5-other-functions-provided-by-slimr)
6. [Conclusion](#6-conclusion)

---

## 1. Preparation
### 1.1 Installation
Option One: [![CRAN Version](https://img.shields.io/cran/v/SlimR?label=CRAN)](https://cran.r-project.org/package=SlimR)

Install SlimR directly from CRAN using: (Stable version)

```r
install.packages("SlimR")
```
*Note: Try adjusting the CRAN image to "Global (CDN)" or use "BiocManager::install("SlimR")" if you encounter a version mismatch during installation.*

Option Two: [![GitHub R package version](https://img.shields.io/github/r-package/v/Zhaoqing-wang/SlimR?label=GitHub&color=green)](https://github.com/Zhaoqing-wang/SlimR/releases)

Install SlimR directly from GitHub using: (Development version, more recommended)

```r
devtools::install_github("Zhaoqing-wang/SlimR")
```

### 1.2 Loading SlimR
Load the package in your R environment:
```r
library(SlimR)
```

### 1.3 Dependencies (if installation fails)
SlimR requires R (≥ 3.5) and depends on the following packages: `cowplot`, `dplyr`, `ggplot2`, `patchwork`, `pheatmap`, `readxl`, `scales`, `Seurat`, `tidyr`, `tools`. Install missing dependencies using:
```r
# Install dependencies if needed:
install.packages(c("cowplot", "dplyr", "ggplot2", "patchwork", 
                   "pheatmap", "readxl", "scales", "Seurat", 
                   "tidyr", "tools"))
```

## 2. Standardized Marker_list Input
SlimR requires a standardized list format for storing marker information, metrics, and corresponding cell types (list names = cell types, first column = markers, subsequent columns = metrics).

### 2.1 From Cellmarker2 Database
Cellmarkers2: A database of cell types and markers covering different species and tissue types. Reference: Hu et al. (2023) <doi:10.1093/nar/gkac947>.

Load the database:
```r
Cellmarker2 <- SlimR::Cellmarker2
```
Optional metadata exploration:
```r
Cellmarker2_table <- SlimR::Cellmarker2_table
View(Cellmarker2_table)
```
Generate marker list:
```r
Markers_list_Cellmarker2 <- Markers_filter_Cellmarker2(
  Cellmarker2,
  species = "Human",
  tissue_class = "Intestine",
  tissue_type = NULL,
  cancer_type = NULL,
  cell_type = NULL
)
```
*Important: Select at least the 'species' and 'tissue_class' parameters to ensure the accuracy of the annotation. Note: Output usable in sections 3.1, 4.1, 4.2, 4.3 and 5.1.*

### 2.2 From PanglaoDB Database
PanglaoDB: Database of cell types and markers covering different species and tissue types. Reference: Franzén et al. (2019) <doi:10.1093/database/baz046>.

Load the database:
```r
PanglaoDB <- SlimR::PanglaoDB
```
Optional metadata exploration:
```r
PanglaoDB_table <- SlimR::PanglaoDB_table
View(PanglaoDB_table)
```
Generate marker list:
```r
Markers_list_panglaoDB <- Markers_filter_PanglaoDB(
  PanglaoDB,
  species_input = 'Human',
  organ_input = 'GI tract'
)
```
*Important: Select the 'species_input' and 'organ_input' parameters to ensure the accuracy of the annotation. Note: Output usable in sections 3.1, 4.1, 4.2, 4.3 and 5.2.*

### 2.3 From scIBD Database
scIBD: A database of human intestine markers. Reference: Nie et al. (2023) <doi:10.1038/s43588-023-00464-9>.
```r
Markers_list_scIBD <- SlimR::Markers_list_scIBD
```
*Important: This is for human intestinal annotation only. The input Seurat object was ensured to be a human intestinal type to ensure the accuracy of the labeling. Note: The output is available in Sections 3.1, 4.1, 4.2, 4.3 and 5.3*

### 2.4 From TCellSI Database
TCellSI: A database of T cell markers. Reference: Yang et al. (2024) <doi:10.1002/imt2.231>.
```r
Markers_list_TCellSI <- SlimR::Markers_list_TCellSI
```
*Important: This is only for T cell subset annotation. Ensure that the input Seurat object is of T cell type to guarantee the accuracy of the annotation. Note: Output usable in sections 3.1, 4.1, 4.2, 4.3 and 5.4.*

### 2.5 From Seurat Objects
First identify cluster features:
```r
seurat_markers <- FindAllMarkers(
  sce.all, 
  group.by = "Cell_type", 
  only.pos = TRUE
)
```
Then generate marker list:
```r
Markers_list_Seurat <- read_seurat_markers(
  seurat_markers,
  sort_by = "avg_log2FC",
  gene_filter = 10
)
```
*Note: Output usable in sections 3.1, 4.1, 4.2, 4.3 and 5.3.*

### 2.6 From Excel Tables
**Format Requirements**:  
- Each sheet name = cell type  
- First row = column headers  
- First column = markers  
- Subsequent columns = metrics  
```r
Markers_list_Excel <- read_excel_markers("D:/Laboratory/Marker_load.xlsx")
```
*Note: Output usable in sections 3.1, 4.1, 4.2, 4.3 and 5.4.*

## 3. Automated Annotation Workflow
### 3.1 Calculate Celltype
Uses "marker_list" to calculate probability, prediction results and generate heatmap for cell annotation.
```r
SlimR_anno_result <- Celltype_Calculate(seurat_obj = sce,
    gene_list = Markers_list,
    species = "Human",
    cluster_col = "seurat_clusters",
    assay = "RNA",
    min_expression = 0.1,
    specificity_weight = 3,
    compute_AUC = TRUE,
    plot_AUC = TRUE,
    AUC_correction = TRUE
    )
```
*Important: The parameter "cluster_col" in the function "Celltype_Calculate" and the function "Celltype_Annotation" must be strictly the same to avoid false matches.*

### 3.2 Annotation Celltype
Assigns SlimR predicted cell types to the Seurat object based on cluster annotations, and stores the results in the meta.data slot.
```r
sce <- Celltype_Annotation(seurat_obj = sce,
    cluster_col = "seurat_clusters",
    SlimR_anno_result = SlimR_anno_result,
    plot_UMAP = TRUE
    )
```
*Important: The parameter "cluster_col" in the function "Celltype_Calculate" and the function "Celltype_Annotation" must be strictly the same to avoid false matches.*


## 4. Semi-Automated Annotation Workflow
### 4.1 Annotation Heatmap
Generate a heatmap to estimate the likelihood that various cell clusters exhibited similarity to control cell types:
```r
Celltype_annotation_Heatmap(
  seurat_obj = sce,
  gene_list = Markers_list,
  species = "Human",
  cluster_col = "seurat_cluster",
  min_expression = 0.1,
  specificity_weight = 3
)
```

### 4.2 Annotation Dot Plot
Generates per-cell-type expression dot plots:
```r
Celltype_annotation_Dotplot(
  seurat_obj = sce,
  gene_list = Markers_list,
  gene_list_type = "Cellmarker2",
  species = "Human",
  save_path = "./annotations/"
  )
```

### 4.3 Annotation Box Plot
Generates per-cell-type expression box plots:
```r
Celltype_annotation_Box(
  seurat_obj = sce,
  gene_list = Markers_list, 
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_annotation_Bar/"
)
```

## 5. Other functions provided by SlimR
### 5.1 Dot plot With Cellmarker2 Database
```r
Celltype_annotation_Cellmarker2(
  seurat_obj = sce,
  gene_list = Markers_list_Cellmarker2,
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_annotation_Cellmarkers2/"
)
```

### 5.2 Dot plot With PanglaoDB Database
```r
Celltype_annotation_PanglaoDB(
  seurat_obj = sce,
  gene_list = Markers_list_panglaoDB,
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_annotation_PanglaoDB/"
)
```

### 5.3 Dot plot With Seurat-Based Marker Lists
```r
Celltype_annotation_Seurat(
  seurat_obj = sce,
  gene_list = Markers_list_Seurat,
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_annotation_Seurat/"
)
```

### 5.4 Dot plot With Excel-Based Marker Lists
Generates integrated dot plots and metric heatmaps:
```r
Celltype_annotation_Excel(
  seurat_obj = sce,
  gene_list = Markers_list_Excel,
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_annotation_Excel/"
)
```
## 6. Conclusion
Thank you for using SlimR. For questions, issues, or suggestions, please contact:

**Zhaoqing Wang**  
📧 851091628@qq.com; zhaoqingwang@mail.sdu.edu.cn
