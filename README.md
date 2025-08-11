# SlimR: Marker-Based R Package for Single-Cell and Spatial-Transcriptomic Annotation

[![CRAN Version](https://img.shields.io/cran/v/SlimR?label=CRAN)](https://cran.r-project.org/package=SlimR)
[![GitHub License](https://img.shields.io/github/license/Zhaoqing-wang/SlimR?label=License)](https://github.com/Zhaoqing-wang/SlimR/blob/main/LICENSE)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/SlimR?text=Downloads)](https://cran.r-project.org/package=SlimR)
[![GitHub R package version](https://img.shields.io/github/r-package/v/Zhaoqing-wang/SlimR?label=GitHub&color=green)](https://github.com/Zhaoqing-wang/SlimR/releases)
[![GitHub commit activity](https://img.shields.io/github/commit-activity/w/Zhaoqing-wang/SlimR?label=Commit%20activity)](https://github.com/Zhaoqing-wang/SlimR/commits/main/)

## Overview

<img width="233.28" height="270" alt="Sticker" src="inst/Sticker.png" align="right">

SlimR is an R package designed for annotating single-cell and spatial-transcriptomics (ST) datasets. It supports the creation of a unified marker list, Markers_list, using sources including: the package's built-in curated species-specific cell type and marker reference databases (e.g., 'Cellmarker2', 'PanglaoDB', 'scIBD', 'TCellSI'), Seurat objects containing cell label information, or user-provided Excel tables mapping cell types to markers.

Based on the Markers_list, SlimR can calculate gene expression of different cell types and predict annotation information and calculate AUC ('Celltype_Calculate'), and annotate it ('Celltype_Annotation'), then verify it ('Celltype_Verification'). At the same time, it can calculate gene expression corresponding to the cell type to generate the corresponding annotation reference map for manual annotation (e.g., 'Heatmap', 'Combined Plot', 'Box plot').

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
   - [3.1 Calculate cell types](#31-calculate-cell-types)  
   - [3.2 Annotate cell types](#32-annotate-cell-types) 
   - [3.3 Verify cell types](#33-verify-cell-types) 

4. [Semi-Automated Annotation Workflow](#4-semi-automated-annotation-workflow)  
   - [4.1 Annotation Heatmap](#41-annotation-heatmap)  
   - [4.2 Annotation Combined Plot](#42-annotation-combined-plot) 
   - [4.3 Annotation Box Plot](#43-annotation-box-plot)

5. [Other functions provided by SlimR](#5-other-functions-provided-by-slimr)
6. [Conclusion](#6-conclusion)

---

## 1. Preparation
### 1.1 Installation
Option One: [![CRAN Version](https://img.shields.io/cran/v/SlimR?label=CRAN)](https://cran.r-project.org/package=SlimR)

Install SlimR directly from CRAN using: (Stable version, recommended when the version equivalent to GitHub package version)

```r
install.packages("SlimR")
```
*Note: Try adjusting the CRAN image to "Global (CDN)" or use "BiocManager::install("SlimR")" if you encounter a version mismatch during installation.*

Option Two: [![GitHub R package version](https://img.shields.io/github/r-package/v/Zhaoqing-wang/SlimR?label=GitHub&color=green)](https://github.com/Zhaoqing-wang/SlimR/releases)

Install SlimR directly from GitHub using: (Development version, more recommended when the version is higher than CRAN package version)

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
SlimR requires a standardized list format for storing marker information, metrics (can be omitted), and corresponding cell types (list names = cell types (necessary), first column = markers (necessary), subsequent columns = metrics (can be omitted)).

### 2.1 From Cellmarker2 Database
Cellmarkers2: A database of cell types and markers covering different species and tissue types. 

Reference: Hu et al. (2023) <doi:10.1093/nar/gkac947>.

#### 2.1.1 Load the database:
```r
Cellmarker2 <- SlimR::Cellmarker2
```
#### 2.1.2 Optional metadata exploration:
```r
Cellmarker2_table <- SlimR::Cellmarker2_table
View(Cellmarker2_table)
```
#### 2.1.3 Generate marker list:
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
**Important: Select at least the 'species' and 'tissue_class' parameters to ensure the accuracy of the annotation.**

*Note: Output usable in sections 3.1, 4.1, 4.2, 4.3 and 5.1.*

### 2.2 From PanglaoDB Database
PanglaoDB: Database of cell types and markers covering different species and tissue types. 

Reference: Franzén et al. (2019) <doi:10.1093/database/baz046>.

#### 2.2.1 Load the database:
```r
PanglaoDB <- SlimR::PanglaoDB
```
#### 2.2.2 Optional metadata exploration:
```r
PanglaoDB_table <- SlimR::PanglaoDB_table
View(PanglaoDB_table)
```
#### 2.2.3 Generate marker list:
```r
Markers_list_panglaoDB <- Markers_filter_PanglaoDB(
  PanglaoDB,
  species_input = 'Human',
  organ_input = 'GI tract'
)
```
**Important: Select the 'species_input' and 'organ_input' parameters to ensure the accuracy of the annotation.**

*Note: Output 'Markers_list' usable in sections 3.1, 4.1, 4.2, 4.3 and 5.2.*

### 2.3 From scIBD Database
scIBD: A database of human intestine markers. 

Reference: Nie et al. (2023) <doi:10.1038/s43588-023-00464-9>.
```r
Markers_list_scIBD <- SlimR::Markers_list_scIBD
```
**Important: This is for human intestinal annotation only. The input Seurat object was ensured to be a human intestinal type to ensure the accuracy of the labeling.**

*Note: Output 'Markers_list' usable in sections 3.1, 4.1, 4.2, 4.3 and 5.3*

### 2.4 From TCellSI Database
TCellSI: A database of T cell markers. 

Reference: Yang et al. (2024) <doi:10.1002/imt2.231>.
```r
Markers_list_TCellSI <- SlimR::Markers_list_TCellSI
```
**Important: This is only for T cell subset annotation. Ensure that the input Seurat object is of T cell type to guarantee the accuracy of the annotation.**

*Note: Output 'Markers_list' usable in sections 3.1, 4.1, 4.2, 4.3 and 5.4.*

### 2.5 From Seurat Objects
#### 2.5.1 First identify cluster features:
```r
seurat_markers <- FindAllMarkers(
  sce.all, 
  group.by = "Cell_type", 
  only.pos = TRUE
)
```
#### 2.5.2 Then generate marker list:
```r
Markers_list_Seurat <- read_seurat_markers(
  seurat_markers,
  sort_by = "avg_log2FC",
  gene_filter = 10
)
```
*Note: Output 'Markers_list' usable in sections 3.1, 4.1, 4.2, 4.3 and 5.3.*

### 2.6 From Excel Tables
**Format Requirements**:  

- Each sheet name = cell type  (necessary)

- First row = column headers  (necessary)

- First column = markers  (necessary)

- Subsequent columns = metrics  (can be omitted)

```r
Markers_list_Excel <- read_excel_markers("D:/Laboratory/Marker_load.xlsx")
```
*Note: Output 'Markers_list' usable in sections 3.1, 4.1, 4.2, 4.3 and 5.4.*

## 3. Automated Annotation Workflow
### 3.1 Calculate cell types
#### 3.1.1 Calculate Cell types (Core)
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
**Important: The parameter "cluster_col" in the function "Celltype_Calculate" and the function "Celltype_Annotation" must be strictly the same to avoid false matches.**

*Note: Using the parameter `AUC_correction = TRUE` takes a little longer to compute, but it is recommended to correct the predicted cell type this way.*

#### 3.1.2 Plot Heatmap (Optional)
Check the annotation probability of the cell type to be annotated in the input 'cluster_col' column and the cell type in 'Markers_list' with the following code.
```r
print(SlimR_anno_result$Heatmap_plot)
```
*Note: If the heatmap is not generated properly, please run the function "library(pheatmap)" first.*


#### 3.1.3 View prediction results (Optional)
Cell type information results predicted by SlimR can be viewed with the following code.
```r
View(SlimR_anno_result$Prediction_results)
```

#### 3.1.4 Plot ROC curve and AUC value (Optional)
Furthermore, the ROC curve and AUC value of the corresponding "cluster_col" and the predicted cell type can be viewed by the following code.
```r
print(SlimR_anno_result$AUC_plot)
```
**Improtant: This feature depends on 'plot_AUC = TRUE'**

*Note: If the heatmap is not generated properly, please run the function "library(ggplot2)" first.*

#### 3.1.5 Correction for predicted cell type (Alternative)
After viewing the list of predicted cell types and the corresponding AUC values, the predicted cell types can be corrected with the following code. (For example, cluster "15" in "cluster_col" corresponds to the predicted cell type "Intestinal stem cell")
```r
# For example, cluster "15" in "cluster_col" corresponds to the predicted cell type "Intestinal stem cell"
SlimR_anno_result$Prediction_results$Predicted_cell_type[
  SlimR_anno_result$Prediction_results$cluster_col == 15
] <- "Intestinal stem cell"
```

### 3.2 Annotate cell types
Assigns SlimR predicted cell types information in `SlimR_anno_result$Prediction_results$Predicted_cell_type` to the Seurat object based on cluster annotations, and stores the results into `seurat_obj@meta.data$annotation_col`.

```r
sce <- Celltype_Annotation(seurat_obj = sce,
    cluster_col = "seurat_clusters",
    SlimR_anno_result = SlimR_anno_result,
    plot_UMAP = TRUE,
    annotation_col = "Cell_type_SlimR"
    )
```
**Important: The parameter "cluster_col" in the function "Celltype_Calculate" and the function "Celltype_Annotation" must be strictly the same to avoid false matches. And the parameter "annotation_col" in the function "Celltype_Annotation" and the function "Celltype_Verification" must be strictly the same to avoid false matches.**

### 3.3 Verify cell types
By using the highly variable genes in `SlimR_anno_result$Expression_list` corresponding to predicted cell type information in `SlimR_anno_result$Prediction_results$Predicted_cell_type`, generate dotplot based on clusters as `seurat_obj@meta.data$annotation_col`.
```r
Celltype_Verification(seurat_obj = sce,
    SlimR_anno_result = SlimR_anno_result,
    gene_number = 5,
    colour_low = "white",
    colour_high = "navy",
    annotation_col = "Cell_type_SlimR"
    )
```
**Important: The parameter "annotation_col" in the function "Celltype_Annotation" and the function "Celltype_Verification" must be strictly the same to avoid false matches.**

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
*Note: Now this function has been incorporated into "Celltype_Calculate()", and it is recommended to use "Celltype_Calculate()" instead.*

### 4.2 Annotation Combined Plot
Generates per-cell-type expression dot plots with metric heatmap (when the metric information exists):
```r
Celltype_annotation_Combined(
  seurat_obj = sce,
  gene_list = Markers_list,
  gene_list_type = "Cellmarker2",
  species = "Human",
  save_path = "./SlimR/Celltype_annotation_Combined/"
  )
```
Each resulting combined image consists of a dot plot above and a heat map below (if mertic information present). Dot plots show the expression level and expression ratio relationship between cell type and corresponding markers. Below it, there is a metric heatmap for the corresponding marker (if the metric information exists).

### 4.3 Annotation Box Plot
Generates per-cell-type expression box plots:
```r
Celltype_annotation_Box(
  seurat_obj = sce,
  gene_list = Markers_list, 
  species = "Human",
  cluster_col = "seurat_cluster",
  assay = "RNA",
  save_path = "./SlimR/Celltype_annotation_Box/"
)
```
Each generated boxplot shows the box plot of the expression levels of the corresponding markers for that cell type, with the colors corresponding to the average expression levels of the markers.

## 5. Other functions provided by SlimR
Functions in section 5.1, 5.2, 5.3 and 5.4 has been incorporated into "Celltype_annotation_Combined()", and it is recommended to use "Celltype_annotation_Combined() instead.

#### 5.1 Dot plot With Cellmarker2 Database
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
#### 5.2 Dot plot With PanglaoDB Database
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
#### 5.3 Dot plot With Seurat-Based Marker Lists
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
#### 5.4 Dot plot With Excel-Based Marker Lists
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
*Note: This function also works with 'Markers_list' without mertic information or with mertic information generated in other ways.*

## 6. Conclusion
Thank you for using SlimR. For questions, issues, or suggestions, please contact:

**Zhaoqing Wang**  
📧 851091628@qq.com; zhaoqingwang@mail.sdu.edu.cn
