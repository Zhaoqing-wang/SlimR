### SlimR: Marker-Based R Package for Single-Cell and Spatial Transcriptomics Annotation  

#### Overview  
SlimR is an R package designed for annotating single-cell and spatial transcriptomics datasets. It supports the creation of a unified marker list (`Markers_list`) using multiple sources:  
- User-provided Excel tables mapping cell types to markers
- Seurat objects containing cell label information
- Built-in curated species-specific reference databases (e.g., `Cellmarker2`, `PanglaoDB`)  

Leveraging the standardized `Markers_list`, SlimR generates an annotation heatmap (`Annotation_heatmap`) to visualize the relationship between input cell types and reference markers. It further facilitates systematic analysis of each cell type, producing cell-type-specific reference plots (e.g., `Markers_dotplot`, `Metric_heatmap`, `Mean_expression_bar_plot`).

---

#### Table of Contents  
1. [Prerequisites](#1-prerequisites)  
   - [1.1 Install SlimR](#11-install-slimr)  
   - [1.2 Load SlimR](#12-load-slimr)  
   - [1.3 Dependencies](#13-dependencies)  
2. [Standardized `Markers_list` Import](#2-standardized-markers_list-import)  
   - [2.1 Excel Tables](#21-excel-tables-custom-annotation-data)  
   - [2.2 Seurat Objects](#22-seurat-objects-single-cell-to-spatial-mapping)  
   - [2.3 CellMarker2 Database](#23-preprocessed-cellmarker2-database)  
   - [2.4 PanglaoDB Database](#24-preprocessed-panglaodb-database)  
3. [Automated Annotation Workflow](#3-automated-annotation-workflow)  
   - [3.1 Annotation Heatmap](#31-annotation-heatmap)  
   - [3.2 Annotation Bar Plot](#32-annotation-bar-plot)  
4. [Semi-Automated Annotation Workflow](#4-semi-automated-annotation-workflow)  
   - [4.1 Excel-Derived `Markers_list`](#41-using-excel-derived-markers_list)  
   - [4.2 Seurat-Derived `Markers_list`](#42-using-seurat-derived-markers_list)  
   - [4.3 CellMarker2-Derived `Markers_list`](#43-using-cellmarker2-derived-markers_list)  
   - [4.4 PanglaoDB-Derived `Markers_list`](#44-using-panglaodb-derived-markers_list)  
5. [Contact](#5-contact)  

---

### 1. Prerequisites <a name="1-prerequisites"></a>  
#### 1.1 Install SlimR <a name="11-install-slimr"></a>  
```r
devtools::install_github("Zhaoqing-wang/SlimR")
```

#### 1.2 Load SlimR <a name="12-load-slimr"></a>  
```r
library(SlimR)
```

#### 1.3 Dependencies <a name="13-dependencies"></a>  
SlimR requires **R (≥ 3.5)** and depends on:  
`cowplot`, `dplyr`, `ggplot2`, `patchwork`, `pheatmap`, `readxl`, `scales`, `Seurat`, `tidyr`, `tools`  
```r
# Install dependencies if needed:
install.packages(c("cowplot", "dplyr", "ggplot2", "patchwork", 
                   "pheatmap", "readxl", "scales", "Seurat", 
                   "tidyr", "tools"))
```

---

### 2. Standardized `Markers_list` Import <a name="2-standardized-markers_list-import"></a>  

#### 2.1 Excel Tables (Custom Annotation Data) <a name="21-excel-tables-custom-annotation-data"></a>  
**Format Requirements**:  
- Sheet name = Cell type  
- Column 1 = Marker genes  
- Columns 2+ = Metrics  

```r
Markers_list_Excel <- read_excel_markers("path/to/Marker_load.xlsx")
```  
> **Compatibility**: Sections 3.1, 3.2, 4.1  

#### 2.2 Seurat Objects (Single-Cell to Spatial Mapping) <a name="22-seurat-objects-single-cell-to-spatial-mapping"></a>  
```r
seurat_markers <- FindAllMarkers(
  sce.all, 
  group.by = "Cell_type", 
  only.pos = TRUE
)

Markers_list_Seurat <- read_seurat_markers(
  seurat_markers,
  sort_by = "avg_log2FC",
  gene_filter = 10
)
```  
> **Compatibility**: Sections 3.1, 3.2, 4.2  

#### 2.3 Preprocessed CellMarker2 Database <a name="23-preprocessed-cellmarker2-database"></a>  
```r
# Load database
Cellmarker2 <- SlimR::Cellmarker2

# Explore metadata
Cellmarker2_table <- SlimR::Cellmarker2_table
View(Cellmarker2_table)

# Filter markers
Markers_list_Cellmarker2 <- Markers_filter_Cellmarker2(
  Cellmarker2,
  species = "Human",
  tissue_class = "Intestine"
)
```  
> **Compatibility**: Sections 3.1, 3.2, 4.3  

#### 2.4 Preprocessed PanglaoDB Database <a name="24-preprocessed-panglaodb-database"></a>  
```r
# Load database
PanglaoDB <- SlimR::PanglaoDB

# Explore metadata
PanglaoDB_table <- SlimR::PanglaoDB_table
View(PanglaoDB_table)

# Filter markers
Markers_list_panglaoDB <- Markers_filter_PanglaoDB(
  PanglaoDB,
  species_input = 'Human',
  organ_input = 'GI tract'
)
```  
> **Compatibility**: Sections 3.1, 3.2, 4.4  

---

### 3. Automated Annotation Workflow <a name="3-automated-annotation-workflow"></a>  

#### 3.1 Annotation Heatmap <a name="31-annotation-heatmap"></a>  
```r
Celltype_annotation_Heatmap(
  seurat_obj = sce.all,
  gene_list = Markers_list,
  species = "Human",
  cluster_col = "RNA_snn_res.0.4"
)
```  

#### 3.2 Annotation Bar Plot <a name="32-annotation-bar-plot"></a>  
```r
Celltype_annotation_Bar(
  seurat_obj = sce.all,
  gene_list = Markers_list, 
  cluster_col = "seurat_cluster",
  save_path = "./SlimR/Celltype_annotation_Bar/"
)
```  

---

### 4. Semi-Automated Annotation Workflow <a name="4-semi-automated-annotation-workflow"></a>  

#### 4.1 Using Excel-Derived `Markers_list` <a name="41-using-excel-derived-markers_list"></a>  
```r
Celltype_annotation_Excel(
  seurat_obj = sce.all,
  gene_list = Markers_list_Excel,
  save_path = "./SlimR/Celltype_annotation_Excel/"
)
```  

#### 4.2 Using Seurat-Derived `Markers_list` <a name="42-using-seurat-derived-markers_list"></a>  
```r
Celltype_annotation_Seurat(
  seurat_obj = sce.all,
  gene_list = Markers_list_Seurat,
  save_path = "./SlimR/Celltype_annotation_Seurat/"
)
```  

#### 4.3 Using CellMarker2-Derived `Markers_list` <a name="43-using-cellmarker2-derived-markers_list"></a>  
```r
Celltype_annotation_Cellmarker2(
  seurat_obj = sce.all,
  gene_list = Markers_list_Cellmarker2,
  save_path = "./SlimR/Celltype_annotation_Cellmarkers2.0/"
)
```  

#### 4.4 Using PanglaoDB-Derived `Markers_list` <a name="44-using-panglaodb-derived-markers_list"></a>  
```r
Celltype_annotation_PanglaoDB(
  seurat_obj = sce.all,
  gene_list = Markers_list_panglaoDB,
  save_path = "./SlimR/Celltype_annotation_PanglaoDB/"
)
```  

---

### 5. Contact <a name="5-contact"></a>  
For questions or feedback:  
**Zhaoqing Wang**  
📧 851091628@qq.com  

---  
*Thank you for using SlimR!*
