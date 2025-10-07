# SlimR 1.0.8 (2025-10-08)

-   This version adds the function of machine learning (e.g., 'Random Forest', 'Gradient Boosting', 'Support Vector Machine', 'Ensemble Learning') for cell types calculate parameter recognition.
-   Optimize the data filter mode of "Markers_list_scIBD" in the package, and filter through `sort_by = "logFC"` and `gene_filter = 20` parameter.
-   Adjust the calculation process of the 'FSS' value in the `read_seurat_markers()` function when 'resources' is set to 'presto'.
-   Optimize the prompt output during the execution of the `Celltype_Verification()` function.
-   Modify and optimize README and NEWS file.

# SlimR 1.0.7 (2025-08-19)

-   Added new function `Celltype_Verification()` for predicted cell types validation and generate the validation dotplot.
-   Optimize the function 'Read_seurat_markers()'. This is compatible with the 'presto::wilcoxauc()' source tag, and the 'Feature Significance Score' (FSS, product value of `log2FC` and `Expression ratio`) can be calculated and sorted accordingly.
-   Add custom color parameters `colour_low` and `colour_high` to all ploting output functions.
-   Renamed `Celltype_annotation_Dotplot()` to `Celltype_Annotation_Features()`, `Celltype_annotation_Box()` to `Celltype_Annotation_Combined()`, `read_seurat_markers()` to `Read_seurat_markers()`, `read_excel_markers()` to `Read_excel_markers()` for unified function naming structure.
-   Enhanced README with detailed process descriptions.
-   Optimized message output system for cleaner console feedback.
-   Resolved various code bugs reported by users.
-   Modified codebase to meet CRAN standards and policies.

# SlimR 1.0.6 (2025-08-06)

-   Integrated "scIBD" human intestine reference database.
-   Added AUC calculation and visualization to `Celltype_Calculate()`.
-   Implemented AUC-based prediction correction in cell typing.
-   Streamlined code output formatting.
-   Fixed critical bugs in prediction pipeline.
-   Modified code to meet CRAN submission requirements.

# SlimR 1.0.5 (2025-08-05)

-   Added "TCellSI" T-cell reference database.
-   Introduced `Celltype_Calculate()` for automated scoring.
-   Added `Celltype_Annotation()` for end-to-end cell typing.
-   Improved message output system.
-   Resolved multiple code errors.
-   Modified code to meet CRAN standards and policies.

# SlimR 1.0.4 (2025-07-30)

-   Optimized `Celltype_annotation_Heatmap()` performance.
-   Enhanced probability calculation in `calculate_probability()`.
-   Modified code to meet CRAN submission requirements.
-   Change the License type from "GPL-3" to "MIT".

# SlimR 1.0.3 (2025-07-28)

-   Replaced `calculate_mean_expression()` with `calculate_probability()` in `Celltype_annotation_Heatmap()`.
-   Modified code to meet CRAN standards and policies.

# SlimR 1.0.1 (2025-07-19)

-   Changed `Celltype_annotation_Bar()` to `Celltype_annotation_Box()` with improved visualization capabilities.

# SlimR 1.0.0 (2025-07-07)

-   Initial release of SlimR package with core cell type annotation framework and basic visualization functions.