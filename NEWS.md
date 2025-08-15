# SlimR 1.0.7 (2025-08-16)

-   Added new function `Celltype_Verification()` for predicted cell types validation and generate the validation dotplot.
-   Optimize the function `read_seurat_markers()`. This can be compatible with the `presto::wilcoxauc()` source markers.
-   Renamed `Celltype_annotation_Dotplot()` to `Celltype_Annotation_Features()` for unified function naming structure.
-   Renamed `Celltype_annotation_Box()` to `Celltype_Annotation_Combined()` for unified function naming structure.
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