SlimR NEWS

== SlimR 1.0.7 (2025-08-11) ==

* MAJOR UPDATES:
  - Added new function `Celltype_Verification()` for predicted cell type validation
    using high-variability genes selection and dotplot visualization
  - Renamed `Celltype_annotation_Dotplot()` to `Celltype_annotation_Combined()`
    for unified function naming structure

* IMPROVEMENTS:
  - Enhanced README with detailed process descriptions
  - Optimized message output system for cleaner console feedback
  - Resolved various code bugs reported by users

* CRAN COMPLIANCE:
  - Modified codebase to meet CRAN standards and policies


== SlimR 1.0.6 (2025-08-06) ==

* NEW FEATURES:
  - Integrated "scIBD" human intestine reference database
  - Added AUC calculation and visualization to `Celltype_Calculate()`
  - Implemented AUC-based prediction correction in cell typing

* MAINTENANCE:
  - Streamlined code output formatting
  - Fixed critical bugs in prediction pipeline
  - CRAN standards compliance modifications


== SlimR 1.0.5 (2025-08-05) ==

* DATABASE EXPANSION:
  - Added "TCellSI" T-cell reference database

* NEW FUNCTIONALITY:
  - Introduced `Celltype_Calculate()` for automated scoring
  - Added `Celltype_Annotation()` for end-to-end cell typing

* OPTIMIZATIONS:
  - Improved message output system
  - Resolved multiple code errors
  - CRAN standards compliance enhancements


== SlimR 1.0.4 (2025-07-30) ==

* ALGORITHM IMPROVEMENTS:
  - Optimized `Celltype_annotation_Heatmap()` performance
  - Enhanced probability calculation in `calculate_probability()`

* CRAN COMPLIANCE:
  - Modified code to meet CRAN submission requirements


== SlimR 1.0.3 (2025-07-28) ==

* FUNCTIONAL UPGRADE:
  - Replaced `calculate_mean_expression()` with 
    `calculate_probability()` in `Celltype_annotation_Heatmap()`

* CRAN RELEASE:
  - Initial CRAN-compliant version with code modifications


== SlimR 1.0.1 (2025-07-19) ==

* FUNCTION RENAMING:
  - Changed `Celltype_annotation_Bar()` to `Celltype_annotation_Box()`
    with improved visualization capabilities


== SlimR 1.0.0 (2025-07-07) ==

* INITIAL RELEASE:
  - First stable version of SlimR package
  - Core cell type annotation framework
  - Basic visualization functions
