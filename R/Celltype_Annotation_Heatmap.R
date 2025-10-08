#' Uses "marker_list" to generate heatmap for cell annotation
#'
#' @param seurat_obj Enter the Seurat object with annotation columns such as
#'     "seurat_cluster" in meta.data to be annotated.
#' @param gene_list A list of cells and corresponding gene controls, the name of
#'     the list is cell type, and the first column of the list corresponds to markers.
#'     Lists can be generated using functions such as "Markers_filter_Cellmarker2 ()",
#'     "Markers_filter_PanglaoDB ()", "read_excel_markers ()", "read_seurat_markers ()", etc.
#' @param species This parameter selects the species "Human" or "Mouse" for standard
#'     gene format correction of markers entered by "Marker_list".
#' @param cluster_col Enter annotation columns such as "seurat_cluster" in meta.data
#'     of the Seurat object to be annotated. Default parameters use "cluster_col =
#'     'seurat_clusters'".
#' @param assay Enter the assay used by the Seurat object, such as "RNA". Default
#'     parameters use "assay = 'RNA'".
#' @param min_expression The min_expression parameter defines a threshold value to
#'     determine whether a cell's expression of a feature is considered "expressed"
#'     or not. It is used to filter out low-expression cells that may contribute
#'     noise to the analysis. Default parameters use "min_expression = 0.1".
#' @param specificity_weight The specificity_weight parameter controls how much the
#'     expression variability (standard deviation) of a feature within a cluster
#'     contributes to its "specificity score." It amplifies or suppresses the impact
#'     of variability in the final score calculation.Default parameters use
#'     "specificity_weight = 3".
#' @param colour_low Color for lowest probability level in Heatmap visualization of
#'     probability matrix. (default = "navy")
#' @param colour_high Color for highest probability level Heatmap visualization of
#'     probability matrix. (default = "firebrick3")
#'
#' @returns The heatmap of the comparison between "cluster_col" in the
#'     Seurat object and the given gene set "gene_list" needs to be annotated.
#' @export
#' @family Section_4_Semi_Automated_Annotation
#'
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' \dontrun{
#' Celltype_Annotation_Heatmap(seurat_obj = sce,
#'     gene_list = Markers_list,
#'     species = "Human",
#'     cluster_col = "seurat_clusters",
#'     assay = "RNA",
#'     min_expression = 0.1,
#'     specificity_weight = 3,
#'     colour_low = "navy",
#'     colour_high = "firebrick3"
#'     )
#'     }
#'
Celltype_Annotation_Heatmap <- function(
    seurat_obj,
    gene_list,
    species,
    cluster_col = "seurat_clusters",
    assay = "RNA",
    min_expression = 0.1,
    specificity_weight = 3,
    colour_low = "navy",
    colour_high = "firebrick3"
) {
  required_packages <- c("ggplot2", "patchwork", "dplyr", "scales", "tidyr", "gridExtra", "gtable", "grid")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Please install the required package: %s", pkg))
    }
    library(pkg, character.only = TRUE)
  }

  if (!inherits(seurat_obj, "Seurat")) stop("Input object must be a Seurat object!")
  if (!is.list(gene_list)) stop("Gene list must be a list of data.frames!")
  if (species != "Human" && species != "Mouse") stop("species must be 'Human' or 'Mouse'")

  cluster_scores_list <- list()
  cell_types <- names(gene_list)

  message(paste0("SlimR: The 'Celltype_annotation_Heatmap()' function has now been replaced by the 'Celltype_Calculate()' function. You can still use it, but this function is no longer actively updated. It is recommended to use 'Celltype_Calculate()' instead."))

  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]

    current_df <- gene_list[[cell_type]]

    if (ncol(current_df) < 1) {
      warning(paste("Skipping", cell_type, ": Requires at least a gene column"))
      next
    }

    genes <- current_df[[1]]
    genes_processed <- if (species == "Human") {
      toupper(genes)
    } else {
      paste0(toupper(substr(genes, 1, 1)), tolower(substr(genes, 2, nchar(genes))))
    }

    valid_idx <- genes_processed %in% rownames(seurat_obj[[assay]])
    if (sum(valid_idx) == 0) {
      warning(paste("No valid genes for", cell_type))
      next
    }

    valid_data <- data.frame(
      original = genes[valid_idx],
      processed = genes_processed[valid_idx],
      stringsAsFactors = FALSE
    ) %>% distinct(processed, .keep_all = TRUE)

    gene_order_processed <- valid_data$processed
    gene_order_original <- valid_data$original

    prob_expression <- calculate_probability(object = seurat_obj,
                                             cluster_col = cluster_col,
                                             assay = assay,
                                             features = gene_order_processed,
                                             min_expression = min_expression,
                                             specificity_weight = specificity_weight)
    cluster_scores_list[[cell_type]] <- prob_expression$cluster_scores
  }

  expr_matrix <- do.call(rbind, cluster_scores_list)

  normalize_row <- function(x) {
    if (diff(range(x)) == 0) return(rep(0, length(x)))
    (x - min(x)) / (max(x) - min(x))
  }

  normalize_matrix <- apply(expr_matrix, 2, normalize_row)

  result_matrix <- t(normalize_matrix)

  p <- pheatmap::pheatmap(result_matrix,
                          main = "Cell annotation heatmap | SlimR",
                          color = colorRampPalette(c(colour_low, "white", colour_high))(50),
                          fontsize = 12,
                          cluster_rows = T,
                          cluster_cols = T,
                          legend_breaks = c(0,1),
                          legend_labels = c("Low probability","High probability"))
  return(p)
}
