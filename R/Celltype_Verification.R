#' Perform cell type verification using top variable genes
#'
#' @description This function performs verification of predicted cell types by selecting
#'     high-variability genes and generating a dotplot visualization.
#'
#' @param seurat_obj A Seurat object containing single-cell data.
#' @param SlimR_anno_result A list containing SlimR annotation results with:
#'     Expression_list - List of expression matrices for each cell type.
#'     Prediction_results - Data frame with cluster annotations.
#' @param gene_number Integer specifying number of top genes to select per cell type.
#' @param annotation_col Character string specifying the column in meta.data to use for grouping.
#' @param colour_low Color for lowest expression level. (default = "white")
#' @param colour_high Color for highest expression level. (default = "black")
#'
#' @return A ggplot object showing expression of top variable genes.
#'
#' @export
#' @family Celltype_annotation
#'
#' @importFrom Seurat DotPlot FetchData
#' @importFrom dplyr distinct bind_rows arrange desc top_n
#' @importFrom stats sd
#' @importFrom ggplot2 theme_bw element_blank element_text guide_legend scale_color_gradientn
#' @importFrom ggplot2 ggtitle
#'
#' @examples
#' \dontrun{
#' Celltype_Verification(seurat_obj = sce,
#'     SlimR_anno_result = SlimR_anno_result,
#'     gene_number = 5,
#'     colour_low = "white",
#'     colour_high = "navy",
#'     annotation_col = "Cell_type_SlimR"
#'     )
#'     }
#'
Celltype_Verification <- function(
    seurat_obj,
    SlimR_anno_result,
    gene_number = 5,
    colour_low = "white",
    colour_high = "navy",
    annotation_col = "Cell_type_SlimR"
) {
  if (!inherits(seurat_obj, "Seurat")) stop("seurat_obj must be a Seurat object")
  if (!is.list(SlimR_anno_result)) stop("SlimR_anno_result must be a list")
  if (!"Prediction_results" %in% names(SlimR_anno_result)) stop("Prediction_results not found in SlimR_anno_result")
  if (!"Expression_list" %in% names(SlimR_anno_result)) stop("Expression_list not found in SlimR_anno_result")
  if (!is.numeric(gene_number) || gene_number < 1) stop("gene_number must be a positive integer")
  if (!(annotation_col %in% colnames(seurat_obj@meta.data))) stop(paste0(annotation_col, " not found in seurat_obj meta.data, please run Celltype_Annotation() first."))

  predicted_types <- unique(SlimR_anno_result$Prediction_results$Predicted_cell_type)
  predicted_types <- predicted_types[!is.na(predicted_types)]
  cv <- NULL

  feature_list <- list()
  for (cell_type in predicted_types) {
    if (!(cell_type %in% names(SlimR_anno_result$Expression_list))) next

    expr_df <- SlimR_anno_result$Expression_list[[cell_type]]
    if (nrow(expr_df) == 0) next

    cv_values <- apply(expr_df, 2, function(x) {
      if (mean(x) == 0) return(0)
      sd(x) / mean(x)
    })

    cv_df <- data.frame(
      gene = names(cv_values),
      cv = cv_values
    ) %>%
      arrange(desc(cv)) %>%
      top_n(gene_number, cv)

    feature_list[[cell_type]] <- cv_df$gene
  }

  all_features <- unique(unlist(feature_list))
  if (length(all_features) == 0) stop("No valid features found for verification")

  cluster_order <- unique(seurat_obj@meta.data[[annotation_col]])
  cluster_order <- cluster_order[!is.na(cluster_order)]

  dotplot <- Seurat::DotPlot(
    object = seurat_obj,
    features = all_features,
    group.by = annotation_col
  ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = "Celltype verification dotplot | SlimR"
    ) +
    ggplot2::guides(size = ggplot2::guide_legend(order = 3)) +
    ggplot2::scale_color_gradientn(
      colours = c(colour_low, colour_high),
      values = seq(0, 1, length.out = 2))

  message(paste0("SlimR verification: By using the highly variable genes in 'SlimR_anno_result$Expression_list' corresponding to predicted cell type information in 'SlimR_anno_result$Prediction_results$Predicted_cell_type', the dotplot is generated based on clusters as 'seurat_obj@meta.data$",annotation_col,"."))
  return(dotplot)
}
