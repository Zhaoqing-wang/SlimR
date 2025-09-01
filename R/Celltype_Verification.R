#' Perform cell type verification and generate the validation dotplot
#'
#' @description This function performs verification of predicted cell types by selecting
#'     high log2FC and high expression proportion genes and generates and generate the
#'     validation dotplot.
#'
#' @param seurat_obj A Seurat object containing single-cell data.
#' @param SlimR_anno_result A list containing SlimR annotation results with:
#'     Expression_list - List of expression matrices for each cell type.
#'     Prediction_results - Data frame with cluster annotations.
#' @param assay Enter the assay used by the Seurat object, such as "RNA". Default
#'     parameters use "assay = 'RNA'".
#' @param gene_number Integer specifying number of top genes to select per cell type.
#' @param annotation_col Character string specifying the column in meta.data to use for grouping.
#' @param colour_low Color for lowest expression level. (default = "white")
#' @param colour_high Color for highest expression level. (default = "navy")
#'
#' @return A ggplot object showing expression of top variable genes.
#'
#' @export
#' @family Section_3_Automated_Annotation_Workflow
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
#'     assay = "RNA",
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
    assay = "RNA",
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

  Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[annotation_col]]
  assay <- if (is.null(assay)) DefaultAssay(seurat_obj) else assay

  predicted_types <- unique(names(table(seurat_obj@meta.data[[annotation_col]])))
  predicted_types <- predicted_types[!is.na(predicted_types)]

  feature_list <- list()
  cell_types <- predicted_types
  total <- length(cell_types)
  cycles <- 0

  message(paste0("SlimR verification: The input idents: '",annotation_col,"' has ",total," cell types to be verify."))

  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    message(paste0("\n","[", i, "/", total, "] Verifying cell type: ", cell_type))

    if (cell_type %in% names(SlimR_anno_result$Expression_list)) {
      message(paste0("[", i, "/", total, "] ",cell_type," was verified using the markers information from the 'SlimR_anno_result$Expression_list'."))

      markers <- unique(colnames(SlimR_anno_result$Expression_list[[cell_type]]))

      prob_expression <- calculate_probability(object = seurat_obj,
                                               cluster_col = annotation_col,
                                               assay = assay,
                                               features = markers)

      if (is.null(prob_expression) || is.null(prob_expression$cluster_expr)) {
        next
      }

      if (nrow(prob_expression$cluster_expr) == 0) {
        next
      }

      expr_df <- prob_expression$cluster_expr
      frac_df <- prob_expression$cluster_frac

      current_expr <- as.numeric(expr_df[cell_type, ])
      names(current_expr) <- colnames(expr_df)

      other_means <- colMeans(expr_df[setdiff(cell_types, cell_type), , drop = FALSE])
      other_means[other_means == 0] <- 1e-5
      log2fc <- log2(current_expr / other_means)

      current_frac <- as.numeric(frac_df[cell_type, ])
      scores <- log2fc * current_frac
      scores[is.nan(scores)] <- 0

      sorted_genes <- names(sort(scores, decreasing = TRUE))
      top_genes <- head(sorted_genes, gene_number)
      feature_list[[cell_type]] <- top_genes

      cycles <- cycles + 1

    } else if (cell_type %in% predicted_types) {
      message(paste0("[", i, "/", total, "] ",cell_type," not found in 'SlimR_anno_result$Expression_list', using marker information from the function 'FindMarkers()' to verify."))

      Seurat::DefaultAssay(seurat_obj) <- assay
      markers <- Seurat::FindMarkers(seurat_obj, ident.1 = cell_type, only.pos = TRUE)

      markers$gene <- row.names(markers)
      markers$score <- markers$avg_log2FC * markers$pct.1

      if (nrow(markers) == 0) next

      top_markers <- markers %>%
        top_n(gene_number, score) %>%
        dplyr::pull(gene)

      feature_list[[cell_type]] <- top_markers

      cycles <- cycles + 1

    } else {
      next
    }
  }

  message(paste0(" "))
  
  all_features <- unique(unlist(feature_list))
  if (length(all_features) == 0) stop("No valid features found for verification")

  cluster_order <- unique(seurat_obj@meta.data[[annotation_col]])
  cluster_order <- cluster_order[!is.na(cluster_order)]

  dotplot <- Seurat::DotPlot(
    object = seurat_obj,
    features = all_features,
    assay = assay,
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

  message(paste0("\n","SlimR verification: Use the cell group identity information in 'seurat_obj@meta.data$",annotation_col," and use the product value of 'log2FC' and 'expression ratio' as the ranking basis."))

  return(dotplot)
}
