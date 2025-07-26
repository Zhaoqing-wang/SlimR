#' Counts average expression of gene set and plots Boxplot (Use in package)
#'
#' @param object Enter a Seurat object.
#' @param features Enter one or a set of markers.
#' @param assay Enter the assay used by the Seurat object, such as "RNA". Default
#'     parameters use "assay = NULL".
#' @param cluster_col Enter the meta.data column in the Seurat object to be
#'     annotated, such as "seurat_cluster". Default parameters use "cluster_col = NULL".
#'
#' @returns  Average expression box plot of genes in the input "Seurat" object
#'     given "cluster_col" and given "features".
#'
#' @importFrom Seurat `%||%`
#' @importFrom Seurat DefaultAssay DefaultAssay<- CellsByIdentities FetchData
#' @importFrom dplyr group_by summarise left_join
#' @importFrom ggplot2 geom_boxplot geom_point position_dodge scale_color_gradient
#'
#'
plot_mean_expression <- function(object, features, assay = NULL, cluster_col = NULL) {
  if (!is.null(cluster_col) && !(cluster_col %in% colnames(object@meta.data))) {
    stop("cluster_col not found in meta.data")
  }

  assay <- assay %||% DefaultAssay(object)
  DefaultAssay(object) <- assay

  cells <- unlist(CellsByIdentities(object = object, cells = colnames(object[[assay]])))

  data.features <- FetchData(object = object, vars = features, cells = cells)

  if (!is.null(cluster_col)) {
    data.features$id <- object@meta.data[cells, cluster_col, drop = TRUE]
  } else {
    data.features$id <- Idents(object = object)[cells]
  }

  id.levels <- levels(factor(data.features$id))
  data.features$id <- factor(data.features$id, levels = id.levels)

  cluster_expr_list <- lapply(unique(data.features$id), function(ident) {
    data.use <- data.features[data.features$id == ident, features, drop = FALSE]
    avg.exp <- apply(data.use, 2, function(x) mean(expm1(x)))
    return(avg.exp)
  })

  expr_matrix <- do.call(rbind, cluster_expr_list)
  rownames(expr_matrix) <- unique(data.features$id)

  expr_df <- as.data.frame(expr_matrix) %>%
    tibble::rownames_to_column("cluster") %>%
    tidyr::pivot_longer(
      cols = -cluster,
      names_to = "gene",
      values_to = "mean_expression"
    )

  expr_df$cluster <- factor(expr_df$cluster, levels = id.levels)

  cluster_avg_exp <- expr_df %>%
    group_by(cluster) %>%
    summarise(avg_cluster_exp = mean(mean_expression), .groups = "drop")

  expr_df <- left_join(expr_df, cluster_avg_exp, by = "cluster")

  p <- ggplot(expr_df, aes(x = cluster, y = mean_expression, color = avg_cluster_exp, fill = avg_cluster_exp)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_point(position = position_dodge(width = 0.6), size = 2, stroke = 0) +
    scale_color_gradient(low = "lightgrey", high = "dodgerblue") +
    scale_fill_gradient(low = "lightgrey", high = "dodgerblue") +
    labs(
      title = "Gene Expression per Cluster",
      subtitle = paste("Features:", paste(features, collapse = ", ")),
      x = "Cell Cluster",
      y = "Mean Expression"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      panel.grid = element_blank()
    )

  return(p)
}
