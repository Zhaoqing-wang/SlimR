#' "Use in package" Counts average expression of gene set and plots barplot
#'
#' @param object Enter a Seurat object.
#' @param features Enter one or a set of markers.
#' @param assay Enter the assay used by the Seurat object, such as "RNA". Default
#'     parameters use "assay = NULL".
#' @param cluster_col Enter the meta.data column in the Seurat object to be
#'     annotated, such as "seurat_cluster". Default parameters use "cluster_col = NULL".
#'
#' @returns  Average expression bar plot of genes in the input "Seurat" object
#'     given "cluster_col" and given "features".
#'
#' @examples
#' \donttest{plot_mean_expression(sce.all,
#'           features = c("CD19","CD79A","MS4A1")
#'           )
#'           }
#'
plot_mean_expression <- function(object, features, assay = NULL, cluster_col = NULL) {
  if (!is.null(cluster_col) && !(cluster_col %in% colnames(object@meta.data))) {
    stop("cluster_col not found in meta.data")
  }

  assay <- assay %||% DefaultAssay(object)
  DefaultAssay(object) <- assay

  cells <- unlist(CellsByIdentities(object = object,
                                    cells = colnames(object[[assay]])))

  data.features <- FetchData(object = object,
                             vars = features,
                             cells = cells)

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

  df <- data.frame(
    cluster = rownames(expr_matrix),
    Mean_exp = rowMeans(expr_matrix)
  )
  df$cluster <- factor(df$cluster, levels = id.levels)

  p <- ggplot(df, aes(x = Mean_exp, y = cluster, fill = Mean_exp)) +
    geom_col(width = 0.7) +
    scale_fill_gradient(low = "lightgrey", high = "dodgerblue") +
    labs(
      title = "Mean expression between clusters",
      subtitle = paste("Features:", paste(features, collapse = ", ")),
      x = "Mean Expression",
      y = "Cell Cluster"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10, hjust = 0.5, vjust = 0.5),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      panel.grid = element_blank()
    )

  return(p)
}
