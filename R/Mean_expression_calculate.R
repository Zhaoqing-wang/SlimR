#' Title
#'
#' @param object
#' @param features
#' @param assay
#' @param cluster_col
#'
#' @returns
#'
#' @examples
calculate_mean_expression <- function(object, features,assay = NULL, cluster_col = NULL) {
  assay <- if (is.null(assay)) DefaultAssay(object) else assay
  DefaultAssay(object) <- assay

  cells <- unlist(CellsByIdentities(object = object))

  data.features <- FetchData(
    object = object,
    vars = features,
    cells = cells
  )

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

  return(rowMeans(expr_matrix))
}
