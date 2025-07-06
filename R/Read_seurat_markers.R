#' Title
#'
#' @param df
#' @param sort_by
#' @param gene_fliter
#'
#' @returns
#' @export
#'
#' @examples
read_seurat_markers <- function(df, sort_by = "avg_log2FC", gene_fliter = 20) {
  if (!sort_by %in% c("avg_log2FC", "p_val_adj")) {
    stop("sort_by must be either 'avg_log2FC' or 'p_val_adj'")
  }

  clusters <- split(df, df$cluster)

  processed <- lapply(clusters, function(cluster_df) {
    if (sort_by == "avg_log2FC") {
      sorted_df <- cluster_df[order(-cluster_df$avg_log2FC), , drop = FALSE]
    } else {
      sorted_df <- cluster_df[order(cluster_df$p_val_adj), , drop = FALSE]
    }

    filtered_df <- head(sorted_df, gene_fliter)
    desired_order <- c("gene", "avg_log2FC", "p_val_adj", "p_val", "pct.1", "pct.2")
    reordered_df <- filtered_df[, desired_order, drop = FALSE]

    return(reordered_df)
  })

  names(processed) <- names(clusters)

  return(processed)
}
