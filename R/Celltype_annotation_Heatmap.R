#' Title
#'
#' @param seurat_obj
#' @param gene_list
#' @param species
#' @param cluster_col
#' @param assay
#'
#' @returns
#' @export
#'
#' @examples
Celltype_annotation_Heatmap <- function(
    seurat_obj,
    gene_list,
    species,
    cluster_col = "seurat_clusters",
    assay = "RNA"
) {
  required_packages <- c("ggplot2", "patchwork", "dplyr", "scales", "tidyr", "gridExtra", "gtable", "grid")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
  }

  if (!inherits(seurat_obj, "Seurat")) stop("Input object must be a Seurat object!")
  if (!is.list(gene_list)) stop("Gene list must be a list of data.frames!")
  if (species != "Human" && species != "Mouse") stop("species must be 'Human' or 'Mouse'")

  cluster_expr_list <- list()
  for (cell_type in names(gene_list)) {
    message("Processing cell type:", cell_type)
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

    mean_expression <- calculate_mean_expression(object = seurat_obj,
                                                 cluster_col = cluster_col,
                                                 assay = assay,
                                                 features = gene_order_processed)
    cluster_expr_list[[cell_type]] <- scales::rescale(mean_expression , na.rm = TRUE)
    message(paste0(cell_type," mean expression calculated","\n"))
  }

  expr_matrix <- do.call(rbind, cluster_expr_list)

  result_matrix <- t(expr_matrix)

  p <- pheatmap::pheatmap(result_matrix,
                main = "Cell annotation heatmap | SlimR",
                color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                fontsize = 12,
                cluster_rows = T,
                cluster_cols = T,
                legend_breaks = c(0,0.5,1),
                legend_labels = c("Low", "Median", "High"))
  return(p)
}
