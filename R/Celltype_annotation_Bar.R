#' Title
#'
#' @param seurat_obj
#' @param gene_list
#' @param species
#' @param cluster_col
#' @param assay
#' @param save_path
#' @param metric_names
#'
#' @returns
#' @export
#'
#' @examples
Celltype_annotation_Bar <- function(
    seurat_obj,
    gene_list,
    species,
    cluster_col = "seurat_clusters",
    assay = "RNA",
    save_path = "./Celltype_annotation_Excel/",
    metric_names = NULL
) {
  required_packages <- c("ggplot2", "patchwork", "dplyr", "scales", "tidyr", "gridExtra", "gtable", "grid")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
  }

  if (!inherits(seurat_obj, "Seurat")) stop("Input object must be a Seurat object!")
  if (!is.list(gene_list)) stop("Gene list must be a list of data.frames!")
  if (species != "Human" && species != "Mouse") stop("species must be 'Human' or 'Mouse'")

  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

  for (cell_type in names(gene_list)) {
    cat("Processing cell type:", cell_type, "\n")
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

    num_clusters <- length(unique(Seurat::Idents(seurat_obj)))
    num_genes <- length(gene_order_original)
    plot_height <- max(6, num_clusters * 0.5) + 2
    plot_width <- 5

    bar_plot <- plot_mean_expression(
      object = seurat_obj,
      features = gene_order_processed,
      assay = assay,
      cluster_col = cluster_col
    )

    total_width <- plot_width
    ggsave(
      filename = file.path(save_path, paste0(cell_type, " mean expression.png")),
      plot = bar_plot,
      height = plot_height,
      width = plot_width,
      limitsize = FALSE
    )
    cat("Bar plot saved for", cell_type, "\n\n")
  }

  message("Visualization saved to:", normalizePath(save_path))
}
