#' Uses "marker_list" to generate Box plot for cell annotation
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
#'     "seurat_clusters"".
#' @param assay Enter the assay used by the Seurat object, such as "RNA". Default
#'     parameters use "assay = "RNA"".
#' @param save_path The output path of the cell annotation picture. Default parameters
#'     use "save_path = "./SlimR/Celltype_annotation_Bar/"".
#' @param metric_names Warning: Do not enter information. This parameter is used to
#'     check if "Marker_list" conforms to the SlimR output.
#'
#' @returns The cell annotation picture is saved in "save_path".
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom Seurat Idents
#' @importFrom ggplot2 ggsave
#' @importFrom dplyr distinct
#'
#' @examples
#' \dontrun{Celltype_annotation_Box(seurat_obj = sce,
#'          gene_list = Markers_list,
#'          species = "Human",
#'          cluster_col = "seurat_clusters",
#'          assay = "RNA",
#'          save_path = file.path(tempdir(),"SlimR_Celltype_annotation_Box")
#'          )
#'          }
#'
Celltype_annotation_Box <- function(
    seurat_obj,
    gene_list,
    species,
    cluster_col = "seurat_clusters",
    assay = "RNA",
    save_path = NULL,
    metric_names = NULL
) {
  if (!inherits(seurat_obj, "Seurat")) stop("Input object must be a Seurat object!")
  if (!is.list(gene_list)) stop("Gene list must be a list of data.frames!")
  if (species != "Human" && species != "Mouse") stop("species must be 'Human' or 'Mouse'")
  if (missing(save_path)) {stop("Output path must be explicitly specified")}
  if (!interactive() && !grepl(tempdir(),save_path, fixed = TRUE)) {
    warning("Writing to non-temporary locations is restricted", immediate. = TRUE)
    path <- file.path(tempdir(), "fallback_output")
  }

  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

  for (cell_type in names(gene_list)) {
    message("Processing cell type:", cell_type, "\n")
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
    message("Bar plot saved for", cell_type, "\n\n")
  }

  message("Visualization saved to:", normalizePath(save_path))
}
