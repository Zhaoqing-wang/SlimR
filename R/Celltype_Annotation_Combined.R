#' Uses "marker_list" to generate combined plot for cell annotation
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
#' @param save_path The output path of the cell annotation picture. Example parameters
#'     use "save_path = './SlimR/Celltype_annotation_Bar/'".
#' @param colour_low Color for lowest expression level. (default = "white")
#' @param colour_high Color for highest expression level. (default = "black")
#'
#' @returns The cell annotation picture is saved in "save_path".
#' @export
#' @family Semi_Automated_Annotation_Workflow
#'
#' @importFrom Seurat Idents
#' @importFrom ggplot2 ggsave
#' @importFrom dplyr distinct %>%
#' @importFrom ggplot2 geom_boxplot geom_point position_dodge scale_color_gradient
#'
#' @examples
#' \dontrun{
#' Celltype_Annotation_Combined(seurat_obj = sce,
#'     gene_list = Markers_list,
#'     species = "Human",
#'     cluster_col = "seurat_clusters",
#'     assay = "RNA",
#'     save_path = file.path(tempdir(),"SlimR_Celltype_Annotation_Combined"),
#'     colour_low = "white",
#'     colour_high = "navy"
#'     )
#'     }
#'
Celltype_Annotation_Combined <- function(
    seurat_obj,
    gene_list,
    species,
    cluster_col = "seurat_clusters",
    assay = "RNA",
    save_path = NULL,
    colour_low = "white",
    colour_high = "navy"
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

  cell_types <- names(gene_list)
  total <- length(cell_types)
  cycles <- 0

  message(paste0("SlimR: The input 'Markers_list' has ",total," cell types to be processed."))

  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    message(paste0("\n","[", i, "/", total, "] Processing cell type: ", cell_type))

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
    plot_height <- 7
    plot_width <- max(8, num_clusters * 0.8) + 2

    split_and_format_features <- function(features, n_per_line = 10) {
      if (length(features) == 0) return("")
      feature_groups <- split(features, ceiling(seq_along(features) / n_per_line))
      formatted_lines <- lapply(feature_groups, function(group) {
        paste(group, collapse = ", ")
      })
      paste(formatted_lines, collapse = "\n")
    }

    expr_df <- calculate_expression(
      object = seurat_obj,
      features = gene_order_processed,
      assay = assay,
      cluster_col = cluster_col,
      colour_low = colour_low,
      colour_high = colour_high
    )

    p <- ggplot(expr_df, aes(x = cluster, y = mean_expression, color = Avg_exp, fill = Avg_exp)) +
      geom_boxplot(outlier.shape = NA, width = 0.6) +
      geom_point(position = position_dodge(width = 0.6), size = 2, stroke = 0) +
      scale_color_gradient(low = colour_low, high = colour_high) +
      scale_fill_gradient(low = colour_low, high = colour_high) +
      labs(
        title = paste0(cell_type," combined features expression per Cluster | SlimR"),
        subtitle = paste("Features:", split_and_format_features(gene_order_processed, n_per_line = 15)),
        x = "Cell Cluster",
        y = "Mean Expression"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        panel.grid = element_blank()
      )

    ggsave(
      filename = file.path(save_path, paste0(cell_type, " mean expression.png")),
      plot = p,
      height = plot_height,
      width = plot_width,
      limitsize = FALSE
    )
    cycles <- cycles + 1
    message(paste0("[", i, "/", total, "] Combined plot saved for: ", cell_type))
  }

  message(paste0("\n","SlimR: Out of the ",total," cell types in 'Markers_list', ",cycles," cell types have been processed. You can see the reason for not processing cell types by 'warnings()'."))

  message(paste0("\n","SlimR: Visualization saved to: ", normalizePath(save_path)))
}
