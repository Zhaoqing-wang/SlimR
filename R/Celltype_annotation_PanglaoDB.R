#' Uses "marker_list" from PanglaoDB for cell annotation
#'
#' @param seurat_obj Enter the Seurat object with annotation columns such as
#'     "seurat_cluster" in meta.data to be annotated.
#' @param gene_list Enter the standard "Marker_list" generated by the PanglaoDB
#'     database for the SlimR package, generated by the "Markers_filter_PanglaoDB ()"
#'     function.
#' @param species This parameter selects the species "Human" or "Mouse" for standard
#'     gene format correction of markers entered by "Marker_list".
#' @param cluster_col Enter annotation columns such as "seurat_cluster" in meta.data
#'     of the Seurat object to be annotated. Default parameters use "cluster_col =
#'     "seurat_clusters"".
#' @param assay Enter the assay used by the Seurat object, such as "RNA". Default
#'     parameters use "assay = "RNA"".
#' @param save_path The output path of the cell annotation picture. Default parameters
#'     use "save_path = "./SlimR/Celltype_annotation_PanglaoDB/"".
#' @param metric_names Warning: Do not enter information. This parameter is used to
#'     check if "Marker_list" conforms to the PanglaoDB database output.
#'
#' @returns The cell annotation picture is saved in "save_path".
#' @export
#'
#' @importFrom dplyr all_of
#'
#' @examples
#' \dontrun{
#' Celltype_annotation_PanglaoDB(seurat_obj = sce,
#'     gene_list = Markers_list_panglaoDB,
#'     species = "Human",
#'     cluster_col = "seurat_clusters",
#'     assay = "RNA",
#'     save_path = file.path(tempdir(),"SlimR_Celltype_annotation_PanglaoDB")
#'     )
#'     }
#'
Celltype_annotation_PanglaoDB <- function(
    seurat_obj,
    gene_list,
    species,
    cluster_col = "seurat_clusters",
    assay = "RNA",
    save_path = NULL,
    metric_names = NULL
) {
  required_packages <- c("ggplot2", "patchwork", "dplyr", "scales")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Please install the required package: %s", pkg))
    }
    library(pkg, character.only = TRUE)
  }

  if (!inherits(seurat_obj, "Seurat")) stop("Input object must be a Seurat object!")
  if (!is.list(gene_list)) stop("Gene list must be a list of data.frames!")
  if (species != "Human" && species != "Mouse") stop("species must be 'Human' or 'Mouse'")
  if (missing(save_path)) {stop("Output path must be explicitly specified")}
  if (!interactive() && !grepl(tempdir(), save_path, fixed = TRUE)) {
    warning("Writing to non-temporary locations is restricted", immediate. = TRUE)
    path <- file.path(tempdir(), "fallback_output")
  }

  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

  common_theme <- function(base_size = 10) {
    ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, family = "sans"),
        axis.title = ggplot2::element_text(family = "sans"),
        plot.title = ggplot2::element_text(hjust = 0, face = "bold", size = 12),
        legend.position = "right",
        panel.grid = ggplot2::element_blank()
      )
  }

  for (cell_type in names(gene_list)) {
    message("Processing cell type:", cell_type, "\n")
    current_df <- gene_list[[cell_type]]

    if (ncol(current_df) < 4) {
      warning(paste("Skipping", cell_type, ": Insufficient columns"))
      next
    }

    metric_cols <- if (!is.null(metric_names)) {
      if (length(metric_names) != 3) stop("Need exactly 3 metric names")
      metric_names
    } else {
      colnames(current_df)[2:4]
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

    filtered_data <- data.frame(
      gene = genes[valid_idx],
      processed = genes_processed[valid_idx],
      current_df[valid_idx, 2:4]
    )
    colnames(filtered_data)[3:5] <- metric_cols

    gene_order <- unique(filtered_data$gene)

    num_clusters <- length(unique(Seurat::Idents(seurat_obj)))
    num_genes <- length(gene_order)
    plot_height <- max(6, num_clusters * 0.5) + 2
    plot_width <- max(10, num_genes * 0.4)

    dp <- Seurat::DotPlot(
      seurat_obj,
      features = gene_order,
      assay = assay,
      group.by = cluster_col,
      cols = c("white", "dodgerblue")
    ) +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          family = "sans",
          size = 10
        ),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family = "sans")
      ) +
      labs(
        title = paste("Cell Type:", cell_type, "| Database: PanglaoDB | SlimR"),
        subtitle = "Dot size: Expression percentage | Color: Normalized expression level"
      )

    metric_data <- filtered_data %>%
      tidyr::pivot_longer(
        cols = all_of(metric_cols),
        names_to = "metric",
        values_to = "score"
      ) %>%
      dplyr::group_by(metric) %>%
      dplyr::mutate(
        scaled_score = scales::rescale(score, na.rm = TRUE)
      ) %>%
      dplyr::ungroup()

    hp <- ggplot2::ggplot(
      metric_data,
      ggplot2::aes(
        x = factor(gene, levels = gene_order),
        y = metric,
        fill = scaled_score
      )
    ) +ggplot2::geom_tile(color = NA) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradientn(
        colors = c("white", "dodgerblue"),
        na.value = "white",
        limits = c(0, 1)
      ) +
      ggplot2::labs(
        title = "Normalized metrics in database PanglaoDB",
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1
        ),
        axis.title.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank()
      )

    base_height <- 5
    heatmap_height_ratio <- 0.3

    combined_plot <- cowplot::plot_grid(
      dp,
      hp,
      ncol = 1,
      align = "v",
      rel_heights = c(1, heatmap_height_ratio)
    )

    ggsave(
      filename = file.path(save_path, paste0(cell_type, ".png")),
      plot = combined_plot,
      height = plot_height * (1 + heatmap_height_ratio),
      width = plot_width,
      limitsize = FALSE
    )
    message("Combined plot saved for", cell_type, "\n\n")
  }

  message("Visualization saved to:", normalizePath(save_path))
}
