#' Uses "marker_list" to calculate probability, prediction results, AUC and generate heatmap for cell annotation
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
#' @param min_expression The min_expression parameter defines a threshold value to
#'     determine whether a cell's expression of a feature is considered "expressed"
#'     or not. It is used to filter out low-expression cells that may contribute
#'     noise to the analysis. Default parameters use "min_expression = 0.1".
#' @param specificity_weight The specificity_weight parameter controls how much the
#'     expression variability (standard deviation) of a feature within a cluster
#'     contributes to its "specificity score." It amplifies or suppresses the impact
#'     of variability in the final score calculation.Default parameters use
#'     "specificity_weight = 3".
#' @param threshold This parameter refers to the normalized similarity between the
#'     "alternative cell type" and the "predicted cell type" in the returned results
#'     (the default parameter is 0.8).
#' @param compute_AUC Logical indicating whether to calculate AUC values for predicted
#'     cell types. AUC measures how well the marker genes distinguish the cluster from
#'     others. When TRUE, adds an AUC column to the prediction results. (default: TRUE)
#' @param AUC_correction Logical value controlling AUC-based correction (default = TRUE).
#'     When set to TRUE:
#'     - Computes AUC values for candidate cell types (probability > threshold)
#'     - Selects the cell type with the highest AUC as the final predicted type
#'     - Records the selected type's AUC value in the "AUC" column.
#'
#' @returns A list containing:
#' \itemize{
#'   \item Expression_list: List of expression matrices for each cell type
#'   \item Expression_scores_matrix: Matrix of expression scores
#'   \item Probability_matrix: Matrix of normalized probabilities
#'   \item Prediction_results: Data frame with cluster annotations including:
#'     \itemize{
#'       \item cluster_col: Cluster identifier
#'       \item Predicted_cell_type: Primary predicted cell type
#'       \item AUC: Area Under the Curve value (when compute_AUC = TRUE)
#'       \item Alternative_cell_types: Semi-colon separated alternative cell types
#'     }
#'   \item Heatmap_plot: Heatmap visualization of probability matrix
#' }
#'
#' @export
#' @family Celltype_annotation
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom utils tail
#'
#' @examples
#' \dontrun{
#' SlimR_anno_result <- Celltype_Calculate(seurat_obj = sce,
#'     gene_list = Markers_list,
#'     species = "Human",
#'     cluster_col = "seurat_clusters",
#'     assay = "RNA",
#'     min_expression = 0.1,
#'     specificity_weight = 3,
#'     compute_AUC = TRUE,
#'     AUC_correction = TRUE
#'     )
#'     }
#'
Celltype_Calculate <- function(
    seurat_obj,
    gene_list,
    species,
    cluster_col = "seurat_clusters",
    assay = "RNA",
    min_expression = 0.1,
    specificity_weight = 3,
    threshold = 0.8,
    compute_AUC = TRUE,
    AUC_correction = TRUE
) {
  required_packages <- c("ggplot2", "patchwork", "dplyr", "scales", "tidyr", "gridExtra", "gtable", "grid", "pheatmap")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Please install the required package: %s", pkg))
    }
    library(pkg, character.only = TRUE)
  }

  if (AUC_correction) compute_AUC <- TRUE

  if (!inherits(seurat_obj, "Seurat")) stop("Input object must be a Seurat object!")
  if (!is.list(gene_list)) stop("Gene list must be a list of data.frames!")
  if (species != "Human" && species != "Mouse") stop("species must be 'Human' or 'Mouse'")

  cluster_scores_list <- list()
  cluster_mean_list <- list()
  valid_genes_list <- list()

  cell_types <- names(gene_list)
  total <- length(cell_types)

  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    message(paste0("[", i, "/", total, "] Processing cell type: ", cell_type))

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

    valid_genes_list[[cell_type]] <- gene_order_processed

    prob_expression <- calculate_probability(object = seurat_obj,
                                             cluster_col = cluster_col,
                                             assay = assay,
                                             features = gene_order_processed,
                                             min_expression = min_expression,
                                             specificity_weight = specificity_weight)
    cluster_scores_list[[cell_type]] <- prob_expression$cluster_scores
    cluster_mean_list[[cell_type]] <- prob_expression$cluster_expr
    message(paste0("[", i, "/", total, "] ", cell_type)," characteristic genes expression calculated. \n")
  }

  expr_list <- cluster_mean_list
  scores_matrix <- do.call(rbind, cluster_scores_list)

  normalize_row <- function(x) {
    if (diff(range(x)) == 0) return(rep(0, length(x)))
    (x - min(x)) / (max(x) - min(x))
  }

  normalize_matrix <- apply(scores_matrix, 2, normalize_row)
  result_matrix <- t(normalize_matrix)

  p <- pheatmap::pheatmap(result_matrix,
                          main = "Cell annotation heatmap | SlimR",
                          color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                          fontsize = 12,
                          cluster_rows = T,
                          cluster_cols = T,
                          legend_breaks = c(0,1),
                          legend_labels = c("Low probability","High probability"))

  generate_prediction_table <- function(df, threshold = threshold) {
    clusters <- rownames(df)
    predicted_cell_types <- vector("character", length = length(clusters))
    alternative_cell_types <- vector("character", length = length(clusters))
    candidate_types_list <- list()

    for (i in seq_along(clusters)) {
      cluster <- clusters[i]
      row_values <- as.numeric(unlist(df[i, ]))
      cell_types <- names(df[i, ])
      max_index <- which.max(row_values)
      predicted <- cell_types[max_index]
      candidate_types <- cell_types[row_values > threshold]
      candidate_types_list[[cluster]] <- candidate_types

      alt <- candidate_types[candidate_types != predicted]
      alt_str <- if (length(alt) > 0) paste(alt, collapse = "; ") else NA_character_
      predicted_cell_types[i] <- predicted
      alternative_cell_types[i] <- alt_str
    }
    result_df <- data.frame(
      cluster_col = clusters,
      Predicted_cell_type = predicted_cell_types,
      Alternative_cell_types = alternative_cell_types,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    attr(result_df, "candidate_types") <- candidate_types_list
    return(result_df)
  }

  scores_matrix <- as.data.frame(t(scores_matrix))
  probability_matrix <- as.data.frame(result_matrix)
  prediction_results <- generate_prediction_table(probability_matrix, threshold = threshold)
  candidate_types_list <- attr(prediction_results, "candidate_types")

  fastAUC <- function(predictions, labels) {
    ord <- order(predictions, decreasing = TRUE)
    labels <- labels[ord]
    predictions <- predictions[ord]

    tpr <- cumsum(labels) / sum(labels)
    fpr <- cumsum(!labels) / sum(!labels)

    tpr <- c(0, tpr, 1)
    fpr <- c(0, fpr, 1)

    auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
    return(auc)
  }

  if (compute_AUC) {
    if (AUC_correction) {
      message(paste0("Performing AUC correction for all candidate cell types (threshold > ",threshold,")."))
      new_predicted <- character(nrow(prediction_results))
      new_aucs <- numeric(nrow(prediction_results))
      new_alt_list <- character(nrow(prediction_results))

      for (i in seq_len(nrow(prediction_results))) {
        cluster_id <- prediction_results$cluster_col[i]
        candidate_types <- candidate_types_list[[cluster_id]]

        if (length(candidate_types) == 0) {
          new_predicted[i] <- NA
          new_aucs[i] <- NA
          new_alt_list[i] <- NA
          next
        }

        auc_vals <- numeric(length(candidate_types))
        names(auc_vals) <- candidate_types

        for (j in seq_along(candidate_types)) {
          cell_type <- candidate_types[j]
          features <- valid_genes_list[[cell_type]]

          all_cells <- colnames(seurat_obj)
          expr_data <- FetchData(seurat_obj, vars = features, cells = all_cells)
          cell_scores <- rowMeans(expr_data, na.rm = TRUE)

          labels <- seurat_obj@meta.data[all_cells, cluster_col] == cluster_id
          if (length(unique(labels)) < 2) {
            auc_vals[j] <- NA
            warning(paste("Skipping AUC for cluster", cluster_id, "and cell type", cell_type, ": Only one class present"))
          } else {
            auc_vals[j] <- fastAUC(cell_scores, labels)
          }
        }

        if (all(is.na(auc_vals))) {
          best_idx <- 1
          best_auc <- NA
        } else {
          best_idx <- which.max(auc_vals)
          best_auc <- auc_vals[best_idx]
        }

        best_type <- candidate_types[best_idx]
        new_predicted[i] <- best_type
        new_aucs[i] <- best_auc

        alt_types <- candidate_types[-best_idx]
        alt_aucs <- auc_vals[-best_idx]
        alt_strs <- character(0)

        for (k in seq_along(alt_types)) {
          alt_strs[k] <- paste0(alt_types[k], " (",round(alt_aucs[k], digits = 7), ")")
        }
        new_alt_list[i] <- paste(alt_strs, collapse = " ; ")
      }

      prediction_results$Predicted_cell_type <- new_predicted
      prediction_results$AUC <- new_aucs
      prediction_results$Alternative_cell_types <- new_alt_list

      message(paste0("\n","The predicted cell types were corrected by AUC values."))

    } else {
      message("Calculating AUC values for predicted cell type.")
      auc_values <- numeric(nrow(prediction_results))

      for (i in seq_len(nrow(prediction_results))) {
        cluster_id <- prediction_results$cluster_col[i]
        cell_type <- prediction_results$Predicted_cell_type[i]

        if (!cell_type %in% names(valid_genes_list)) {
          warning(paste("Skipping AUC for cluster", cluster_id, ": No valid genes for", cell_type))
          auc_values[i] <- NA
          next
        }
        features <- valid_genes_list[[cell_type]]

        all_cells <- colnames(seurat_obj)
        expr_data <- FetchData(seurat_obj, vars = features, cells = all_cells)
        cell_scores <- rowMeans(expr_data, na.rm = TRUE)

        labels <- seurat_obj@meta.data[all_cells, cluster_col] == cluster_id
        if (length(unique(labels)) < 2) {
          warning(paste("Skipping AUC for cluster", cluster_id, ": Only one class present"))
          auc_values[i] <- NA
        } else {
          auc_values[i] <- fastAUC(cell_scores, labels)
        }
      }
      prediction_results$AUC <- auc_values
    }
  } else {
    prediction_results$AUC <- NA
  }

  prediction_results <- prediction_results[, c("cluster_col", "Predicted_cell_type", "AUC", "Alternative_cell_types")]

  heatmap_plot <- p

  return(list(Expression_list = expr_list,
              Expression_scores_matrix = scores_matrix,
              Probability_matrix = probability_matrix,
              Prediction_results = prediction_results,
              Heatmap_plot = heatmap_plot))
}
