#' Adaptive Parameter Tuning for Single-Cell Data Annotation in SlimR
#'
#' This function automatically determines optimal min_expression and specificity_weight
#' parameters for single-cell data analysis based on dataset characteristics using
#' adaptive algorithms derived from empirical analysis of single-cell datasets.
#'
#' @param seurat_obj A Seurat object containing single-cell data
#' @param features Character vector of feature names (genes) to analyze. If NULL,
#'        will use highly variable features from the Seurat object.
#' @param assay Name of assay to use (default: default assay)
#' @param cluster_col Column name in metadata containing cluster information
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item min_expression: Recommended expression threshold
#'   \item specificity_weight: Recommended specificity weight
#'   \item dataset_features: Extracted dataset characteristics
#'   \item parameter_rationale: Explanation of parameter choices
#' }
#'
#' @export
#' @family Section_3_Automated_Annotation
#' 
#' @importFrom stats dist median sd aggregate quantile
#'
#' @examples
#' \dontrun{
#' SlimR_params <- Parameter_Calculate(
#'   seurat_obj = sce,
#'   features = c("CD3E", "CD4", "CD8A"),
#'   assay = "RNA",
#'   cluster_col = "seurat_clusters",
#'   verbose = TRUE
#'   )
#' }
#'
Parameter_Calculate <- function(
    seurat_obj,
    features = NULL,
    assay = NULL,
    cluster_col = NULL,
    verbose = TRUE
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input object must be a Seurat object")
  }
  
  assay <- if (is.null(assay)) Seurat::DefaultAssay(seurat_obj) else assay
  Seurat::DefaultAssay(seurat_obj) <- assay
  
  if (is.null(features) || length(features) == 0) {
    if (verbose) message("SlimR parameter calculate: No features provided, using variable features.")
    features <- Seurat::VariableFeatures(seurat_obj)
    if (length(features) == 0) {
      features <- head(rownames(seurat_obj[[assay]]), 2000)
    }
    features <- head(features, 500)
  }
  
  valid_features <- features[features %in% rownames(seurat_obj[[assay]])]
  if (length(valid_features) < 3) {
    warning("Fewer than 3 valid features found, results may be unreliable")
  }
  
  if (verbose) message("SlimR parameter calculate: Extracting dataset features from ", length(valid_features), " genes.")
  
  dataset_features <- extract_dataset_features(seurat_obj, valid_features, assay, cluster_col)
  
  if (verbose) message("SlimR parameter calculate: Computing adaptive parameters.")
  
  optimal_params <- compute_adaptive_parameters(dataset_features)
  
  if (verbose) {
    message("SlimR parameter calculate: Parameter recommendation: ")
    message("  min_expression: ", round(optimal_params$min_expression, 4))
    message("  specificity_weight: ", round(optimal_params$specificity_weight, 4))
    message("  Rationale: ", optimal_params$rationale)
  }
  
  result <- list(
    min_expression = optimal_params$min_expression,
    specificity_weight = optimal_params$specificity_weight,
    dataset_features = dataset_features,
    parameter_rationale = optimal_params$rationale
  )
  
  return(result)
}

#' Extract Dataset Characteristics for Machine Learning (Use in package)
#'
#' Computes various statistical features from single-cell data that are used
#' as input for the parameter prediction model.
#'
#' @param seurat_obj Seurat object
#' @param features Features to analyze
#' @param assay Assay name
#' @param cluster_col Cluster column name
#'
#' @return List of dataset characteristics including expression statistics,
#'         variability measures, and cluster properties
#'
#' @family Section_1_Functions_Use_in_Package
#' 
#' @importFrom stats dist median sd aggregate
#' 
extract_dataset_features <- function(seurat_obj, features, assay = NULL, cluster_col = NULL) {
  assay <- if (is.null(assay)) Seurat::DefaultAssay(seurat_obj) else assay
  Seurat::DefaultAssay(seurat_obj) <- assay
  
  cells <- unlist(Seurat::CellsByIdentities(object = seurat_obj))
  data.features <- Seurat::FetchData(object = seurat_obj, vars = features, cells = cells)
  
  # Assign cluster identities
  if (!is.null(cluster_col)) {
    data.features$id <- seurat_obj@meta.data[cells, cluster_col, drop = TRUE]
  } else {
    data.features$id <- Seurat::Idents(object = seurat_obj)[cells]
  }
  
  features <- setdiff(features, "id")
  
  # Compute gene expression statistics
  expression_stats <- sapply(features, function(gene) {
    expr <- data.features[[gene]]
    c(
      mean_expr = mean(expr),
      sd_expr = stats::sd(expr),
      zero_frac = mean(expr == 0),
      median_expr = stats::median(expr),
      cv_expr = stats::sd(expr) / (mean(expr) + 1e-6)  # Coefficient of variation
    )
  })
  
  dataset_features <- list(
    n_genes = length(features),
    n_cells = nrow(data.features),
    n_clusters = length(unique(data.features$id)),
    
    global_mean_expression = mean(as.matrix(data.features[, features])),
    global_zero_fraction = mean(as.matrix(data.features[, features]) == 0),
    expression_sparsity = mean(apply(data.features[, features], 2, function(x) mean(x == 0))),
    
    mean_gene_cv = mean(expression_stats["cv_expr", ], na.rm = TRUE),
    sd_gene_cv = stats::sd(expression_stats["cv_expr", ], na.rm = TRUE),
    
    cluster_variability = calculate_cluster_variability(data.features, features),
    
    expression_skewness = calculate_expression_skewness(data.features[, features]),
    batch_effect_score = estimate_batch_effect(seurat_obj, assay),
    
    expression_quantiles = stats::quantile(
      as.matrix(data.features[, features])[as.matrix(data.features[, features]) > 0], 
      probs = c(0.1, 0.25, 0.5, 0.75, 0.9), 
      na.rm = TRUE
    )
  )
  
  return(dataset_features)
}

#' Calculate Cluster Variability (Use in package)
#'
#' Measures the degree of separation between different cell clusters
#' based on expression patterns.
#'
#' @param data.features Data frame containing expression data and cluster labels
#' @param features Feature names to include in analysis
#'
#' @return Numeric value representing cluster separation strength
#'
#' @family Section_1_Functions_Use_in_Package
#' 
#' @importFrom stats dist aggregate
#' 
calculate_cluster_variability <- function(data.features, features) {
  cluster_means <- stats::aggregate(data.features[, features], 
                            by = list(cluster = data.features$id), 
                            mean)
  
  cluster_matrix <- as.matrix(cluster_means[, -1])
  
  if (nrow(cluster_matrix) > 1) {
    # Compute mean distance between cluster centroids
    dist_matrix <- stats::dist(cluster_matrix)
    variability <- mean(as.matrix(dist_matrix))
  } else {
    variability <- 0  # Only one cluster
  }
  
  return(variability)
}

#' Calculate Expression Distribution Skewness (Use in package)
#'
#' Computes the average skewness of gene expression distributions
#' across all features.
#'
#' @param expression_matrix Matrix of expression values
#'
#' @return Mean absolute skewness across all genes
#'
#' @family Section_1_Functions_Use_in_Package
#' 
calculate_expression_skewness <- function(expression_matrix) {
  skew_vals <- apply(expression_matrix, 2, function(x) {
    if (stats::sd(x) == 0) return(0)
    mean((x - mean(x))^3) / (stats::sd(x)^3)  # Fisher-Pearson coefficient of skewness
  })
  return(mean(abs(skew_vals), na.rm = TRUE))
}

#' Estimate Batch Effect Strength (Use in package)
#'
#' Roughly estimates the potential impact of batch effects
#' using available metadata.
#'
#' @param seurat_obj Seurat object
#' @param assay Assay name
#'
#' @return Batch effect score (0 indicates no detectable batch effect)
#'
#' @family Section_1_Functions_Use_in_Package
#' 
estimate_batch_effect <- function(seurat_obj, assay) {
  if ("batch" %in% colnames(seurat_obj@meta.data)) {
    batch_groups <- unique(seurat_obj@meta.data$batch)
    if (length(batch_groups) > 1) {
      return(length(batch_groups) * 0.1)
    }
  }
  return(0)
}

#' Compute Adaptive Parameters Based on Dataset Features (Use in package)
#'
#' Calculates optimal min_expression and specificity_weight parameters
#' using empirically-derived adaptive algorithms based on dataset characteristics.
#'
#' @param dataset_features List of dataset characteristics from extract_dataset_features()
#'
#' @return List containing min_expression, specificity_weight, and rationale
#'
#' @family Section_1_Functions_Use_in_Package
#'
#' @importFrom stats quantile
#'
compute_adaptive_parameters <- function(dataset_features) {
  
  sparsity <- dataset_features$global_zero_fraction
  mean_expr <- dataset_features$global_mean_expression
  cv <- dataset_features$mean_gene_cv
  cluster_var <- dataset_features$cluster_variability
  n_clusters <- dataset_features$n_clusters
  skewness <- dataset_features$expression_skewness
  expr_quantiles <- dataset_features$expression_quantiles
  
  rationale_parts <- character(0)
  
  if (!is.null(expr_quantiles) && length(expr_quantiles) >= 2) {
    base_min_expr <- expr_quantiles[1]
    if (base_min_expr < 0.01) base_min_expr <- 0.05
  } else {
    base_min_expr <- 0.1
  }
  
  if (sparsity > 0.85) {
    min_expression <- base_min_expr * 0.5
    rationale_parts <- c(rationale_parts, "High sparsity detected, lowering threshold")
  } else if (sparsity > 0.7) {
    min_expression <- base_min_expr * 0.8
    rationale_parts <- c(rationale_parts, "Moderate sparsity")
  } else if (sparsity < 0.3) {
    min_expression <- base_min_expr * 1.5
    rationale_parts <- c(rationale_parts, "Dense expression, raising threshold")
  } else {
    min_expression <- base_min_expr
  }
  
  if (skewness > 3) {
    min_expression <- min_expression * 1.2
    rationale_parts <- c(rationale_parts, "High skewness adjustment")
  }
  
  min_expression <- max(0.01, min(0.5, min_expression))
  
  if (cluster_var > 5) {
    base_weight <- 1.5
    rationale_parts <- c(rationale_parts, "Good cluster separation, moderate weight")
  } else if (cluster_var > 2) {
    base_weight <- 2.5
    rationale_parts <- c(rationale_parts, "Moderate cluster separation")
  } else if (cluster_var > 0.5) {
    base_weight <- 4.0
    rationale_parts <- c(rationale_parts, "Weak cluster separation, higher weight needed")
  } else {
    base_weight <- 5.0
    rationale_parts <- c(rationale_parts, "Poor cluster separation, high weight applied")
  }
  
  if (cv > 2) {
    specificity_weight <- base_weight * 0.8
    rationale_parts <- c(rationale_parts, "High gene variability, reducing weight")
  } else if (cv < 0.8) {
    specificity_weight <- base_weight * 1.3
    rationale_parts <- c(rationale_parts, "Low gene variability, increasing weight")
  } else {
    specificity_weight <- base_weight
  }
  
  if (n_clusters > 15) {
    specificity_weight <- specificity_weight * 1.1
  } else if (n_clusters < 5) {
    specificity_weight <- specificity_weight * 0.9
  }
  
  specificity_weight <- max(0.5, min(8, specificity_weight))
  
  rationale <- paste(rationale_parts, collapse = "; ")
  if (rationale == "") rationale <- "Standard parameters applied"
  
  return(list(
    min_expression = round(min_expression, 4),
    specificity_weight = round(specificity_weight, 4),
    rationale = rationale
  ))
}
