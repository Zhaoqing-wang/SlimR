#' Adaptive Parameter Tuning for Single-Cell Data Annotation in SlimR
#'
#' This function uses machine learning to automatically determine optimal 
#' min_expression and specificity_weight parameters for single-cell data analysis
#' based on dataset characteristics.
#'
#' @param object A Seurat object containing single-cell data
#' @param features Character vector of feature names (genes) to analyze
#' @param assay Name of assay to use (default: default assay)
#' @param cluster_col Column name in metadata containing cluster information
#' @param method Machine learning method: "rf" (random forest), "gbm" (gradient boosting),
#'               "svm" (support vector machine), or "ensemble" (default)
#' @param n_models Number of models for ensemble learning (default: 3)
#' @param return_model Whether to return trained model (default: FALSE)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item min_expression: Recommended expression threshold
#'   \item specificity_weight: Recommended specificity weight
#'   \item performance: Model performance metric (R-squared)
#'   \item dataset_features: Extracted dataset characteristics
#'   \item model: Trained model (if return_model = TRUE)
#' }
#'
#' @export
#' @family Section_3_Automated_Annotation_Workflow
#' 
#' @importFrom stats dist median predict runif sd aggregate
#' @importFrom utils installed.packages
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' params <- calculate_parameter(
#'   object = seurat_obj,
#'   features = c("CD3E", "CD4", "CD8A")
#' )
#' 
#' # Use with custom method
#' params <- calculate_parameter(
#'   object = seurat_obj,
#'   features = variable_genes,
#'   method = "rf"
#' )
#' }
#'
#' @export
calculate_parameter <- function(
    object,
    features,
    assay = NULL,
    cluster_col = NULL,
    method = "ensemble",
    n_models = 3,
    return_model = FALSE,
    verbose = TRUE
) {
  # Input validation
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Please install caret package: install.packages('caret')")
  }
  
  if (!inherits(object, "Seurat")) {
    stop("Input object must be a Seurat object")
  }
  
  if (length(features) < 5) {
    warning("Using fewer than 5 features may affect parameter tuning accuracy")
  }
  
  if (verbose) message("Extracting dataset features...")
  
  # Extract dataset characteristics for ML model
  dataset_features <- extract_dataset_features(object, features, assay, cluster_col)
  
  if (verbose) message("Generating training data...")
  
  # Generate synthetic training data based on empirical rules
  training_data <- generate_training_data(dataset_features)
  
  if (verbose) message("Training machine learning model...")
  
  # Train ML model to predict optimal parameters
  model_results <- train_parameter_model(
    training_data, 
    method = method, 
    n_models = n_models,
    verbose = verbose
  )
  
  # Predict optimal parameters for current dataset
  predicted_params <- predict_optimal_parameters(model_results$model, dataset_features)
  
  # Apply constraints and adjustments to predicted parameters
  final_params <- postprocess_parameters(predicted_params, dataset_features)
  
  if (verbose) {
    message("Parameter recommendation:")
    message("  min_expression: ", round(final_params$min_expression, 3))
    message("  specificity_weight: ", round(final_params$specificity_weight, 3))
    message("  Model performance (R-squared): ", round(model_results$performance, 3))
  }
  
  # Prepare results
  result <- list(
    min_expression = final_params$min_expression,
    specificity_weight = final_params$specificity_weight,
    performance = model_results$performance,
    dataset_features = dataset_features
  )
  
  if (return_model) {
    result$model <- model_results$model
  }
  
  return(result)
}

#' Extract Dataset Characteristics for Machine Learning
#'
#' Computes various statistical features from single-cell data that are used
#' as input for the parameter prediction model.
#'
#' @param object Seurat object
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
extract_dataset_features <- function(object, features, assay = NULL, cluster_col = NULL) {
  assay <- if (is.null(assay)) Seurat::DefaultAssay(object) else assay
  Seurat::DefaultAssay(object) <- assay
  
  cells <- unlist(Seurat::CellsByIdentities(object = object))
  data.features <- Seurat::FetchData(object = object, vars = features, cells = cells)
  
  # Assign cluster identities
  if (!is.null(cluster_col)) {
    data.features$id <- object@meta.data[cells, cluster_col, drop = TRUE]
  } else {
    data.features$id <- Seurat::Idents(object = object)[cells]
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
  
  # Compute comprehensive dataset characteristics
  dataset_features <- list(
    # Basic dataset properties
    n_genes = length(features),
    n_cells = nrow(data.features),
    n_clusters = length(unique(data.features$id)),
    
    # Global expression characteristics
    global_mean_expression = mean(as.matrix(data.features[, features])),
    global_zero_fraction = mean(as.matrix(data.features[, features]) == 0),
    expression_sparsity = mean(apply(data.features[, features], 2, function(x) mean(x == 0))),
    
    # Gene variability measures
    mean_gene_cv = mean(expression_stats["cv_expr", ], na.rm = TRUE),
    sd_gene_cv = stats::sd(expression_stats["cv_expr", ], na.rm = TRUE),
    
    # Cluster separation metrics
    cluster_variability = calculate_cluster_variability(data.features, features),
    
    # Distribution characteristics
    expression_skewness = calculate_expression_skewness(data.features[, features]),
    batch_effect_score = estimate_batch_effect(object, assay)
  )
  
  return(dataset_features)
}

#' Calculate Cluster Variability
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

#' Calculate Expression Distribution Skewness
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

#' Estimate Batch Effect Strength
#'
#' Roughly estimates the potential impact of batch effects
#' using available metadata.
#'
#' @param object Seurat object
#' @param assay Assay name
#'
#' @return Batch effect score (0 indicates no detectable batch effect)
#'
#' @family Section_1_Functions_Use_in_Package
#' 
estimate_batch_effect <- function(object, assay) {
  # Simple batch effect estimation
  if ("batch" %in% colnames(object@meta.data)) {
    # If batch information is available, use simplified approach
    # without requiring bluster package
    batch_groups <- unique(object@meta.data$batch)
    if (length(batch_groups) > 1) {
      # Return a simple metric based on batch group count
      return(length(batch_groups) * 0.1)
    }
  }
  return(0)  # Default: no batch effect detected
}

#' Generate Training Data for Machine Learning Model
#'
#' Creates synthetic training data based on empirical rules about
#' optimal parameter relationships with dataset characteristics.
#'
#' @param dataset_features List of actual dataset characteristics
#' @param n_samples Number of synthetic samples to generate
#'
#' @return Data frame with synthetic features and optimal parameter targets
#'
#' @family Section_1_Functions_Use_in_Package
#' 
#' @importFrom stats runif
#' 
generate_training_data <- function(dataset_features, n_samples = 1000) {
  set.seed(123)  # For reproducible synthetic data generation
  
  training_data <- data.frame(
    # Simulate diverse dataset characteristics
    n_genes = stats::runif(n_samples, 50, 5000),
    n_cells = stats::runif(n_samples, 100, 50000),
    n_clusters = sample(2:20, n_samples, replace = TRUE),
    global_mean_expression = stats::runif(n_samples, 0.1, 5),
    global_zero_fraction = stats::runif(n_samples, 0.1, 0.9),
    expression_sparsity = stats::runif(n_samples, 0.1, 0.95),
    mean_gene_cv = stats::runif(n_samples, 0.5, 3),
    sd_gene_cv = stats::runif(n_samples, 0.1, 1),
    cluster_variability = stats::runif(n_samples, 0.1, 10),
    expression_skewness = stats::runif(n_samples, 0.5, 5),
    batch_effect_score = stats::runif(n_samples, -1, 1)
  )
  
  # Generate target parameters using empirical relationships
  training_data$optimal_min_expression <- with(training_data, {
    base_min <- 0.05 + 0.15 * global_zero_fraction
    adj_min <- base_min * (1 + 0.2 * expression_skewness)
    pmin(pmax(adj_min, 0.01), 0.3)  # Constrain to reasonable range
  })
  
  training_data$optimal_specificity_weight <- with(training_data, {
    base_weight <- 1 + 2 * cluster_variability / (1 + expression_sparsity)
    adj_weight <- base_weight * (1 + 0.5 * mean_gene_cv)
    pmin(pmax(adj_weight, 0.5), 8)  # Constrain to reasonable range
  })
  
  # Add realistic noise to simulate real-world variation
  training_data$optimal_min_expression <- training_data$optimal_min_expression * 
    stats::runif(n_samples, 0.8, 1.2)
  training_data$optimal_specificity_weight <- training_data$optimal_specificity_weight * 
    stats::runif(n_samples, 0.8, 1.2)
  
  return(training_data)
}

#' Train Parameter Prediction Model
#'
#' Trains machine learning models to predict optimal parameters
#' based on dataset characteristics.
#'
#' @param training_data Data frame with features and target parameters
#' @param method Machine learning method to use
#' @param n_models Number of models for ensemble learning
#' @param verbose Whether to print training progress
#'
#' @family Section_1_Functions_Use_in_Package
#' 
#' @return List containing trained model and performance metrics
#' 
train_parameter_model <- function(training_data, method = "ensemble", n_models = 3, verbose = TRUE) {
  set.seed(123)
  
  feature_columns <- setdiff(colnames(training_data), 
                            c("optimal_min_expression", "optimal_specificity_weight"))
  
  if (method == "ensemble") {
    # Ensemble approach: combine multiple model types
    models <- list()
    performances <- numeric(n_models)
    
    for (i in 1:n_models) {
      if (verbose) message("Training ensemble model ", i, "/", n_models)
      
      # Use different data subsets for diversity
      train_idx <- sample(nrow(training_data), nrow(training_data) * 0.8)
      train_subset <- training_data[train_idx, ]
      
      # Alternate between different model types
      if (i %% 3 == 1) {
        # Random Forest
        model <- caret::train(
          x = train_subset[, feature_columns],
          y = data.frame(
            min_expression = train_subset$optimal_min_expression,
            specificity_weight = train_subset$optimal_specificity_weight
          ),
          method = "rf",
          trControl = caret::trainControl(method = "cv", number = 5),
          verbose = FALSE
        )
      } else if (i %% 3 == 2) {
        # Gradient Boosting Machine
        model <- caret::train(
          x = train_subset[, feature_columns],
          y = data.frame(
            min_expression = train_subset$optimal_min_expression,
            specificity_weight = train_subset$optimal_specificity_weight
          ),
          method = "gbm",
          trControl = caret::trainControl(method = "cv", number = 5),
          verbose = FALSE
        )
      } else {
        # Support Vector Machine
        model <- caret::train(
          x = train_subset[, feature_columns],
          y = data.frame(
            min_expression = train_subset$optimal_min_expression,
            specificity_weight = train_subset$optimal_specificity_weight
          ),
          method = "svmRadial",
          trControl = caret::trainControl(method = "cv", number = 5),
          verbose = FALSE
        )
      }
      
      models[[i]] <- model
      performances[i] <- max(model$results$Rsquared, na.rm = TRUE)
    }
    
    # Select best performing model
    best_idx <- which.max(performances)
    final_model <- models[[best_idx]]
    performance <- performances[best_idx]
    
  } else {
    # Single model approach
    if (verbose) message("Training ", method, " model...")
    
    model_method <- switch(method,
      "rf" = "rf",
      "gbm" = "gbm", 
      "svm" = "svmRadial",
      stop("Unsupported machine learning method: ", method)
    )
    
    final_model <- caret::train(
      x = training_data[, feature_columns],
      y = data.frame(
        min_expression = training_data$optimal_min_expression,
        specificity_weight = training_data$optimal_specificity_weight
      ),
      method = model_method,
      trControl = caret::trainControl(method = "cv", number = 5),
      verbose = FALSE
    )
    
    performance <- max(final_model$results$Rsquared, na.rm = TRUE)
  }
  
  return(list(model = final_model, performance = performance))
}

#' Predict Optimal Parameters Using Trained Model
#'
#' Applies the trained machine learning model to predict optimal
#' parameters for the current dataset.
#'
#' @param model Trained machine learning model
#' @param dataset_features Extracted characteristics of current dataset
#'
#' @return List containing predicted min_expression and specificity_weight
#'
#' @family Section_1_Functions_Use_in_Package
#' 
#' @importFrom stats predict
#' 
predict_optimal_parameters <- function(model, dataset_features) {
  # Convert features to data frame format expected by model
  feature_df <- as.data.frame(dataset_features)
  feature_df <- feature_df[names(dataset_features)]  # Maintain consistent order
  
  # Generate predictions
  predicted <- stats::predict(model, newdata = feature_df)
  
  return(list(
    min_expression = predicted[1, "min_expression"],
    specificity_weight = predicted[1, "specificity_weight"]
  ))
}

#' Post-process Predicted Parameters
#'
#' Applies constraints and dataset-specific adjustments to ensure
#' predicted parameters are within reasonable ranges.
#'
#' @param predicted_params List of raw predicted parameters
#' @param dataset_features Characteristics of current dataset
#'
#' @family Section_1_Functions_Use_in_Package
#' 
#' @return List of finalized parameters after post-processing
#' 
postprocess_parameters <- function(predicted_params, dataset_features) {
  min_expression <- predicted_params$min_expression
  specificity_weight <- predicted_params$specificity_weight
  
  # Apply dataset-specific adjustments
  if (dataset_features$global_zero_fraction > 0.8) {
    # Higher threshold for sparse datasets
    min_expression <- min_expression * 1.2
  }
  
  if (dataset_features$cluster_variability < 1) {
    # Higher weight for datasets with poor cluster separation
    specificity_weight <- specificity_weight * 1.3
  }
  
  # Enforce reasonable parameter bounds
  min_expression <- pmin(pmax(min_expression, 0.01), 0.3)
  specificity_weight <- pmin(pmax(specificity_weight, 0.5), 8)
  
  return(list(
    min_expression = min_expression,
    specificity_weight = specificity_weight
  ))
}

#' Enhanced Probability Calculator with Adaptive Parameter Tuning
#'
#' Wrapper function that automatically determines optimal parameters
#' before calculating cluster probabilities.
#'
#' @param object Seurat object
#' @param features Features to analyze
#' @param assay Assay name
#' @param cluster_col Cluster column name  
#' @param auto_tune Whether to automatically tune parameters (default: TRUE)
#' @param min_expression Manual expression threshold (optional)
#' @param specificity_weight Manual specificity weight (optional)
#' @param tuning_method Machine learning method for parameter tuning
#' @param verbose Whether to print progress messages
#'
#' @return Results from calculate_probability with optional tuning information
#'
#' @family Section_1_Functions_Use_in_Package
#' 
#' @examples
#' \dontrun{
#' # Automatic parameter tuning
#' result <- calculate_probability_adaptive(
#'   object = seurat_obj,
#'   features = marker_genes
#' )
#' 
#' # Manual parameters (no tuning)
#' result <- calculate_probability_adaptive(
#'   object = seurat_obj,
#'   features = marker_genes,
#'   auto_tune = FALSE,
#'   min_expression = 0.15,
#'   specificity_weight = 2
#' )
#' }
#'
#' @export
calculate_probability_adaptive <- function(
    object,
    features,
    assay = NULL,
    cluster_col = NULL,
    auto_tune = TRUE,
    min_expression = NULL,
    specificity_weight = NULL,
    tuning_method = "ensemble",
    verbose = TRUE
) {
  
  if (auto_tune || is.null(min_expression) || is.null(specificity_weight)) {
    if (verbose) message("Performing automatic parameter tuning...")
    
    tuned_params <- calculate_parameter(
      object = object,
      features = features,
      assay = assay,
      cluster_col = cluster_col,
      method = tuning_method,
      verbose = verbose
    )
    
    # Use tuned parameters unless manually specified
    if (is.null(min_expression)) {
      min_expression <- tuned_params$min_expression
    }
    if (is.null(specificity_weight)) {
      specificity_weight <- tuned_params$specificity_weight
    }
  }
  
  if (verbose) {
    message("Using parameters: min_expression = ", round(min_expression, 3), 
            ", specificity_weight = ", round(specificity_weight, 3))
  }
  
  # Call original function with determined parameters
  result <- calculate_probability(
    object = object,
    features = features,
    assay = assay,
    cluster_col = cluster_col,
    min_expression = min_expression,
    specificity_weight = specificity_weight
  )
  
  # Add tuning information if automatic tuning was performed
  if (auto_tune) {
    result$tuning_parameters <- list(
      min_expression = min_expression,
      specificity_weight = specificity_weight,
      tuning_method = tuning_method
    )
  }
  
  return(result)
}