#' Create Marker_list from the Cellmarkers2 database
#'
#' @param df Standardized Cellmarkers2 database. It is read as data(Cellmarkers2)
#'     in the SlimR library.
#' @param species Species information in Cellmarkers2 database. The default
#'     input is "Human" or "Mouse".The input can be retrieved by "Cellmarkers2_table".
#'     For more information,please refer to http://117.50.127.228/CellMarker/ on
#'     Cellmarkers2's official website.
#' @param tissue_class Tissue_class information in Cellmarkers2 database.
#'     The input can be retrieved by "Cellmarkers2_table". For more information,
#'     please refer to http://117.50.127.228/CellMarker/ on Cellmarkers2's official
#'     website.
#' @param tissue_type Tissue_type information in Cellmarkers2 database.
#'     The input can be retrieved by "Cellmarkers2_table". For more information,
#'     please refer to http://117.50.127.228/CellMarker/ on Cellmarkers2's official
#'     website.
#' @param cancer_type Cancer_type information in Cellmarkers2 database.
#'     The input can be retrieved by "Cellmarkers2_table". For more information,
#'     please refer to http://117.50.127.228/CellMarker/ on Cellmarkers2's official
#'     website.
#' @param cell_type Cell_type information in Cellmarkers2 database.
#'     The input can be retrieved by "Cellmarkers2_table". For more information,
#'     please refer to http://117.50.127.228/CellMarker/ on Cellmarkers2's official
#'     website.
#'
#' @returns The standardized "Marker_list" in the SlimR package
#' @export
#' @family Standardized_Marker_list_Input
#'
#' @importFrom stats aggregate
#'
#' @examples
#' Cellmarker2 <- SlimR::Cellmarker2
#' Markers_list_Cellmarker2 <- Markers_filter_Cellmarker2(
#'     Cellmarker2,
#'     species = "Human",
#'     tissue_class = "Intestine",
#'     tissue_type = NULL,
#'     cancer_type = NULL,
#'     cell_type = NULL
#'     )
#'
Markers_filter_Cellmarker2 <- function(df,
                                       species = NULL,
                                       tissue_class = NULL,
                                       tissue_type = NULL,
                                       cancer_type = NULL,
                                       cell_type = NULL) {

  required_columns <- c("species", "tissue_class", "tissue_type", "cancer_type",
                        "cell_type")

  if (!all(required_columns %in% colnames(df))) {
    stop("Data frame missing necessary columns! Please make sure to include the following: ",
         paste(required_columns, collapse = ", "))
  }

  filters <- list()
  if (!is.null(species)) filters$species <- species
  if (!is.null(tissue_class)) filters$tissue_class <- tissue_class
  if (!is.null(tissue_type)) filters$tissue_type <- tissue_type
  if (!is.null(cancer_type)) filters$cancer_type <- cancer_type
  if (!is.null(cell_type)) filters$cell_type <- cell_type

  filtered_df <- df
  for (col in names(filters)) {
    if (col %in% colnames(filtered_df)) {
      filtered_df <- filtered_df[filtered_df[[col]] %in% filters[[col]], ]
    } else {
      warning(paste("Column", col, "does not exist in the data frame, filter ignored"))
    }
  }

  if (nrow(filtered_df) == 0) {
    stop("Filter result is empty, check criteria")
  }

  result_df <- filtered_df[, c("cell_name", "marker", "counts")]

  if (!is.null(species) && species %in% c("Human", "Mouse")) {
    if (species == "Human") {
      result_df$marker <- toupper(result_df$marker)
    } else {
      result_df$marker <- gsub("^(.)", "\\U\\1", tolower(result_df$marker), perl = TRUE)
    }
  }

  aggregated_df <- aggregate(counts ~ cell_name + marker, data = result_df, sum)
  split_result <- split(aggregated_df, aggregated_df$cell_name)
  clean_result <- lapply(split_result, function(sub_df) sub_df[, -1])

  return(clean_result)
}
