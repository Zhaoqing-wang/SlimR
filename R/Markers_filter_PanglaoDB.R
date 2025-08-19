#' Create Marker_list from the PanglaoDB database
#'
#' @param df Standardized PanglaoDB database. It is read as data(PanglaoDB)
#'     in the SlimR library.
#' @param species_input Species information in PanglaoDB database. The default
#'     input is "Human" or "Mouse".The input can be retrieved by "PanglaoDB_table".
#'     For more information,please refer to https://panglaodb.se/ on PanglaoDB's
#'     official website.
#' @param organ_input Organ type information in the PanglaoDB database.
#'     The input can be retrieved by "PanglaoDB_table".For more information,
#'     please refer to https://panglaodb.se/ on PanglaoDB's official website.
#'
#' @returns The standardized "Marker_list" in the SlimR package
#' @export
#' @family Section_2_Standardized_Markers_List_Input
#'
#' @examples
#' PanglaoDB <- SlimR::PanglaoDB
#' Markers_list_panglaoDB <- Markers_filter_PanglaoDB(
#'     PanglaoDB,
#'     species_input = 'Human',
#'     organ_input = 'GI tract'
#'     )
#'
Markers_filter_PanglaoDB <- function(df, species_input, organ_input) {
  required_columns <- c("species", "official.gene.symbol", "cell.type",
                        "ubiquitousness.index", "organ", "sensitivity_human",
                        "sensitivity_mouse", "specificity_human", "specificity_mouse")
  if (!all(required_columns %in% colnames(df))) {
    stop("Data frame missing necessary columns! Make sure to include the following: ", paste(required_columns, collapse = ", "))
  }

  if (!(species_input %in% c("Human", "Mouse"))) {
    stop("species_input must be 'Human' or 'Mouse'")
  }

  species_values <- if (species_input == "Human") {
    c("Hs", "Mm Hs")
  } else {
    c("Mm", "Mm Hs")
  }

  filtered_df <- df[df$species %in% species_values & df$organ == organ_input, ]

  if (nrow(filtered_df) == 0) {
    stop("The filter result is empty, please check the filter criteria.")
  }

  result_columns <- if (species_input == "Human") {
    c("cell.type", "official.gene.symbol", "ubiquitousness.index", "sensitivity_human", "specificity_human")
  } else {
    c("cell.type", "official.gene.symbol", "ubiquitousness.index", "sensitivity_mouse", "specificity_mouse")
  }

  result_df <- filtered_df[, result_columns]

  if (species_input == "Human") {
    result_df$`official.gene.symbol` <- toupper(result_df$`official.gene.symbol`)
  } else {
    result_df$`official.gene.symbol` <- gsub("(.)(.*)", "\\U\\1\\L\\2",
                                             tolower(result_df$`official.gene.symbol`),
                                             perl = TRUE)
  }

  split_result <- split(result_df, result_df$cell.type)
  clean_result <- lapply(split_result, function(sub_df) {
    sub_df[, -1]
  })

  return(clean_result)
}
