#' Create "Marker_list" from Excel files ".xlsx"
#'
#' @param path The path information of Marker files stored in ".xlsx" format.
#'     The Sheet name in the file is filled with cell type. The first line of
#'     each Sheet is the table head, the first column is filled with markers
#'     information, and the following column is filled with mertic information.
#' @param has_colnames Logical value indicating whether the first row contains 
#'     column names. If FALSE, the first column will be named "Markers" and 
#'     subsequent columns will be named "Col1", "Col2", etc.
#'
#' @returns The standardized "Marker_list" in the SlimR package.
#' @export
#' @family Section_2_Standardized_Markers_List
#'
#' @importFrom readxl excel_sheets
#' @importFrom tools file_ext file_path_sans_ext
#'
#' @examples
#' \dontrun{
#' Markers_list_Excel <- Read_excel_markers(
#'     "D:/Laboratory/Marker_load.xlsx"
#'     )
#' }
#'
Read_excel_markers <- function(path, has_colnames = TRUE) {
  if (!file.exists(path)) stop("Path does not exist:")
  file_info <- file.info(path)

  if (file_info$isdir) {
    files <- list.files(path, pattern = "\\.xlsx$", full.names = TRUE)
    if (length(files) == 0) return(list())
  } else {
    file_ext <- tolower(tools::file_ext(path))
    if (file_ext != "xlsx") stop("File must be in .xlsx format")
    files <- path
    
  }

  process_file <- function(file) {
    file_name <- tools::file_path_sans_ext(basename(file))
    sheets <- readxl::excel_sheets(file)
    
    sheet_dfs <- lapply(seq_along(sheets), function(i) {
      if (has_colnames) {
        df <- readxl::read_excel(file, sheet = i, col_names = TRUE, progress = TRUE)
      } else {
        df <- readxl::read_excel(file, sheet = i, col_names = FALSE, progress = TRUE)
        if (ncol(df) >= 1) {
          colnames(df)[1] <- "Markers"
        }
        if (ncol(df) > 1) {
          colnames(df)[2:ncol(df)] <- paste0("Col", 1:(ncol(df)-1))
        }
      }
      return(df)
    })
    
    names(sheet_dfs) <- sheets
    return(sheet_dfs)

  }

  all_dfs_list <- lapply(files, process_file)
  Marker_list <- unlist(all_dfs_list, recursive = FALSE)

  return(Marker_list)
}
