#' Create "Marker_list" from Excel files ".xlsx"
#'
#' @param path The path information of Marker files stored in ".xlsx" format.
#'     The Sheet name in the file is filled with cell type. The first line of
#'     each Sheet is the table head, the first column is filled with markers
#'     information, and the following column is filled with mertic information.
#'
#' @returns The standardized "Marker_list" in the SlimR package.
#' @export
#' @family Standardized_Marker_list_Input
#'
#' @importFrom readxl excel_sheets
#'
#' @examples
#' \dontrun{
#' Markers_list_Excel <- Read_excel_markers(
#'     "D:/Laboratory/Marker_load.xlsx"
#'     )
#'     }
#'
Read_excel_markers <- function(path) {
  if (!file.exists(path)) stop("Path does not exist:")
  file_info <- file.info(path)

  if (file_info$isdir) {
    files <- list.files(path,pattern="\\.xlsx$",full.names=TRUE)
    if (length(files)==0) return(list())

  } else {
    file_ext <- tolower(tools::file_ext(path))
    if (file_ext!="xlsx") stop("File must be in .xlsx format")
    files <- path

  }

  process_file <- function(file) {
    file_name <- tools::file_path_sans_ext(basename(file))
    sheets <- excel_sheets(file)
    sheet_dfs <- lapply(seq_along(sheets),function(i) readxl::read_excel(file,sheet=i,progress = TRUE))
    names(sheet_dfs) <- paste0(sheets)
    return(sheet_dfs)

  }

  all_dfs_list <- lapply(files,process_file)
  Marker_list <- unlist(all_dfs_list,recursive=FALSE)

  return(Marker_list)
}
