#' Title
#'
#' @param path
#'
#' @returns
#' @export
#'
#' @examples
read_excel_markers <- function(path) {
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
    sheet_dfs <- lapply(seq_along(sheets),function(i) read_excel(file,sheet=i,progress = TRUE))
    names(sheet_dfs) <- paste0(sheets)
    return(sheet_dfs)

  }

  all_dfs_list <- lapply(files,process_file)
  Marker_list <- unlist(all_dfs_list,recursive=FALSE)

  return(Marker_list)
}
