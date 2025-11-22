.onAttach <- function(libname, pkgname) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    message <- paste0(
      "Please cite: Wang Z (2025). ",
      crayon::italic("SlimR: Machine Learning-Assisted, Marker-Based Tool for Single-Cell and Spatial Transcriptomics Annotation."),
      " R package version", crayon::bold(" 1.0.9."),
      " Available at: https://github.com/Zhaoqing-wang/SlimR"
    )
  } else {
    message <- paste0(
      "Please cite: Wang Z (2025). SlimR: Machine Learning-Assisted, Marker-Based Tool for Single-Cell and Spatial Transcriptomics Annotation.",
      "R package version 1.0.9. Available at: https://github.com/Zhaoqing-wang/SlimR"
    )
  }

  packageStartupMessage(message)
}
