.onAttach <- function(libname, pkgname) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    message <- paste(
      "Please cite: Wang Z (2025). ",
      crayon::italic("SlimR: Marker-Based Package for Single-Cell and Spatial-Transcriptomic Annotation."),
      " R package version", crayon::bold("1.0.7."),
      " Available at: https://github.com/Zhaoqing-wang/SlimR"
    )
  } else {
    message <- paste(
      "Please cite: Wang Z (2025). SlimR: Marker-Based Package for Single-Cell and Spatial-Transcriptomic Annotation.",
      "R package version 1.0.7. Available at: https://github.com/Zhaoqing-wang/SlimR"
    )
  }

  packageStartupMessage(message)
}
