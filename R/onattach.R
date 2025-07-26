.onAttach <- function(libname, pkgname) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    message <- paste(
      "Please cite: Wang Z (2025). ",
      crayon::italic("Marker-Based Package for Single-Cell and Spatial-Transcriptomic Annotation"),
      ". R package version ", crayon::bold("v1.0.3"),
      ". Available at: https://github.com/Zhaoqing-wang/SlimR"
    )
  } else {
    message <- paste(
      "Please cite: Wang Z (2025). Marker-Based Package for Single-Cell and Spatial-Transcriptomic Annotation. ",
      "R package version v1.0.3. Available at: https://github.com/Zhaoqing-wang/SlimR"
    )
  }

  packageStartupMessage(message)
}
