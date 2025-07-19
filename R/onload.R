.onLoad <- function(libname, pkgname) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    message(paste(
      "Please cite: Wang Z (2025). ",
      crayon::italic("SlimR: Marker-Based Package for Single-Cell and ST Annotation"),
      ". R package version ", crayon::bold("v1.0.2"),
      ". Available at: ", crayon::url("https://github.com/Zhaoqing-wang/SlimR")
    ))
  } else {
    message(paste(
      "Please cite: Wang Z (2025). SlimR: Marker-Based Package for Single-Cell and ST Annotation. ",
      "R package version v1.0.2 Available at: https://github.com/Zhaoqing-wang/SlimR"
    ))
  }
}
