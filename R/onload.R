.onLoad <- function(libname, pkgname) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    message(paste("Please cite: Wang Z (2025). ",
                  crayon::italic("SlimR: Marker-Based Package for Single-Cell and ST Annotation"),
                  ". R package version 1.0.2"))
  } else {
    message("Please cite: Wang Z (2025). SlimR: Marker-Based Package for Single-Cell and ST Annotation. R package version 1.0.2")
  }
}
