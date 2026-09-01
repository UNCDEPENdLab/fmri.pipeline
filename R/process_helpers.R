rscript_executable <- function() {
  executable <- if (identical(.Platform$OS.type, "windows")) {
    "Rscript.exe"
  } else {
    "Rscript"
  }
  file.path(R.home("bin"), executable)
}
