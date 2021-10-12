#' @export
runExample <- function() {
  appDir <- system.file("shiny-examples", "guloMAM", package = "guloMAM")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `guloMAM`.", call. = FALSE)
  }

  options(shiny.autoreload = TRUE)
  shiny::runApp(appDir, display.mode = "normal")
}

# enableBookmarking = "url"