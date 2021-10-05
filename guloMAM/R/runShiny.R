#' @export
runExample <- function() {
  appDir <- system.file("shiny-examples", "guloMAM", package = "guloMAM")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `guloMAM`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}