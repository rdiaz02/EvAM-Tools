#' @export
runShiny <- function() {
  appDir <- system.file("shiny-examples", "guloMAM", package = "guloMAM")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `guloMAM`.", call. = FALSE)
  }

  options(shiny.autoreload = TRUE, browser = "/usr/bin/google-chrome")
  shiny::runApp(appDir, display.mode = "normal")
}

# enableBookmarking = "url"