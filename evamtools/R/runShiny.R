#' @export
runShiny <- function() {
  appDir <- system.file("shiny-examples", "evamtools", package = "evamtools")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `evamtools`.", call. = FALSE)
  }

  options(shiny.autoreload = TRUE, browser = "/usr/bin/google-chrome")
  shiny::runApp(appDir, port = 3000, host = "0.0.0.0", display.mode = "normal")
}

# enableBookmarking = "url"