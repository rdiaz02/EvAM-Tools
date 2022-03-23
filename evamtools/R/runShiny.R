## Copyright 2021, 2022 Pablo Herrera Nieto

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' @title Run the web application of evamtools
#' 
#' @description  Launch the server with the web based app
runShiny <- function(host="0.0.0.0", port = 3000) {
  appDir <- system.file("shiny-examples", "evamtools", package = "evamtools")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `evamtools`.", call. = FALSE)
  }

  # options(shiny.autoreload = TRUE, browser = "/usr/bin/google-chrome")
  shiny::runApp(appDir, port = port, host = host, display.mode = "normal")
}
