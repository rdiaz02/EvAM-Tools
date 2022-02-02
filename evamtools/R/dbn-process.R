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

## Commented out to avoid warning about fitCPN, until this
## gets fixed in Nicol's original code: they should be using
## importFrom: https://github.com/phillipnicol/OncoBN/issues/2

#' Run DBN on data
#' 
#' @param data data.frame with cross sectional data
#' @return list$edges Data.frame with From, To, Edge and Theta.
#' @return list$thetas Thetas of the best fit
#' @return list$likelihood Float. Likelihood of the fit
do_DBN <- function(data) {
  invisible(capture.output(out <- fitCPN(data, algorithm = "DP")))
  ## thetas <- OncoBN:::inferTheta(data, out)
  thetas <- evam_inferTheta(data, out)
  dbn_out <- igraph::as_data_frame(out$graph)
  dbn_out$from[dbn_out$from == "WT"] <- "Root"
  dbn_out$edges <- mapply(function(x,y){paste(x, "->", y, sep="")}, dbn_out$from, dbn_out$to)
  dbn_out$thetas <- rep(unlist(thetas$thetas), vapply(thetas$parents, length, integer(1)))
  colnames(dbn_out) <- c("From", "To", "Edge", "Thetas")
  return(list(
    edges = dbn_out
              , thetas = thetas
              , likelihood = out$score
              ))
}
