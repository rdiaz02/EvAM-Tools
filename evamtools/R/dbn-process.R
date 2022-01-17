## library(OncoBN)
## source("../External-code/DBN/dbn-external.R")


#' Run DBN on data
#' 
#' @param data data.frame with cross sectional data
#' @return list$edges Data.frame with From, To, Edge and Theta.
#' @return list$thetas Thetas of the best fit
#' @return list$likelihood Float. Likelihood of the fit
do_DBN <- function(data) {
  invisible(capture.output(out <- fitCPN(data, algorithm = "DP")))
  thetas <- inferTheta(data, out)
  dbn_out <- create_data_frame_from_theta(thetas, colnames(data))
  return(list(
    edges = dbn_out
              , thetas = thetas
              , likelihood = out$score
              ))
}


create_data_frame_from_theta <- function(thetas, gene_names) {  
  ## Ugly patch to get a clean data.frame from the output of inferTheta
  dbn_out <- data.frame(From = character(),
                        To = character(),
                        Edge = character(),
                        Theta = numeric(),
                        stringsAsFactors = FALSE)
    thetas$parents <- lapply(thetas$parents,
                             function(i) str_replace(i, "WT", "Root"))

  for (i in c(1:(length(gene_names)))) {
    parent <- thetas$parents[[i]]
    gene <- gene_names[[i]]
    theta <- thetas$thetas[[i]]

    if (length(parent) == 1) { ## May cause problems
      edge <- paste(parent, gene, sep = "->")
      new_row <- list(parent, gene, edge, theta)
      dbn_out <- rbind(dbn_out, new_row)
    } else {
      for (i in seq_len(parent)) {
        edge <- paste(parent[i], gene, sep = "->")
        new_row <- list(parent[i], gene, edge, theta)
        dbn_out <- rbind(dbn_out, new_row)
      }
    }
  }
  colnames(dbn_out) <- c("From", "To", "Edge", "Thetas")
  return(dbn_out)
}

old_create_data_frame_from_theta <- function(data) {
  ## Ugly patch to get a clean data.frame from the output of inferTheta
  dbn_out <- data.frame(From = character(),
                        To = character(),
                        Edge = character(),
                        Theta = numeric(),
                        stringsAsFactors = FALSE)
    data$parents <- lapply(data$parents,
                           function(i) str_replace(i, "WT", "Root"))
  for (i in c(1:(length(data$parents)))) {
    parent <- data$parents[[i]]
    child <- data$childs[[i]]
    theta <- data$thetas[[i]]

    if (length(parent) == 1) { ## May cause problems
      edge <- paste(parent, child, sep = "->")
      new_row <- list(parent, child, edge, theta)
      rbind(dbn_out, new_row) -> dbn_out
    } else {
      for (i in seq_len(parent)) {
        edge <- paste(parent[i], child[i], sep = "->")
        new_row <- list(parent[i], child[i], edge, theta)
        rbind(dbn_out, new_row) -> dbn_out
      }

    }
  }
  colnames(dbn_out) <- c("From", "To", "Edge", "Thetas")
  return(dbn_out)
}
