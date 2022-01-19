# library(OncoBN)
# source("../External-code/DBN/dbn-external.R")


#' Run DBN on data
#' 
#' @param data data.frame with cross sectional data
#' @return list$edges Data.frame with From, To, Edge and Theta.
#' @return list$thetas Thetas of the best fit
#' @return list$likelihood Float. Likelihood of the fit
do_DBN <- function(data) {
  invisible(capture.output(out <- fitCPN(data, algorithm = "DP")))
  thetas <- inferTheta(data, out)
  dbn_out <- create_data_frame_from_theta(thetas, out)
  return(list(
    edges = dbn_out
              , thetas = thetas
              , likelihood = out$score
              ))
}


create_data_frame_from_theta <- function(thetas, out) {  
  ## Ugly patch to get a clean data.frame from the output of inferTheta
  out$edgelist[out$edgelist == "WT"] <- "Root"
  from_node <- out$edgelist[seq(1, length(out$edgelist), by=2)]
  to_node <- out$edgelist[seq(2, length(out$edgelist), by=2)]
  edges <- mapply(function(x,y){paste(x, "->", y, sep="")}, from_node, to_node)
  names(edges) <- NULL

  # gene_names <- list()
  # theta_values <- list()
  # counter <- 1

  # for (parents in thetas$parents){
  #   for(parent in parents){
  #     if(!(parent %in% gene_names)){
  #       theta_values <- c(theta_values, thetas$thetas[[counter]])
  #       gene_names <- c(gene_names, parent)
  #     }
  #   }
  #   counter <- counter + 1
  # }

  # gene_names <- list()
  theta_values <- list()
  counter <- 1
  for (parents in thetas$parents){
    for(parent in parents){
        theta_values <- c(theta_values, thetas$thetas[[counter]])
        # gene_names <- c(gene_names, parent)
    }
    counter <- counter + 1
  }

  theta_values <- unlist(theta_values)
  dbn_out <- data.frame(from_node, to_node, edges, theta_values)
  # gene_names[gene_names == "WT"] <- "Root"
  # names(theta_values) <- gene_names

  # dbn_out <- data.frame(from_node, to_node, edges, unlist(theta_values[from_node]))
    
  colnames(dbn_out) <- c("From", "To", "Edge", "Thetas")
  return(dbn_out)
  # theta_values <- vapply(from_node 
  #   , function(x){}
  #   , float(1))

  # dbn_out <- data.frame(From = character(),
  #                       To = character(),
  #                       Edge = character(),
  #                       Theta = numeric(),
  #                       stringsAsFactors = FALSE)

  
  # thetas$parents <- lapply(thetas$parents,
  #                           function(i) str_replace(i, "WT", "Root"))

  # for (i in c(1:(length(gene_names)))) {
  #   parent <- thetas$parents[[i]]
  #   gene <- gene_names[[i]]
  #   theta <- thetas$thetas[[i]]

  #   if (length(parent) == 1) { ## May cause problems
  #     edge <- paste(parent, gene, sep = "->")
  #     new_row <- list(parent, gene, edge, theta)
  #     dbn_out <- rbind(dbn_out, new_row)
  #   } else {
  #     for (i in seq_len(parent)) {
  #       edge <- paste(parent[i], gene, sep = "->")
  #       new_row <- list(parent[i], gene, edge, theta)
  #       dbn_out <- rbind(dbn_out, new_row)
  #     }
  #   }
  # }
  # return(dbn_out)
}

# old_create_data_frame_from_theta <- function(data) {
#   ## Ugly patch to get a clean data.frame from the output of inferTheta
#   dbn_out <- data.frame(From = character(),
#                         To = character(),
#                         Edge = character(),
#                         Theta = numeric(),
#                         stringsAsFactors = FALSE)
#     data$parents <- lapply(data$parents,
#                            function(i) str_replace(i, "WT", "Root"))
#   for (i in c(1:(length(data$parents)))) {
#     parent <- data$parents[[i]]
#     child <- data$childs[[i]]
#     theta <- data$thetas[[i]]

#     if (length(parent) == 1) { ## May cause problems
#       edge <- paste(parent, child, sep = "->")
#       new_row <- list(parent, child, edge, theta)
#       rbind(dbn_out, new_row) -> dbn_out
#     } else {
#       for (i in seq_len(parent)) {
#         edge <- paste(parent[i], child[i], sep = "->")
#         new_row <- list(parent[i], child[i], edge, theta)
#         rbind(dbn_out, new_row) -> dbn_out
#       }

#     }
#   }
#   colnames(dbn_out) <- c("From", "To", "Edge", "Thetas")
#   return(dbn_out)
# }
