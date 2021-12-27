#' Return theta from 
#' Copied from https://github.com/phillipnicol/OncoBN
#' The package does not export this, but it would be perfect
#' 
#' @param df output from DBN
#' @param fit output from DBN
#' @return List with parent and thetas
inferTheta <- function(df, fit) {
  df <- cbind(df, 1)
  colnames(df)[ncol(df)] <- "WT" 
  
  g <- make_directed_graph(fit$edgelist) 
  parent_list <- list()
  theta_list <- list()
  for(i in 1:(ncol(df) - 1)) {
    parents <- incident(g, colnames(df)[i], mode = "in")
    parent_list[[i]] <- ends(g, parents)[,1]

    p_on <- as.data.frame(df[, parent_list[[i]] ])
    theta_list[[i]] <- length(which(df[rowSums(p_on) > 0, i] == 1))/nrow(as.data.frame(p_on[rowSums(p_on) > 0, ]))
    if(theta_list[[i]] > 1) {
      stop(theta_list[[i]])
    }
  }
  out <- list()
  out$parents <- parent_list
  out$thetas <- theta_list
  return(out)
}
