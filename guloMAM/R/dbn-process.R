library(OncoBN)

#' Run DBN on data
#' 
#' @param data data.frame with cross sectional data
#' @return list$edges Data.frame with From, To, Edge and Theta.
#' @return list$thetas Thetas of the best fit
#' @return list$likelihood Float. Likelihood of the fit
do_DBN <- function(data) {
  invisible(capture.output(out <- fitCPN(data, algorithm = "DP")))

  thetas <- inferTheta(data, out) 
  dbn_out <- create_data_frame_from_theta(thetas)

  return(list(edges = dbn_out
              , thetas = thetas
              , likelihood = out$score
              ))
}

#' Return theta from 
#' Modified from OncoDB package https://github.com/phillipnicol/OncoBN
#' 
#' @param df output from DBN
#' @param fit output from DBN
#' @return List with parent, children and thetas
inferTheta <- function(df, fit){
  df <- cbind(df, 1)
  colnames(df)[ncol(df)] <- "WT" 
  
  g <- make_directed_graph(fit$edgelist) 
  parent_list <- list()
  child_list <- list()
  theta_list <- list()
  for(i in 1:(ncol(df) - 1)) {
    parents <- incident(g, colnames(df)[i], mode = "in")
    parent_list[[i]] <- ends(g, parents)[,1]
    child_list[[i]] <- ends(g, parents)[,2]
     
    p_on <- as.data.frame(df[, parent_list[[i]] ])
    theta_list[[i]] <- length(which(df[rowSums(p_on) > 0, i] == 1))/nrow(as.data.frame(p_on[rowSums(p_on) > 0, ]))
    if(theta_list[[i]] > 1) {
      stop(theta_list[[i]])
    }
  }

  return(list(
    parents = parent_list
    , thetas = theta_list
    , childs = child_list
  ))
}

create_data_frame_from_theta <- function(data){  
  ## Ugly patch to get a clean data.frame from the output of inferTheta
  dbn_out <- data.frame(From = character(),
                        To = character(),
                        Edge = character(),
                        Theta = numeric(),
                        stringsAsFactors = FALSE)
  data$parents <- lapply(data$parents, function(i) str_replace(i, "WT", "Root"))
  for(i in c(1:(length(data$parents)))){
    parent<- data$parents[[i]]
    child <- data$childs[[i]]
    theta <- data$thetas[[i]]

    if (length(parent) == 1){ ## May cause problems 
      edge <- paste(parent, child, sep="->")
      new_row <- list(parent, child, edge, theta)
      rbind(dbn_out,new_row) -> dbn_out
    } else {
      for(i in 1:length(parent)){
        edge <- paste(parent[i], child[i], sep="->")
        new_row <- list(parent[i], child[i], edge, theta)
        rbind(dbn_out,new_row) -> dbn_out
      }

    }
  }
  colnames(dbn_out) <- c("From", "To", "Edge", "Thetas")
  return(dbn_out)
}