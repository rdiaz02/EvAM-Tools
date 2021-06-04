library(OncoBN)


do_DBN <- function(data) {
  time_dbn <- system.time(
    out <- fitCPN(data, algorithm = "DP")
    )["elapsed"]

  thetas <- inferTheta(data, out) 
  dbn_out <- create_data_frame_from_theta(thetas)

  return(list(edges = dbn_out,
              thetas = thetas,
              likelihood = out$score,
              time = time_dbn
              ))
}

inferTheta <- function(df, fit){
  ## Modified from OncoDB package https://github.com/phillipnicol/OncoBN
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

  out <- list()
  out$parents <- parent_list
  out$thetas <- theta_list
  out$childs <- child_list
  return(out)
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

    if (length(parent) == 1){
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