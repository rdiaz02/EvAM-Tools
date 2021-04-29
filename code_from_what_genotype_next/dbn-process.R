library(OncoBN)


do_DBN <- function(data){
  out <- fitCPN(data, algorithm="DP") 
  out$edgelist <- replace(out$edgelist, out$edgelist=="WT", "Root")
  parents <- out$edgelist[seq(1, length(out$edgelist), 2)]
  childs <- out$edgelist[seq(2, length(out$edgelist), 2)]
  edges <- paste(parents, childs, sep= " -> ")
  thetas <- unlist(inferTheta(data, out)$thetas)

  dbn_out <- data.frame(From=parents, To=childs, 
                        Edges=edges, Theta=thetas)
  
  tmp <- try(df_2_access_genots_and_graph(dbn_out[, c("From", "To")]))
  
  if(inherits(tmp, "try-error")) {
    stop("how is this happening? there was edges component!")
  } else {
    accessible_genots <- tmp$accessible_genots
    fgraph <- unrestricted_fitness_graph_sparseM(accessible_genots)
  }
  
  weights <- dbn_out[c("To", "Theta")]
  rownames(weights) <- dbn_out[, "To"]
  weighted_fgraph <- transition_fg_sparseM(fgraph, weights)
  
  ##TODO: include thetas for out-of-the-path mutations?
  
  trans_mat_genots <- rowScaleMatrix(weighted_fgraph)
  
  return(list(edges = dbn_out,
              weighted_fgraph = weighted_fgraph,
              trans_mat_genots = trans_mat_genots,
              likelihood = out$score
              ))
}


inferTheta <- function(df, fit){
  ## Copied from OncoDB package https://github.com/phillipnicol/OncoBN
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

