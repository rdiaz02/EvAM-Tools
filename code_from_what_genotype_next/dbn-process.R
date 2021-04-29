library(OncoBN)


do_DBN <- function(data){
  out <- fitCPN(data, algorithm="GA") 
  ## This version of fitCPN is modified in my repo to return theta
  ## TODO: modify to also use DP
  
  out$edgelist <- replace(out$edgelist, out$edgelist=="WT", "Root")
  parents <- out$edgelist[seq(1, length(out$edgelist), 2)]
  childs <- out$edgelist[seq(2, length(out$edgelist), 2)]
  edges <- paste(parents, childs, sep= " -> ")
  thetas <- out$thetas[seq(2, length(out$thetas))]
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


