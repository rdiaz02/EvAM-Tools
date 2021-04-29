library(OncoBN)


do_DBN <- function(data){
  out <- fitCPN(data, algorithm="GA")
  out$edgelist <- replace(out$edgelist, out$edgelist=="WT", "Root")
  parents <- out$edgelist[seq(1, length(out$edgelist), 2)]
  childs <- out$edgelist[seq(2, length(out$edgelist), 2)]
  edges <- paste(parents, childs, sep= " -> ")
  thetas <- out$thetas[seq(2, length(out$thetas))]
  browser()
  dbn_out <- data.frame(From=parents, To=childs, 
                        Edges=edges, Theta=thetas)
  return(list(edges=dbn_out))
}


