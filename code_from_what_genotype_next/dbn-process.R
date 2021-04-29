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
  return(dbn_out)
}

# 
# N <- 50
# na <- N
# nb <- N + round( 10 * runif(1))
# nab <- .6 * N + round( 10 * runif(1))
# n00 <- round( 10 * runif(1))
# dB <- matrix(
#   c(
#     rep(c(1, 0), na) 
#     , rep(c(0, 1), nb)
#     , rep(c(1, 1), nab)
#     , rep(c(0, 0), n00)
#   ), ncol = 2, byrow = TRUE
# )
# colnames(dB) <- LETTERS[1:2]
# 
# z <- fitCPN(dB, model="DBN", algorithm="DP")

