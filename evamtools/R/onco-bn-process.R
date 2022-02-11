## Copyright 2021, 2022 Pablo Herrera Nieto, Ramon Diaz-Uriarte

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

do_OncoBN <- function(data,
                      model = "DBN",
                      algorithm = "DP",
                      epsilon = min(colMeans(data)/2),
                      silent = TRUE) {

    if(silent)
        invisible(capture.output(fit <- fitCPN(data,
                                               model = model,
                                               algorithm = algorithm,
                                               epsilon = epsilon)))
    else
        fit <- fitCPN(data,
                      model = model,
                      algorithm = algorithm,
                      epsilon = epsilon)
    thetas <- fit$theta
    names(thetas) <- colnames(data)
    dbn_out <- igraph::as_data_frame(fit$graph)
    colnames(dbn_out) <- c("From", "To")
    dbn_out$From[dbn_out$From == "WT"] <- "Root"
    dbn_out$Edge <- mapply(function(x, y) {paste0(x, " -> ", y)},
                            dbn_out$From, dbn_out$To)
    dbn_out$Thetas <- thetas[dbn_out$To]
    ## Allow double checking and running through HESBCN algorithm
    dbn_out$Relation <- ifelse(model == "DBN", "OR", "AND")
    dbn_out$Relation[dbn_out$From == "Root"] <- "Single"
    
    ## The parent set is used for plotting
    ps <- aggregate(Relation ~ To,
                    data = dbn_out,
                    FUN = function(x){
                        if(length(x) == 1) return("Single")
                        else return(x[1])}
                    )
    ps_v <- ps[, "Relation"]
    names(ps_v) <- ps[, "To"]

    
    return(list(
        edges = dbn_out
      , thetas = thetas
      , likelihood = fit$score
      , epsilon = fit$epsilon
      , model = model
      , parent_set = ps_v
      , genots_predicted = DBN_prob_genotypes(fit, colnames(data))
    ))
}



## From https://github.com/phillipnicol/OncoBN/issues/3#issuecomment-1033260046
DBN_prob_genotypes <- function(fit, gene_names) {
    n <- length(gene_names)
    genotypes <- expand.grid(replicate(n, 0:1, simplify = FALSE))
    colnames(genotypes) <- gene_names
    genotypes$Prob <- apply(genotypes, 1,
                            function(x) Lik.genotype(fit, x))
    return(genotypes)
}


