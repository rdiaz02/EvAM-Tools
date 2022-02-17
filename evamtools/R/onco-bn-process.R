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



## \arguments{
## \item{data}{Cross secitonal data. Matrix of genes (columns)
## and individuals (rows)}

## \item{model}{One of "DBN" (default) or "CBN".}

## \item{algorithm}{"DP" (default) or "GA". See \code{\link[OncoBN]{fitCPN}}.}

## \item{epsilon}{Penalty term. See \code{\link[OncoBN]{fitCPN}}.}

## \item{silent}{Whether to show messages produced by OncoBN.}
## }
## \value{
## A list with  a data.frame with From-To edges and associated thetas, the
## thetas, the likelihood and epsilon of the model, and the predicted
## probabilities of genotypes.
## }

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

    ## And single ORs are confusing
    dbn_out$Relation <- vapply(
        dbn_out$To,
        function(x) ps_v[x],
        "some_string"
    )
    
    return(list(
        edges = dbn_out
      , thetas = thetas
      , likelihood = fit$score
      , epsilon = fit$epsilon
      , model = model
      , parent_set = ps_v
      , predicted_genotype_freqs = DBN_prob_genotypes(fit, colnames(data))
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


