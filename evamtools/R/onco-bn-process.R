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
                      k = 3,
                      epsilon = min(colMeans(data)/2),
                      silent = TRUE) {

    if (silent)
        invisible(capture.output(fit <- fitCPN(data,
                                               model = model,
                                               algorithm = algorithm,
                                               k = k,
                                               epsilon = epsilon)))
    else
        fit <- fitCPN(data,
                      model = model,
                      algorithm = algorithm,
                      epsilon = epsilon)

    thetas <- fit$theta

    ## No longer need
    names(thetas) <- colnames(data)

    if (is.null(names(thetas))) stop("fit$theta should be named. ",
                                    "You are probably using and old ",
                                    "version of OncoBN.")
    stopifnot(identical(names(thetas), colnames(data)))

    ## No longer needed, though I do not like that
    ## the graph does not have WT first.
    ## Actually, not really fixed.
    ## Possible issues with the way the graph is returned
    ## Fixing it temporarily
    ## See https://github.com/phillipnicol/OncoBN/issues/3#issuecomment-1049074644
    adjm <- igraph::as_adjacency_matrix(
                        igraph::make_directed_graph(fit$edgelist))

    
    gn <- colnames(data)
    adjm <- adjm[c("WT", gn), c("WT", gn)]
    new_graph <- igraph::graph_from_adjacency_matrix(adjm)
    fit$graph <- new_graph

    dbn_out <- igraph::as_data_frame(fit$graph)
    colnames(dbn_out) <- c("From", "To")
    dbn_out$From[dbn_out$From == "WT"] <- "Root"
    dbn_out$edge <- mapply(function(x, y) { paste0(x, " -> ", y) },
                            dbn_out$From, dbn_out$To)
    dbn_out$theta <- thetas[dbn_out$To]

    dbn_out$Relation <- ifelse(model == "DBN", "OR", "AND")
    dbn_out$Relation[dbn_out$From == "Root"] <- "Single"
      
    ## The parent set is used for plotting
    ## Simple with a table? Oh well
    ps <- aggregate(Relation ~ To,
                    data = dbn_out,
                    FUN = function(x) {
                        if (length(x) == 1) return("Single")
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

    ## Resort what will be the edges component with Root on top
    ir <- which(dbn_out$From == "Root")
    inr <- setdiff(seq_len(nrow(dbn_out)), ir)
    dbn_out <- dbn_out[c(ir, inr), ]
    
    est_genots <- DBN_prob_genotypes(fit, colnames(data))
    ## Give a named vector for the predicted freqs of genotypes
    est_genots <- DBN_est_genots_2_named_genotypes(est_genots)

    return(list(
        edges = dbn_out
      , thetas = thetas
      , likelihood = fit$score
      , epsilon = fit$epsilon
      , model = model
      , parent_set = ps_v
      , predicted_genotype_freqs = est_genots
      , fit = fit
    ))
}



## From https://github.com/phillipnicol/OncoBN/issues/3#issuecomment-1033260046
DBN_prob_genotypes <- function(fit, gene_names) {
    n <- length(gene_names)
    genotypes <- expand.grid(replicate(n, 0:1, simplify = FALSE))
    colnames(genotypes) <- gene_names

    ## No longer needed, as Lik.genotype does the right thing now.
    ## Not quite. 
    ## Pre-check same order and proper naming
    ## See https://github.com/phillipnicol/OncoBN/issues/3#issuecomment-1048814030
    G <- t(as.matrix(as_adjacency_matrix(fit$graph)))
    stopifnot(colnames(genotypes) == colnames(G)[-1])
    stopifnot(colnames(G)[1] == "WT")

    ## ## FIXME: Lik.genotype does many identical operations for all genotypes.
    ## ## Extract that code, and call GA_Likelihood.
    ## genotypes$Prob <- apply(genotypes, 1,
    ##                         function(x) old_Lik.genotype(fit, x))

    ## Faster
    genotypes2 <- cbind(1, genotypes)
    theta.in <- c(1, fit$theta)
    genotypes$Prob <- apply(genotypes2, 1,
                            function(x) evam_GA_Likelihood(x, G, theta.in,
                                                      fit$epsilon,
                                                      fit$model))
    return(genotypes)
}

## Give a named vector for the predicted freqs of genotypes
DBN_est_genots_2_named_genotypes <- function(odt) {
    ## There is no column called Root, unlike OT
    gpnfr <- which(colnames(odt) == "Prob")
    gpn_names <- genot_matrix_2_vector(odt[, -gpnfr])
    odt <- as.vector(odt[, "Prob"])
    names(odt) <- gpn_names
    return(reorder_to_standard_order(odt))
}



## Lik.genotype in commit aa039d3
## old_Lik.genotype <- function(fit, genotype) {
##   x <- c(1,genotype)
##   G <- t(as.matrix(as_adjacency_matrix(fit$graph)))
##   theta.in <- c(1,fit$theta)
##   val <- GA_Likelihood(x, G, theta.in, fit$epsilon, fit$model)
##   return(val)
## }

