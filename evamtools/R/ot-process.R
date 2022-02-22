## Copyright 2016, 2017, 2018 Ramon Diaz-Uriarte

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


options(mc.cores = 1L)
options(boot.parallel = "no")
options(boot.ncpus = 1L)

ot_consensus_sb <- function(x) {
    ## Taking the code from plot.boottree
    child <- x$original$child
    parent.num <- as.numeric(x$consensus)
    parent <- character()
    nmut <- length(child)
    for (i in 1:nmut) {
        number <- parent.num[i]
        if (i == 1) {
            parent[i] <- ""
        }
        else {
            parent[i] <- child[number]
        }
    }

    ## why was this here??
    ## mostfreq <- list(child = child, parent = parent,
    ##                  parent.num = parent.num)

    ## Now, extract the From and To in my usual way
    edges.matrix <- cbind(parent = parent,
                          child = child)
    ## Just removing the "" -> Root
    if(identical(edges.matrix[1, ], c(parent = "", child = "Root"))){
        edges.matrix <- edges.matrix[-1, , drop = FALSE]
    } else {
        stop("this is not expected; where is 1, 2 as Root?")
    }
    
    return(data.frame(From = edges.matrix[, "parent"],
                      To = edges.matrix[, "child"],
                      edge = paste(edges.matrix[, "parent"],
                                   edges.matrix[, "child"],
                                   sep = " -> "),
                      stringsAsFactors = FALSE))
}


## ## You can verify this function by running examples
## ## and doing plot(boot, draw.consensus = TRUE, minfreq = 100)
## ot_consensus_mine <- function(fit, boot) {
##     cn <- colnames(fit$data)
##     child <- boot$original$child
##     stopifnot(child[1] == "Root")
##     child <- child[-1]
##     consensus <- boot$consensus
##     stopifnot(consensus[1] == 0)
##     stopifnot(sum(consensus == 0) == 1)
##     parent <-cn[consensus[-1]]
##     return(data.frame(
##         From = parent,
##         To = child,
##         edge = paste(parent, child, sep = " -> "),
##         stringsAsFactors = FALSE
##     ))
## }
## ## Some checks
## check_ot_consensus <- function(data, nb) {
##     otf <- oncotree.fit(data)
##     otb <- bootstrap.oncotree(otf, R = nb, type = "nonparametric")
##     par(mfrow = c(2, 1))
##     plot(otb, draw.consensus = TRUE, minfreq = 100)
##     cat("\n mine\n")
##     print(mine <- ot_consensus_mine(otf, otb))
##     cat("\n sb\n")
##     print(otsb <- ot_consensus_sb(otb))
##     identical(otsb, mine)
## }
## ## and do, for instance
## ## check_ot_consensus(ov.cgh, 5)

ot_consensus <- ot_consensus_sb

ot_proc <- function(datax, nboot = 1000,
                    distribution.oncotree = TRUE,
                    with_errors_dist_ot = TRUE) {

    error.fun <- "std"
    message(" Starting ot.fit ", date())
    ot.fit <- try(oncotree.fit(datax))
    if(inherits(ot.fit, "try-error")) {
        error.fun <- "tryNULL"
        ot.fit <- try(oncotree.fit(datax, error.fun = NULL))
    }
    message(" Done ot.fit ", date())
    edges.matrix <- cbind(parent = ot.fit$parent$parent,
                          child = ot.fit$parent$child)
    edge.weights <- ot.fit$parent$est.weight
    
    if(identical(edges.matrix[1, ], c(parent = "", child = "Root"))){
        edges.matrix <- edges.matrix[-1, , drop = FALSE]
        edge.weights <- edge.weights[-1]
        if(is.null(edge.weights)) {
            ## no est.weight, probably because error.fun is NULL
            ## use observed
            ## Was not in version for Baseline, but was never needed as
            ## we never faced this issue.
            edge.weights <- ot.fit$parent$obs.weight[-1]
            error.fun <- paste(error.fun, "observed.weights", sep = ",")
        }
    } else {
        stop("this is not expected; where is 1, 2 as Root?")
    }

    ## Marginal freqs, observed and fitted.
    ## We will add them in the output associated with the "To" node,
    ## as these are the predicted frequencies of the receiving node.
    obs_marginal <- colSums(datax)/nrow(datax)
    pred_marginal <- try(marginal.distr(ot.fit))
    if(inherits(pred_marginal, "try-error")) {
        pred_marginal <- marginal.distr(ot.fit, with.errors = FALSE)
        error.fun <- paste(error.fun, "pred.marginal.no.error", sep = ",")
    }
    if(nboot > 0) {
        message(" Starting bootstrap.oncotree ", date())
        ot.boot <- bootstrap.oncotree(ot.fit,
                                      R = nboot, 
                                      type = 'nonparametric')
        message(" Done bootstrap.oncotree ", date())
        ot.boot.freqs <- ot.boot$parent.freq
        ## From print.boottree, to return freq of original tree
        orig.string <- paste(ot.boot$original$parent.num, collapse = ".")
        boot.idx <- match(orig.string, as.character(ot.boot$tree.list$Tree))
        ot.boot.original <- ot.boot$tree.list$Freq[boot.idx]/nboot

        ## Two paranoid checks
        if (!all(edges.matrix[, 1, drop = FALSE] %in% colnames(ot.boot.freqs)))
            stop("colnames ot.boot.freqs weird")
        if (!all(edges.matrix[, 2, drop = FALSE] %in% rownames(ot.boot.freqs)))
            stop("rownames ot.boot.freqs weird")
        boot.freq <- ot.boot.freqs[edges.matrix]/nboot
        consensus <- try(ot_consensus(ot.boot))
    } else {
        boot.freq <- consensus <- ot.boot.original <- NA
    }
        
    
    if (distribution.oncotree) {
        ## with many genotypes and/or large trees can lead to unexpectedly
        ## very long computing times. Removed for now.
        ## Well, I think that only happens if with.errors = TRUE
        ## which we probably don't want anyway.
        message(" Starting distribution.oncotree ", date())
        ## Observed and expected frequencies of genotypes
        est_genots <- distribution.oncotree(ot.fit,
                                            with.probs = TRUE,
                                            with.errors = with_errors_dist_ot)

        ## tt <- as.data.frame(datax)
        ## obs_genots <- aggregate(tt, by = tt, length)[1:(ncol(tt) + 1)]
        ## colnames(obs_genots)[ncol(obs_genots)] <- "Counts"
        obs_genots <- NA
        message(" Ending distribution.oncotree ", date())

        est_genots <- dist_oncotree_output_2_named_genotypes(est_genots)

        ## ## Give a named vector for the estimated genotypes
        ## ## There is a column called Root, unlike OncoBN
        ## gpnroot <- which(colnames(est_genots) == "Root")
        ## gpnfr <- which(colnames(est_genots) == "Prob")
        ## gpn_names <- genot_matrix_2_vector(est_genots[, -c(gpnroot, gpnfr)])
        ## est_genots <- as.vector(est_genots[, "Prob"])
        ## names(est_genots) <- gpn_names

        ## ## If with.errors = FALSE, there can be missing genotypes
        ## est_genots <- reorder_to_standard_order(est_genots)
        ## if (length(is.na(est_genots)))
        ##     est_genots[is.na(est_genots)] <- 0
        
        ## Not using it now
        ## message(" Starting observed vs expected, oncotree ", date())
        ## ## Observed and expected, 2 events. From vignette
        ## est2way <- t(data.matrix(est_genots[2:(ncol(est_genots) - 1)])) %*% diag(est_genots$Prob) %*%
        ##     data.matrix(est_genots[2:(ncol(est_genots) - 1)])
        ## obs2way <-t(ot.fit$data[,-1]) %*% ot.fit$data[,-1]/nrow(ot.fit$data)
        ## message(" Ending observed vs expected, oncotree ", date())
        est2way <- NA
        obs2way <- NA
    } else {
        est_genots <- NA
        est2way <- NA
        obs_genots <- NA
        obs2way <- NA
    }
    
    return(list(edges = data.frame(From = edges.matrix[, "parent"],
                      To = edges.matrix[, "child"],
                      edge = paste(edges.matrix[, "parent"],
                                   edges.matrix[, "child"],
                                   sep = " -> "),
                      OT_edgeBootFreq = boot.freq,
                      OT_edgeWeight = edge.weights,
                      OT_obsMarginal = obs_marginal[edges.matrix[, "child"]],
                      OT_predMarginal = pred_marginal[edges.matrix[, "child"]],
                      stringsAsFactors = FALSE),
                eps = ot.fit$eps,
                consensus = consensus,
                OT_error.fun  = error.fun,
                ot.boot.original = ot.boot.original, ## Frequency of original tree among boot
                predicted_genotype_freqs = est_genots
                ## , genots_observed = obs_genots,
                ## , two_way_predicted = est2way,
                ## , two_way_observed = obs2way
                      ))
}

# library(codetools)
# checkUsageEnv(env = .GlobalEnv)

## Give a named vector for the estimated genotypes
dist_oncotree_output_2_named_genotypes <- function(odt) {
    ## In the output of distribution.oncotree there
    ## is a column called Root, unlike OncoBN
    gpnroot <- which(colnames(odt) == "Root")
    gpnfr <- which(colnames(odt) == "Prob")
    gpn_names <- genot_matrix_2_vector(odt[, -c(gpnroot, gpnfr)])
    odt <- as.vector(odt[, "Prob"])
    names(odt) <- gpn_names

    ## If with.errors = FALSE, there can be missing genotypes
    odt <- reorder_to_standard_order(odt)
    if (length(is.na(odt)))
        odt[is.na(odt)] <- 0
    return(odt)
}

