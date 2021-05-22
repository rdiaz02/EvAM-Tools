## A simplified version of code-all-methods-2-trans-matrix-max-genes-15.R
## Not using MCCBN, CAPRESE, or CAPRI
## But it gives more output of OT, MHN, CBN
## And we produce plots

## I comment out methods and tests

## Copyright 2016, 2017, 2018, 2020 Ramon Diaz-Uriarte

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



## Given a data set (patients as rows, genes as columns) return the
## transition matrices between genotypes according to all methods.

## Done by function: all_methods_2_trans_mat
## almost at the bottom.


## You need to have
## - Schills code
## - MCCBN installed: from git, then R CMD INSTALL: https://github.com/cbg-ethz/MC-CBN
## - CBN: use my version with fixes, in the repo:
##     cd to ct-cbn-0.1.04b-with-rdu-bug-fix-write-lambda-likelihood
##     then:
##       ./configure
##       make
##       and cp, mv, ln the two binaries to a place in the path
##      BEWARE! sometimes weird things can happen (make not working, etc)
##        if you have stale files. 
##        If things break, the simplest is to checkout a new copy from the repo
##        and do the ./configure , make dance there.
##        Probably you need to have autoconf-archive installed.
## - OT, and CAPRIand CAPRESE (the corresponding packages Oncotree and TRONCO)

date()

## Set it to TRUE if you want to load MCCBN, which requires
## having it installed. It then tests the MCCBN functionality too.
MCCBN_INSTALLED <- TRUE


## Since we are not parallelizing here, you might want to set
thiscores <- 1
## export OPENBLAS_NUM_THREADS=thiscores
## export OMP_NUM_THREADS=thiscores
## with the 36 substituted by the number of your CPUs
## in the script that launches analyses (unless those runs are parallelized)
## I also do it here for R itself

library(RhpcBLASctl)
## RhpcBLASctl::blas_set_num_threads(thiscores)
## RhpcBLASctl::omp_set_num_threads(thiscores)


## library(compiler)
## enableJIT(3)
## setCompilerOptions(optimize = 3)
## ## ## library(parallel)

library(dplyr)
library(Oncotree)
library(readr)
library(data.table)
library(igraph)
## library(TRONCO)
library(parallel)
library(foreach)
## library(doParallel)
library(Rgraphviz)
library(stringr)
library(pryr)
library(OncoSimulR)
library(testthat)
library(Matrix)

library(plot.matrix) ## for the plot of the theta matrix


source("schill-trans-mat.R") ## yes, a few seconds because of the testing

## source("LOD-POM.R", echo = FALSE, max.deparse.length = 0)
source("pre-process.R", echo = FALSE, max.deparse.length = 0)
## source("capri-caprese-process.R", echo = FALSE, max.deparse.length = 0)
source("ot-process.R", echo = FALSE, max.deparse.length = 0)
## source("dip-process.R", echo = FALSE, max.deparse.length = 0)  ## No longer using it
source("cbn-process.R", echo = FALSE, max.deparse.length = 0)
source("dbn-process.R", echo = FALSE, max.deparse.length = 0)
source("hesbcn-process.R", echo = FALSE, max.deparse.length = 0)

if(MCCBN_INSTALLED)
    source("mccbn-process.R", echo = FALSE, max.deparse.length = 0)


registerDoSEQ() ## for TRONCO

boot_data_index <- function(x, boot) {
    ## Used by CBN and DiP
    ## boot is an integer. 0 means no boot
    ## that is because I reuse boot for two purposes
    boot <- as.logical(boot)
    if(boot) {
        ind <- sample(nrow(x), nrow(x), replace = TRUE)
        return(x[ind, , drop = FALSE])
    } else {
        return(x)
    }
}

add_pseudosamples <- function(x, n00 = "auto3") {
    if(n00 == "auto") {
        if(nrow(x) <= 500) {
            n00 <- round(nrow(x) * 0.10)
        } else {
            n00 <- round(nrow(x) * 0.05)
        }
    } else if(n00 == "auto2") {
        ## add only if max. freq. of any gene is > 95%
        fmax <- max(colSums(x))/nrow(x)
        if(fmax > 0.95)
            n00 <- round(nrow(x) * 0.05)
        else
            n00 <- 0
    } else if(n00 == "auto3") {
        ## add only if any gene is 100%
        ## add just 1
        fmax <- max(colSums(x))/nrow(x)
        if(fmax == 1) {
            cat("\n  Added one pseudosample \n ")
            n00 <- 1
        } else { 
            n00 <- 0
        }
    }
    return(rbind(x,
                 matrix(0L, nrow = n00, ncol = ncol(x))
                 ))
    ## cn <- colnames(x)
    
    ## tmp <- rbind(x,
    ##              matrix(0L, nrow = n00, ncol = ncol(x))
    ##              )
    ## colnames(tmp) <- cn
    ## return(tmp)
}

## NOTE: MCCBN allowed to run with arbitrary number of columns
all_methods <- function(x, nboot = 0, nboot_caprese_capri = 0,
                        nboot_cbn = 0, nboot_dip = 0,
                        n00 = "auto3", caprese_capri_minimal = TRUE,
                        caprese_capri_cores.ratio = 0,
                        distribution_oncotree = FALSE,
                        min.freq = 0,
                        cores_cbn = 1,
                        do_MCCBN = FALSE) { ## I think we want to keep min.freq to 0?
    
    if(caprese_capri_minimal) {
        ## Most of the time we only want the graphs, that's it
        evaleloss <- FALSE
        evalprederr <- FALSE
        evalposterr <- FALSE
    } else {
        evaleloss <- TRUE
        evalprederr <- TRUE
        evalposterr <- TRUE
    }
    x000 <- x
    x <- add_pseudosamples(x, n00 = n00)
    ## remove.constant makes no difference IFF we add pseudosamples, as
    ## there can be no constant column when we add pseudosamples
    x_tmp <- pre_process(x, remove.constant = FALSE, min.freq = min.freq)

    cat("\n     Doing OT")
    if(ncol(x_tmp) >= 2) {
        time_ot <- system.time(OT <- try(
                                   suppressMessages(
                                       ot_proc(x_tmp,
                                               nboot = nboot,
                                               distribution.oncotree = distribution_oncotree))))["elapsed"]
    } else {
        OT <- NA
        time_ot <- NA
    }

    ## ## remove.constant I think is here because if we do not add
    ## ## pseudosamples, caprese and capri will fail (fail, not just give a
    ## ## tree with a silly edge) if we pass n00 = 0 to add
    ## ## pseudosamples. And in some cases we have run this without adding
    ## ## pseudosamples. Not with the simulations, though.
    ## x_tmp <- pre_process(x, remove.constant = TRUE, min.freq = min.freq)
    ## if(ncol(x_tmp) >= 2) {
    ##     x_tmp <- tronco_common(x_tmp)
    ##     cat("\n     Doing CAPRESE")
    ##     time_caprese <- system.time(CAPRESE <- try(suppressMessages(
    ##                                     caprese_proc(genotData = x_tmp,
    ##                                 nboot = nboot_caprese_capri,
    ##                                 silent = TRUE,
    ##                                 cores.ratio = caprese_capri_cores.ratio,
    ##                                 evaleloss = evaleloss,
    ##                                 evalprederr = evalprederr,
    ##                                 evalposterr = evalposterr))))["elapsed"]
    ##     cat("\n     Doing CAPRI_BIC")
    ##     time_capri_bic <- system.time(CAPRI_BIC <- try(suppressMessages(
    ##                                   capri_proc(genotData = x_tmp,
    ##                                 nboot = nboot_caprese_capri,
    ##                                 nbootWilcox = 100,
    ##                                 regularization = "bic",
    ##                                 command = "hc",
    ##                                 pvalue = 0.05,
    ##                                 cores.ratio = caprese_capri_cores.ratio,
    ##                                 evaleloss = evaleloss,
    ##                                 evalprederr = evalprederr,
    ##                                 evalposterr = evalposterr))))["elapsed"]
    ##     cat("\n     Doing CAPRI_AIC")
    ##     time_capri_aic <- system.time(CAPRI_AIC <- try(suppressMessages(
    ##                                   capri_proc(genotData = x_tmp,
    ##                                 nboot = nboot_caprese_capri,
    ##                                 nbootWilcox = 100,
    ##                                 regularization = "aic",
    ##                                 command = "hc",
    ##                                 pvalue = 0.05,
    ##                                 cores.ratio = caprese_capri_cores.ratio,
    ##                                 evaleloss = evaleloss,
    ##                                 evalprederr = evalprederr,
    ##                                 evalposterr = evalposterr))))["elapsed"]
    ## } else {
    ##     CAPRESE <- CAPRI_AIC <- CAPRI_BIC <- NA
    ##     time_caprese <- time_capri_aic <- time_capri_bic <- NA
    ## }

    ## The code (ct-cbn.h) has
    ## 	if ((n < 1) || (n > 25))
    ## {
    ## 	fprintf(stderr, "Error:  Number of events is %d.  Supported range is {1, ..., 14}.\n", n);
    ## 	exit(1);
    ## }
    ## So despite the message, we can do up to 25?
    x_tmp <- pre_process(x, remove.constant = FALSE, min.freq = min.freq,
                         max.cols = 25)
    if(ncol(x_tmp) >= 2) {
        ## CBN_linear <- try(cbn_proc(x_tmp, addname = "tmpl",
        ##                            init.poset = "linear",
        ##                            nboot = nboot_cbn, parall = TRUE))
        ## CBN_linear <- NA ## disabled
        ## And note the parall argument has no effect as I am not using mclapply
        ## This was changed in cbn-process.R on commit b9027d5, on 2016-09-26.
        ## Yes, ugly!
        
        cat("\n     Doing CBN")
        time_cbn_ot <- system.time(CBN_ot <- try(cbn_proc(x_tmp, addname = "tmpo",
                                   init.poset = "OT",
                                   nboot = nboot_cbn, parall = TRUE,
                                   cores = cores_cbn)))["elapsed"]
    } else {
       time_cbn_ot <- NA
       CBN_ot <- NA 
    }
    if(do_MCCBN) {
        x_tmp <- pre_process(x, remove.constant = FALSE, min.freq = min.freq,
                             max.cols = NULL)
        cat("\n     Doing MCCBN")
        if(ncol(x_tmp) >= 2) {
            ## MCCBN should be able to run with many more columns. See max.cols above
            time_mccbn <- system.time(MCCBN <- try(mccbn_proc(x_tmp)))["elapsed"]
        } else {
            ## MCCBN <- CBN_linear <- CBN_ot <- NA
            time_mccbn <- NA
            MCCBN <- NA
        }
    } else {
        time_mccbn <- NA
        MCCBN <- NA  
    }
    ## x_tmp <- pre_process(x, remove.constant = FALSE, min.freq = min.freq,
    ##                      max.cols = 20)
    ## if(ncol(x_tmp) >= 2) {
    ##     ## if we are really going to use boottstrap,
    ##     ## it is more reasonable to parallelize over boot
    ##     ## than over p-range
    ##     if(nboot_dip == 0) {
    ##         ## When many runs from shell, we do not want any parall
    ##         ## cpr <- detectCores()
    ##         cpr <- 1
    ##         cpb <- 1
    ##         pb <- FALSE
    ##     } else {
    ##         cpr <- 1
    ##         cpb <- detectCores()
    ##         pb <- TRUE
    ##     }
    ##     time_dip_mpn <- system.time(DIP_MPN <- try(dip_proc(x_tmp, method = "MPN",
    ##                             addname = "dm",
    ##                             nboot = nboot_dip,
    ##                             cores_for_p_range = cpr,
    ##                             cores_for_boot = cpb,
    ##                             parall_boot = pb)))["elapsed"]
    ##     time_dip_smpn <- system.time(DIP_SMPN <- try(dip_proc(x_tmp, method = "SMPN", addname = "ds",
    ##                              nboot = nboot_dip,
    ##                              cores_for_p_range = cpr,
    ##                              cores_for_boot = cpb,
    ##                              parall_boot = pb)))["elapsed"]
    ## } else {
    ##     DIP_MPN <- DIP_SMPN <- NA
    ##     time_dip_mpn <- time_dip_smpn <- NA
    ## }

    cat("\n                  _times_cpm: ot", time_ot,
        " cbn ", time_cbn_ot,
        " mccbn ", time_mccbn,
        ## " caprese ", time_caprese,
        ## " capri_aic ", time_capri_aic,
        ## " capri_bic ", time_capri_bic,
        "\n")
    out <- list(
        OT = OT,
        ## CAPRESE = CAPRESE,
        ## CAPRI_BIC = CAPRI_BIC,
        ## CAPRI_AIC = CAPRI_AIC,        
        ## DIP_MPN = DIP_MPN,
        ## DIP_SMPN = DIP_SMPN,
        CBN_ot = CBN_ot,
        MCCBN = MCCBN,
        time_ot = time_ot,
        ## time_caprese = time_caprese,
        ## time_capri_bic = time_capri_bic,
        ## time_capri_aic = time_capri_aic,        
        ## time_dip_mpn = time_dip_mpn,
        ## time_dip_smpn = time_dip_smpn,
        time_cbn_ot = time_cbn_ot,
        time_mccbn = time_mccbn,
        input_data = x000, ## yes, return this!!!
        input_data_pseudosamples = x
    )
    return(out)
}




## This is coming from
## Cancer_Data_sets/CPM-weighted-paths-biol.R
## in the supplementary material for Diaz-Uriarte and Vasallo.
## We leave in there things we don't really need. Simpler.


## DAG of restrictions (as data frame) -> vector of accessible genotypes and graph of DAG of restrictions
## Under an AND model, such as CBN
## return all the accessible genotypes
##     from a DAG of genes plus the DAG as igraph object
df_2_access_genots_and_graph <- function(x) {
    
    ## minor detail: if x is a sparse adjacency mat.
    ##    from igraph this still works. But don't do that.
    g <- igraph::graph_from_data_frame(x, directed = TRUE)
    
    ## g: an igraph object, with a "Root"
    ## returns a list
    children <- sort(setdiff(V(g)$name, "Root"))
    node_depth <- unlist(lapply(children,
                                function(node)
                                    max(unlist(lapply(all_simple_paths(g,
                                                                       from = "Root",
                                                                       to = node),
                                                      length)))
                                ))
    
    names(node_depth) <- children
    node_depth <- sort(node_depth)
    ## pre-allocate a list.
    ## FIXME: Could be smarter as a function of dim(x)?
    all_gty <- vector(mode = "list", length = 100) 
    i <- 1
    for(j in seq_along(node_depth)) {
         
        tmp_gty_1 <- sort(setdiff(names(subcomponent(g, v = names(node_depth)[j],
                                                     mode = "in")),
                                  "Root"))
        all_gty[[i]] <- tmp_gty_1
        
        ## only do union of those not contained in the genotype
        ## to_unite <-  which(unlist(lapply(all_gty[1:(i-1)],
        ##                                  function(z) length(setdiff(z, tmp_gty_1)) > 0)))
        
        if(i > 1) {
            to_unite <-  which(unlist(lapply(all_gty[1:(i-1)],
                                             function(z) !all(z %in% tmp_gty_1))))
        } else {
            to_unite <- vector(length = 0)
        }
        
        if(length(to_unite)) {
            ## we need unique as some sets are the same
            tmp_gty_2 <- unique(lapply(all_gty[to_unite],
                                       function(u) sort(union(u, tmp_gty_1))))
            ## check remaining space in preallocated list.
            ## expand if needed.
            if(length(all_gty) < (i + 1 + length(tmp_gty_2))) {
                all_gty <- c(all_gty,
                             vector(mode = "list", length = max(length(all_gty),
                                                                2 + length(tmp_gty_2))))
            }
            all_gty[(i + 1):(i + length(tmp_gty_2))] <- tmp_gty_2
            i <- i + length(tmp_gty_2)
        }
        i <- i + 1
    }
    all_gty <- all_gty[1:(i - 1)]
    ng <- unlist(lapply(all_gty, length))
    all_gty <- all_gty[order(ng)]
    return(list(accessible_genots = all_gty,
                graph = g))
}






## DAG of restrictions (as data frame) -> vector of accessible genotypes and graph of DAG of restrictions
##  Under an OR model (not XOR, not AND), such as DBN
## return all the accessible genotypes
##     from a DAG of genes plus the DAG as igraph object
df_2_access_genots_and_graph_OR <- function(x) {
    
    ## minor detail: if x is a sparse adjacency mat.
    ##    from igraph this still works. But don't do that.
    g <- igraph::graph_from_data_frame(x, directed = TRUE)
    
    ## g: an igraph object, with a "Root"
    ## returns a list
    children <- sort(setdiff(V(g)$name, "Root"))
    node_depth <- unlist(lapply(children,
                                function(node)
                                    max(unlist(lapply(all_simple_paths(g,
                                                                       from = "Root",
                                                                       to = node),
                                                      length)))
                                ))
    
    names(node_depth) <- children
    node_depth <- sort(node_depth)
    ## pre-allocate a list.
    all_gty <- vector(mode = "list", length = 100) 
    i <- 1
    for(j in seq_along(node_depth)) {
        ## FIXME: we did this traversing above. Reuse that
        all_tmp_gty_1 <- lapply(
            all_simple_paths(g, from = "Root", to = names(node_depth)[j],
                             mode = "out"),
            names)

        all_tmp_gty_1 <- lapply(all_tmp_gty_1,
                                function(u) sort(setdiff(u, "Root")))
        for(k in seq_along(all_tmp_gty_1)) {
            tmp_gty_1 <- all_tmp_gty_1[[k]]
            all_gty[[i]] <- tmp_gty_1 
            
            
            if(i > 1) {
                to_unite <-  which(unlist(lapply(all_gty[1:(i-1)],
                                                 function(z) !all(z %in% tmp_gty_1))))
            } else {
                to_unite <- vector(length = 0)
            }
            
            if(length(to_unite)) {
                ## we need unique as some sets are the same
                tmp_gty_2 <- unique(lapply(all_gty[to_unite],
                                           function(u) sort(union(u, tmp_gty_1))))
                ## check remaining space in preallocated list.
                ## expand if needed.
                if(length(all_gty) < (i + 1 + length(tmp_gty_2))) {
                    all_gty <- c(all_gty,
                                 vector(mode = "list", length = max(length(all_gty),
                                                                    2 + length(tmp_gty_2))))
                }
                all_gty[(i + 1):(i + length(tmp_gty_2))] <- tmp_gty_2
                i <- i + length(tmp_gty_2)
            }
            i <- i + 1
        }
    }
    all_gty <- all_gty[1:(i - 1)]
    ## FIXME: We got duplicates above. Easy to rm here, but better if we did not
    ## get them to begin with
    all_gty <- unique(all_gty)
    ng <- unlist(lapply(all_gty, length))
    all_gty <- all_gty[order(ng)]
    return(list(accessible_genots = all_gty,
                graph = g))
}

df_2_access_genots_and_graph_XOR <- function(x) {
    ##PHN
    ## minor detail: if x is a sparse adjacency mat.
    ##    from igraph this still works. But don't do that.
    g <- igraph::graph_from_data_frame(x, directed = TRUE)
    
    ## g: an igraph object, with a "Root"
    ## returns a list
    # children <- sort(setdiff(V(g)$name, "Root"))
    # node_depth <- unlist(lapply(children,
    #                             function(node)
    #                                 max(unlist(lapply(all_simple_paths(g,
    #                                                                    from = "Root",
    #                                                                    to = node),
    #                                                   length)))
    #                             ))
    
    # names(node_depth) <- children
    # node_depth <- sort(node_depth)
    # ## pre-allocate a list.
    # ## FIXME: Could be smarter as a function of dim(x)?
    # all_gty <- vector(mode = "list", length = 100) 
    # i <- 1
    # for(j in seq_along(node_depth)) {
         
    #     tmp_gty_1 <- sort(setdiff(names(subcomponent(g, v = names(node_depth)[j],
    #                                                  mode = "in")),
    #                               "Root"))
    #     all_gty[[i]] <- tmp_gty_1
        
    #     ## only do union of those not contained in the genotype
    #     ## to_unite <-  which(unlist(lapply(all_gty[1:(i-1)],
    #     ##                                  function(z) length(setdiff(z, tmp_gty_1)) > 0)))
        
    #     if(i > 1) {
    #         to_unite <-  which(unlist(lapply(all_gty[1:(i-1)],
    #                                          function(z) !all(z %in% tmp_gty_1))))
    #     } else {
    #         to_unite <- vector(length = 0)
    #     }
        
    #     if(length(to_unite)) {
    #         ## we need unique as some sets are the same
    #         tmp_gty_2 <- unique(lapply(all_gty[to_unite],
    #                                    function(u) sort(union(u, tmp_gty_1))))
    #         ## check remaining space in preallocated list.
    #         ## expand if needed.
    #         if(length(all_gty) < (i + 1 + length(tmp_gty_2))) {
    #             all_gty <- c(all_gty,
    #                          vector(mode = "list", length = max(length(all_gty),
    #                                                             2 + length(tmp_gty_2))))
    #         }
    #         all_gty[(i + 1):(i + length(tmp_gty_2))] <- tmp_gty_2
    #         i <- i + length(tmp_gty_2)
    #     }
    #     i <- i + 1
    # }
    # all_gty <- all_gty[1:(i - 1)]
    # ng <- unlist(lapply(all_gty, length))
    # all_gty <- all_gty[order(ng)]
    return(list(accessible_genots = all_gty,
                graph = g))
}



## ### Examples
## x2 <- data.frame(From = c("Root", "Root", "A", "B", "C"),
##                  To = c("A", "B", "C", "C", "D"))

## x3 <- data.frame(From = c("Root", "Root", "A", "D", "B", "E"),
##                  To = c("A", "D", "B", "E", "C", "C"))


## df_2_access_genots_and_graph_OR(x2)
## df_2_access_genots_and_graph_OR(x3)




## From function of same name in ruggify-functions.R

## vector of accessible genotypes -> adjacency matrix of genotypes (fitness graph)
## return maximally connected fitness graph for a given set of accessible genotypes
unrestricted_fitness_graph <- function(gacc, plot = FALSE) {
    
    gs <- unlist(lapply(gacc, function(g) paste0(g, collapse = ", ")))
    gs <- c("WT", gs)
    nmut <- c(0, vapply(gacc, length, 1))
    
    adjmat <- matrix(0L, nrow = length(gs), ncol = length(gs))
    rownames(adjmat) <- colnames(adjmat) <- gs
    
    adjmat["WT", gs[which(nmut == 1)]] <- 1L
    
    for(m in 2:max(nmut)){
        g <- gs[which(nmut == m)]
        for (gn in g) {
            parents <- gacc[which(nmut == m-1)-1]
            gns <- unlist(strsplit(gn, ", "))
            parents <- parents[which(unlist(lapply(parents,
                                                   function(p)
                                                       length(setdiff(gns, p)))) == 1)]
            for (p in parents){
                adjmat[paste0(p, collapse = ", "), gn] <- 1L
            }
        }
    }
    if (plot)
        mccbn::plot_poset(adjmat) ## , title = "G0 (unrestricted)")

    stopifnot(all(adjmat %in% c(0L, 1L) ))
    storage.mode(adjmat) <- "integer"

    return(adjmat)
}



## Trying sparse matrices
## https://stackoverflow.com/questions/23107837/r-sparse-matrix-from-list-of-dimension-names
## https://stackoverflow.com/questions/26207850/create-sparse-matrix-from-a-data-frame
## https://cmdlinetips.com/2019/05/introduction-to-sparse-matrices-in-r/
## https://www.gormanalysis.com/blog/sparse-matrix-construction-and-use-in-r/

## vector of accessible genotypes -> adjacency matrix of genotypes (fitness graph)
## return maximally connected fitness graph for a given set of accessible genotypes
unrestricted_fitness_graph_sparseM <- function(gacc, plot = FALSE) {
    gs <- unlist(lapply(gacc, function(g) paste0(g, collapse = ", ")))
    gs <- c("WT", gs)
    nmut <- c(0, vapply(gacc, length, 1))
    
    ## adjmat <- matrix(0L, nrow = length(gs), ncol = length(gs))
    ## rownames(adjmat) <- colnames(adjmat) <- gs
    
    ## ## adjmat <- Matrix(0L, nrow = length(gs), ncol = length(gs),
    ## ##                  sparse = TRUE, dimnames = list(gs, gs))
    ## ## This works but I don't feel comfortable
    ## ## adjmat["WT", gs[which(nmut == 1)]] <- 1L
    ## adjmat[i = rep(ii, length(jj)), j = jj] <- 1L

    jj <- match(gs[which(nmut == 1)], gs)
    ii <- rep.int(match("WT", gs), length(jj))
    adjmat <- sparseMatrix(i = ii, j = jj, x = 1L,
                           dims = c(length(gs), length(gs)), 
                           dimnames = list(gs, gs))

    for(m in 2:max(nmut)){
        g <- gs[which(nmut == m)]
        for (gn in g) {
            parents <- gacc[which(nmut == m-1)-1]
            gns <- unlist(strsplit(gn, ", "))
            parents <- parents[which(unlist(lapply(parents,
                                                   function(p)
                                                       length(setdiff(gns, p)))) == 1)]
            for (p in parents){
                ## Works but better via indices, I think
                ## adjmat[paste0(p, collapse = ", "), gn] <- 1L
                jjj <- match(gn, gs)
                iii <- rep.int(match(paste0(p, collapse = ", "), gs),
                               length(jjj))
                adjmat[iii, jjj] <- 1L
            }
        }
    }
    ## Probably will not work
    if (plot)
        mccbn::plot_poset(adjmat) ## , title = "G0 (unrestricted)")
    ## stopifnot(all(adjmat %in% c(0L, 1L) ))
    ## storage.mode(adjmat) <- "integer"
    return(adjmat)
}



## list of accessible genotypes -> global maximum
##    Beware: just a single maximum
get_global_max <- function(ag) {
    muts <- unlist(lapply(ag, length))
    max_muts <- max(muts)
    ind_max_muts <- which(muts == max_muts)
    if(length(ind_max_muts) != 1) stop("eh??!! ind_max_muts")
    return(paste0(ag[[ind_max_muts]], collapse = ", "))
}


## < /from CPMs-paths-genotypes-and-comb.R >

## output of CPM analysis, string -> accessible genotypes and paths
##                  string: just to identify errors
##    the _w: weights, so we add probs.
##   Modification of function of same name, without _w in
##   CPMs-paths-genotypes-and-comb.R

cpm_access_genots_paths_w <- function(x, string = NULL,
                                    names_weights_paths =
                                        c("rerun_lambda",
                                          "lambda",
                                          "OT_edgeWeight")) {
    if(inherits(x, "try-error") || is.na(x) || is.null(x)) {
        ## The CPM analysis produced no edges component, so
        ## nothing can be done
        if(inherits(x, "try-error")) likely_error <- "Error_in_run"
        if(is.na(x)) likely_error <- "ncol_x"
        if(is.null(x)) likely_error <- "other_error"
        return(list(accessible_genots = "ERROR_CPM_ANALYSIS",
                    num_accessible_genots = "ERROR_CPM_ANALYSIS",
                    CPM_DAG_as_igraph = "ERROR_CPM_ANALYSIS",
                    fgraph = "ERROR_CPM_ANALYSIS",
                    num_paths_to_max = "ERROR_CPM_ANALYSIS",
                    paths_error = "ERROR_CPM_ANALYSIS",
                    paths_to_max = "ERROR_CPM_ANALYSIS",
                    weighted_fgraph = "ERROR_CPM_ANALYSIS",
                    trans_mat_genots = "ERROR_CPM_ANALYSIS",
                    unweighted_paths_to_max = "ERROR_CPM_ANALYSIS",
                    weighted_paths_to_max = "ERROR_CPM_ANALYSIS",
                    diversity_weighted_paths_to_max = "ERROR_CPM_ANALYSIS",
                    likely_error = likely_error
                    
                ))
    }
                         
    x <- x$edges
    tmp <- try(df_2_access_genots_and_graph(x[, c("From", "To")]))
    
    if(inherits(tmp, "try-error")) {
        stop("how is this happening? there was edges component!")
    } else {
        accessible_genots <- tmp$accessible_genots
        fgraph <- unrestricted_fitness_graph(accessible_genots)
        fgraphi <- igraph::graph_from_adjacency_matrix(fgraph)
        gmax <- get_global_max(accessible_genots)
    }

    ## based on n_pahts_single_max in compute-numpaths-clonal-int-stats.R
    ## I still need this to get the paths from the CPM
    lpaths <- NA
    paths_to_max <- NA
    paths_error <- FALSE
    paths <- try(igraph::all_simple_paths(fgraphi,
                                          from = "WT",
                                          to = gmax,
                                          mode = "out"))
    if(inherits(paths, "try-error")) {
        cat("\n     ERROR_in_paths_in_calling_string = ", string, "\n")
        lpaths <- -99
        paths_error <- TRUE
        paths_to_max <- "PATHS_ERROR" ## this will never appear anywhere FIXME
        ## will have to look for the paths_error variable
    } else {
        lpaths <- length(paths)
        paths_to_max <-  unlist(lapply(paths,
                                       function(x) paste(igraph::as_ids(x),
                                                         collapse = " -> ")))
    }

    
    ## Logic of obtaining weighted fitness graph
    
    ##  - obtain all accessible genotypes and fitness graph (unrestricted
    ##      fitness graph) from CPM

    ##  - from the lambdas/probs. of genes given genes obtain
    ##       probabilities/lambdas of genotypes given genotypes. The call
    ##       to "transition_fg" (where weights is the conditional
    ##       prob./lambda) of descendant gene given parent gene

    ##  - when using CBN, might not be probabilities. Make sure they are
    ##    transition probabilities between genotypes: sweep

    ##  - find the probability of each path: call to
    ##    do_weighted_paths_to_max, that uses the list of all paths_to_max. 

    
    which_col_weights <- which(colnames(x) %in% names_weights_paths)
    if(length(which_col_weights) > 1) {
        stop("more than one column with weights")
    } else if (length(which_col_weights) == 1) {
        stopifnot(colnames(x)[2] == "To")
        weights <- unique(x[, c(2, which_col_weights)])
        if(any(duplicated(weights[, "To"]))) {
            stop("Different lambda/weight for same destination gene.",
                 " This is not allowed with conjunctive DAGs")
        }
        rownames(weights) <- weights[, "To"]
        weighted_fgraph <- transition_fg(fgraph, weights)
    } else {
        ## why would we return something? It is NA
        ## weighted_fgraph <- fgraph
        weighted_fgraph <- NA
    }

    if (length(which_col_weights) == 1) {
        ## weighted_fgraph need not have each row sum to 1, if they are lambdas
        ## from CBN for instance. So make sure they are transition matrices
        ## between genotypes.
        trans_mat_genots <- sweep(weighted_fgraph, 1,
                                  rowSums(weighted_fgraph), FUN = "/")
        trans_mat_genots[is.nan(trans_mat_genots)] <- 0
    } else {
        trans_mat_genots <- NA
    }


    if(length(which_col_weights) == 1) {
        weighted_paths_to_max <- do_weighted_paths_to_max(paths_to_max,
                                                          trans_mat_genots)
        diversity_weighted_paths_to_max <-
            OncoSimulR:::shannonI(weighted_paths_to_max[, "probability"])

    } else {
        weighted_paths_to_max <- NA
        diversity_weighted_paths_to_max <- NA
    }
    
    return(list(accessible_genots = accessible_genots,
                num_accessible_genots = length(accessible_genots),
                CPM_DAG_as_igraph = tmp$graph,
                ## fgraph_AM = fgraph,
                fgraph = fgraph,
                num_paths_to_max = lpaths,
                paths_error = paths_error,
                weighted_fgraph = weighted_fgraph,
                trans_mat_genots = trans_mat_genots,
                unweighted_paths_to_max = paths_to_max,
                weighted_paths_to_max = weighted_paths_to_max,
                diversity_weighted_paths_to_max =
                    diversity_weighted_paths_to_max,
                likely_error = "No_error"
                ))
}





## < /from CPMs-paths-genotypes-and-comb.R >

## output of CPM analysis, string -> accessible genotypes and paths
##                  string: just to identify errors
##    the _w: weights, so we add probs.
##   Modification of function of same name, without _w in
##   CPMs-paths-genotypes-and-comb.R

## Like the one above, but only with necessary output
##  for both speed and size and using sparse matrices.
cpm_access_genots_paths_w_simplified <- function(x, string = NULL,
                                    names_weights_paths =
                                        c("rerun_lambda",
                                          "lambda",
                                          "OT_edgeWeight")) {
    if(inherits(x, "try-error") || is.na(x) || is.null(x)) {
        ## The CPM analysis produced no edges component, so
        ## nothing can be done
        ## if(inherits(x, "try-error")) likely_error <- "Error_in_run"
        ## if(is.na(x)) likely_error <- "ncol_x"
        ## if(is.null(x)) likely_error <- "other_error"
        return(list(## accessible_genots = "ERROR_CPM_ANALYSIS",
                    ## num_accessible_genots = "ERROR_CPM_ANALYSIS",
                    ## CPM_DAG_as_igraph = "ERROR_CPM_ANALYSIS",
                    fgraph = "ERROR_CPM_ANALYSIS",
                    ## num_paths_to_max = "ERROR_CPM_ANALYSIS",
                    ## paths_error = "ERROR_CPM_ANALYSIS",
                    ## paths_to_max = "ERROR_CPM_ANALYSIS",
                    weighted_fgraph = "ERROR_CPM_ANALYSIS",
                    trans_mat_genots = "ERROR_CPM_ANALYSIS"
                    ## unweighted_paths_to_max = "ERROR_CPM_ANALYSIS",
                    ## weighted_paths_to_max = "ERROR_CPM_ANALYSIS",
                    ## diversity_weighted_paths_to_max = "ERROR_CPM_ANALYSIS",
                    ## likely_error = likely_error
                    
                ))
    }
                         
    x <- x$edges
    tmp <- try(df_2_access_genots_and_graph(x[, c("From", "To")]))
     
    if(inherits(tmp, "try-error")) {
        stop("how is this happening? there was edges component!")
    } else {
        accessible_genots <- tmp$accessible_genots
        fgraph <- unrestricted_fitness_graph_sparseM(accessible_genots)
        ## fgraphi <- igraph::graph_from_adjacency_matrix(fgraph)
        ## gmax <- get_global_max(accessible_genots)
    }

    
    ## ## based on n_pahts_single_max in compute-numpaths-clonal-int-stats.R
    ## ## I still need this to get the paths from the CPM
    ## lpaths <- NA
    ## paths_to_max <- NA
    ## paths_error <- FALSE
    ## paths <- try(igraph::all_simple_paths(fgraphi,
    ##                                       from = "WT",
    ##                                       to = gmax,
    ##                                       mode = "out"))
    ## ## if(inherits(paths, "try-error")) {
    ## ##     cat("\n     ERROR_in_paths_in_calling_string = ", string, "\n")
    ## ##     lpaths <- -99
    ## ##     paths_error <- TRUE
    ## ##     paths_to_max <- "PATHS_ERROR" ## this will never appear anywhere FIXME
    ## ##     ## will have to look for the paths_error variable
    ## ## } else {
    ## ##     lpaths <- length(paths)
    ## ##     paths_to_max <-  unlist(lapply(paths,
    ## ##                                    function(x) paste(igraph::as_ids(x),
    ## ##                                                      collapse = " -> ")))
    ## ## }

    
    ## Logic of obtaining weighted fitness graph
    
    ##  - obtain all accessible genotypes and fitness graph (unrestricted
    ##      fitness graph) from CPM

    ##  - from the lambdas/probs. of genes given genes obtain
    ##       probabilities/lambdas of genotypes given genotypes. The call
    ##       to "transition_fg" (where weights is the conditional
    ##       prob./lambda) of descendant gene given parent gene

    ##  - when using CBN, might not be probabilities. Make sure they are
    ##    transition probabilities between genotypes: sweep

    ##  - find the probability of each path: call to
    ##    do_weighted_paths_to_max, that uses the list of all paths_to_max. 

    which_col_weights <- which(colnames(x) %in% names_weights_paths)
    if(length(which_col_weights) > 1) {
        stop("more than one column with weights")
    } else if (length(which_col_weights) == 1) {
        stopifnot(colnames(x)[2] == "To")
        weights <- unique(x[, c(2, which_col_weights)])
        if(any(duplicated(weights[, "To"]))) {
            stop("Different lambda/weight for same destination gene.",
                 " This is not allowed with conjunctive DAGs")
        }
        rownames(weights) <- weights[, "To"]
        weighted_fgraph <- transition_fg_sparseM(fgraph, weights)

    } else {
        ## why would we return something? It is NA
        ## weighted_fgraph <- fgraph
        weighted_fgraph <- NA
    }

    if (length(which_col_weights) == 1) {
        ## weighted_fgraph need not have each row sum to 1, if they are lambdas
        ## from CBN for instance. So make sure they are transition matrices
        ## between genotypes.
        
        ## trans_mat_genots <- sweep(weighted_fgraph, 1,
        ##                           rowSums(weighted_fgraph), FUN = "/")

        trans_mat_genots <- rowScaleMatrix(weighted_fgraph)
    } else {
        trans_mat_genots <- NA
    }


    ## if(length(which_col_weights) == 1) {
    ##     weighted_paths_to_max <- do_weighted_paths_to_max(paths_to_max,
    ##                                                       trans_mat_genots)
    ##     diversity_weighted_paths_to_max <-
    ##         OncoSimulR:::shannonI(weighted_paths_to_max[, "probability"])

    ## } else {
    ##     weighted_paths_to_max <- NA
    ##     diversity_weighted_paths_to_max <- NA
    ## }
    
    return(list(
        ## accessible_genots = accessible_genots,
        ## num_accessible_genots = length(accessible_genots),
        ## CPM_DAG_as_igraph = tmp$graph,
        ## fgraph_AM = fgraph,
        fgraph = fgraph,
        ## num_paths_to_max = lpaths,
        ## paths_error = paths_error,
        weighted_fgraph = weighted_fgraph,
        trans_mat_genots = trans_mat_genots
        ## unweighted_paths_to_max = paths_to_max,
        ## weighted_paths_to_max = weighted_paths_to_max,
        ## diversity_weighted_paths_to_max =  diversity_weighted_paths_to_max,
        ## likely_error = "No_error"
    ))
}







## adjacency matrix genotypes, row, column of that matrix,
##      lambdas/weights of descendant gene given parent gene ->
##     the lambda of descendant genotype
## Return the lambda/weight of a single genotype coming from a single genotype
get_single_lambda <- function(x, row, col, weights) {
    genot_to <- unlist(strsplit(colnames(x)[col], ", ", fixed = TRUE))
    genot_from <- unlist(strsplit(rownames(x)[row], ", ", fixed = TRUE))
    ## WT stuff
    if((length(genot_from) == 1) && (genot_from == "WT")) {
        genot_from == ""
    }
    added_gene <- setdiff(genot_to, genot_from)
    return(weights[added_gene, 2])
}



## fitness graph, weights of probs/lambdas descendant gene given parent gene -> weighted fitness graph
##             the fitness graph with weights (relative to each node)
##             of jumping to each of the descendant genotypes

## ## No longer used. Turning into transition matrix done outside
## ## 
## ##          transition = TRUE: return transition matrix
## ##          if they sumto zero by row this is really the transition matrix
## ##          of genotypes
transition_fg <- function(x, weights) { ## , transition = TRUE) {
    pos_do <- which(x == 1, arr.ind = TRUE)
    wfg <- x
    wfg[] <- 0
    tmp <- unlist(
        Map(function(r, c)
             get_single_lambda(x, r, c, weights),
            pos_do[, 1], pos_do[, 2]))
    wfg[pos_do] <- tmp
        
    ## ## Could use Map but let's use a loop instead
    ## for(p in 1:nrow(pos_do)) {
    ##     wfg[pos_do[p, ] ] <-  get_single_lambda(x,
    ##                                             pos_do[p, 1], pos_do[p, 2],
    ##                                             weights)
    ## }
    return(wfg)
}


## fitness graph, weights of probs/lambdas descendant gene given parent gene -> weighted fitness graph
##             the fitness graph with weights (relative to each node)
##             of jumping to each of the descendant genotypes

## ## No longer used. Turning into transition matrix done outside
## ## 
## ##          transition = TRUE: return transition matrix
## ##          if they sumto zero by row this is really the transition matrix
## ##          of genotypes
transition_fg_sparseM <- function(x, weights) { ## , transition = TRUE) {
    ## pos_do <- which(x == 1, arr.ind = TRUE)
    pos_do <- as.matrix(summary(x)[, c("i", "j")])
    wfg <- x
    wfg[] <- 0
    tmp <- unlist(
        Map(function(r, c)
             get_single_lambda(x, r, c, weights),
            pos_do[, 1], pos_do[, 2]))
    wfg[pos_do] <- tmp
    ## ## Could use Map but let's use a loop instead
    ## for(p in 1:nrow(pos_do)) {
    ##     wfg[pos_do[p, ] ] <-  get_single_lambda(x,
    ##                                             pos_do[p, 1], pos_do[p, 2],
    ##                                             weights)
    ## }
    return(wfg)
}

## vector of paths to max, transition matrix genotypes ->
##     weighted paths to max
do_weighted_paths_to_max <- function(paths, trans_mat) {
    probs <- vapply(paths, function(u) prob_single_path(u, trans_mat),
                    FUN.VALUE = 0.0)
    stopifnot(isTRUE(all.equal(sum(probs), 1)))
    df <- data.frame(path = paths,
                     probability = probs, stringsAsFactors = FALSE)

    rownames(df) <- NULL
    return(df)
}

## path to max, transition matrix genotypes ->
##    probability of path
prob_single_path <- function(path, trans_mat) {
    ## the indices of the sequence of genotypes
    ii <- which(row.names(trans_mat) %in%
                strsplit(path, " -> ", fixed = TRUE)[[1]])
    prod(trans_mat[cbind(ii[-length(ii)], ii[-1])])
}



## ## file with a single entry of cpm output -> accessible genotypes and paths
## ##    file is an ANALYZED__ID* file, as produced by run-CPMs.R
## cpm_out_to_paths_genots_w <- function(thefile,
##                                       methods = c("CBN", "MCCBN")
##                                       ## methods = c("OT", "CAPRESE",
##                                       ##               "CAPRI_BIC", "CAPRI_AIC",
##                                       ##           "CBN_ot", "MCCBN")
##                                       ) {
##     out <- NULL
##     load(thefile)
##     ## this creates the out, that codetools complaints about
##     outn <- out[1:5]
##     if(is.null(out$bootstrap)) {
##         outn$bootstap <- FALSE
##     } else {
##         outn$bootstrap <- out$bootstrap
##     }

##     string0 <- as.data.frame(outn, stringsAsFactors = FALSE)

##     string <- paste(paste(names(string0), string0, sep = "_"), collapse = "__")

##     out_w <- cpm_access_genots_paths_w(out$cbn_out, string = string)

##     theout <- c(outn,
##                 out_paths_genots = list(out_w)
##                 )
    
##     ## saveRDS(theout, file = paste0("paths_genots_index_", index, "_", string, ".rds"))
##     ## naming more consistent with that of true fitness graphs
##     saveRDS(theout, file = paste0("paths_", string, ".rds"))
##     cat("\n       done string = ", string, "\n")

##     likely_error <- out_w$likely_error
##     paths_error <- out_w$paths_error
    
##     m1 <- data.frame(num_paths_to_max = out_w$num_paths_to_max,
##                      num_accessible_genots = out_w$num_accessible_genots,
##                      diversity_weighted_paths_to_max = out_w$diversity_weighted_paths_to_max)
                     
##     m1$replicate <- outn$iter
##     likely_error <- as.data.frame(likely_error)
##     paths_error <- as.data.frame(paths_error)
    
##     m1 <- cbind(m1, likely_error)
##     m1 <- cbind(m1, paths_error)
    
##     rm(theout)
##     dfoutn <- as.data.frame(outn)
##     row.names(dfoutn) <- NULL
    
##     retout <- cbind(dfoutn, m1)
##     return(retout)
## }




## a simple check
any_constant_col <- function(x) {
    nr <- nrow(x)
    mcs <- max(colSums(x))
    any(mcs == nr)
}

## convert data frame to a matrix with 0L and 1L
df_2_mat_integer <- function(x) {
    x1 <- as.matrix(x)
    if(max(abs(x - x1)) != 0) stop("failed conversion to matrix")
    x2 <- x1
    storage.mode(x2) <- "integer"
    if(max(abs(x2 - x1)) != 0) stop("Not in 0L, 1L")
    if(max(abs(x - x2)) != 0) stop("df not in 0L, 1L") ## paranoia
    return(x2)
}

## To make it explicit
## but do not set the last row to NaNs or similar.
## Following same logic as in trans_rate_to_trans_mat (in MHN dir)
rowScaleMatrix <- function(x) {
    tm <- x
    sx <- rowSums(x)
    ii <- which(sx > 0)
    for(i in ii) {
        tm[i, ] <- tm[i, ]/sx[i]
    }
    tm
}


## Pass a data set as a matrix with subjects as rows and genes as columns

all_methods_2_trans_mat <- function(x, cores_cbn = 1, do_MCCBN = FALSE) {
     
    x <- df_2_mat_integer(x)

    ## cat("\n  Number of genes before limiting = ", ncol(x))
    ## Always limit to 15
    x <- pre_process(x, remove.constant = FALSE,
                     min.freq = 0, max.cols = 15)
    
    ## cat("\n  Number of genes after limiting = ", ncol(x), "\n")

    if(do_MCCBN)
        methods <- c("OT", ##"CAPRESE", "CAPRI_BIC", "CAPRI_AIC",
                     "CBN_ot"
                   , "MCCBN")
    else
        methods <- c("OT", "CBN_ot")

    cat("\n     Doing MHN")
    time_schill <- system.time(
        out_schill <- do_MHN2(x, lambda = 1/nrow(x)))["elapsed"]

    cat("\n  time MHN = ", time_schill)
    
    cat("\n     Doing DBN")
    time_dbn <- system.time(
      out_dbn <- do_DBN(x))["elapsed"]
    
    cat("\n  time DBN = ", time_dbn)

    cat("\n     Doing HESBCN")
    time_hesbcn <- system.time(
      out_hesbcn <- do_HESBCN(x))["elapsed"]
    
    cat("\n  time HESBCN = ", time_hesbcn)
     
    cpm_out_others <- all_methods(x, cores_cbn = cores_cbn, do_MCCBN = do_MCCBN)
    pre_trans_mat_others <- lapply(cpm_out_others[methods],
        cpm_access_genots_paths_w_simplified)

    pre_trans_mat_DBN_HESBCN <- mapply(
        cpm_access_genots_paths_w_simplified_OR, 
        list(out_dbn, out_hesbcn),
        c("Thetas", "Lambdas"))
    cat("\n    getting transition matrices for all non-mhn methods \n")

    ## ## Unweighted
    ## uw <- lapply(pre_trans_mat_others, function(x) rowScaleMatrix(x$fgraph))
    ## Weighted
    
    wg <- lapply(pre_trans_mat_others[c("OT", "MCCBN" , "CBN_ot")[c(TRUE, do_MCCBN, TRUE)]],
                 function(x) x$trans_mat_genots)
    ## Diagonal
    browser()
    td <- lapply(pre_trans_mat_others[c("MCCBN", "CBN_ot")[c(do_MCCBN, TRUE)]],
                 function(x) trans_rate_to_trans_mat(x$weighted_fgraph,
                                                     method = "uniformization"))
    ## ## Paranoid check
    ## wg2 <- lapply(pre_trans_mat_others[c("OT", "MCCBN", "CBN_ot")],
    ##               function(x) trans_rate_to_trans_mat(x$weighted_fgraph,
    ##                                                   method = "competingExponentials"))
    ## stopifnot(identical(wg, wg2))

    if(do_MCCBN) {
        MCCBN_model <- cpm_out_others$MCCBN$edges
        MCCBN_f_graph <- pre_trans_mat_others$MCCBN$weighted_fgraph
        MCCBN_trans_mat <- wg$MCCBN
        MCCBN_td_trans_mat <- td$MCCBN
    } else {
        MCCBN_model <- NA
        MCCBN_f_graph <- NA
        MCCBN_trans_mat <- NA
        MCCBN_td_trans_mat <- NA
    }
    ## TODO: return DBN results
    return(list(
        OT_model = cpm_out_others$OT$edges,
        OT_f_graph = pre_trans_mat_others$OT$weighted_fgraph,
        OT_trans_mat = wg$OT,
        ## OT_u = uw$OT,
        CBN_model = cpm_out_others$CBN_ot$edges,
        CBN_f_graph = pre_trans_mat_others$CBN_ot$weighted_fgraph,
        CBN_trans_mat = wg$CBN_ot,
        CBN_td_trans_mat = td$CBN_ot,
        ## CBN_uw = uw$CBN_ot,
        MCCBN_model = MCCBN_model,
        MCCBN_f_graph = MCCBN_f_graph,
        MCCBN_trans_mat = MCCBN_trans_mat,
        MCCBN_td_trans_mat = MCCBN_td_trans_mat,
        ## MCCBN_uw = uw$MCCBN,
        ## CAPRESE = uw$CAPRESE,
        ## CAPRI_BIC = uw$CAPRI_BIC,
        ## CAPRI_AIC = uw$CAPRI_AIC,
        MHN_theta = out_schill$theta,
        MHN_exp_theta = exp(out_schill$theta),
        MHN_transitionRateMatrix = out_schill$transitionRateMatrix, 
        MHN_trans_mat = out_schill$transitionMatrixCompExp,
        MHN_td_trans_mat = out_schill$transitionMatrixTimeDiscretized,
        DBN_model = out_dbn$edges,
        DBN_likelihood = out_dbn$likelihood,
        DBN_f_graph = out_dbn$weighted_fgraph,
        DBN_trans_mat = out_dbn$trans_mat_genots
         ))
}




## use plot matrix to plot the sampled genotypes
plot_sampled_genots <- function(data) {
    d1 <- as.data.frame(sampledGenotypes(data))

    ## Reorder by number mutations and names, but WT always first
    which_wt <- which(d1[, 1] == "WT")
    if(length(which_wt)) {
        dwt <- d1[which_wt, ]
        dnwt <- d1[-which_wt, ]
    } else {
        dwt <- NULL
        dnwt <- d1
    }
    
    oo <- order(str_count(dnwt[, 1], ","), dnwt[, 1])
    dnwt <- dnwt[oo, ]
    d1 <- rbind(dwt, dnwt)
    
    m1 <- matrix(d1[, 2], ncol = 1)
    colnames(m1) <- "Freq."
    rownames(m1) <- d1[, 1]
    op <- par(mar = c(3, 5, 5, 3), las = 1)
    plot(m1, cex = 1.5, digits = 0, key = NULL,
         axis.col = list(side = 3), xlab = "", ylab = "",
         main  = "")
    par(op)
}


## output from all_methods_2_trans_mat, data for methods -> figures
##   OT: tree of restrictions, transition probs between genotypes
##        table of genotype freqs. (for lack of better idea)            
##   CBN: DAG of restrictions, transition probs between genotypes,
##                             transition probs genots TD
##   MCCBN: same as CBN 
##   MHN: theta matrix, transition probs between genotypes,
##                             transition probs genots TD

plot_DAG_fg <- function(x, data, orientation = "vertical", 
                        matrix=TRUE,
                        prune_edges = TRUE) {
    plot_fg <- function(fg) {
        ## Ideas from: https://stackoverflow.com/a/48368540
        lyt <- layout.reingold.tilford(fg)
        opx <- par(mar=c(2, 2, 2, 5))
        plot(fg,
             ## If I plot them sideways, labels in self-transitions
             ## overlap. FIXME. This sucks, I want them sideways!
             ## layout = -lyt[, 2:1],
             layout = lyt,
             edge.label = round(E(fg)$weight, 2),
             vertex.color = "SkyBlue2",
             edge.label.color = "black")
        par(opx)
    }

    ot_tree <- graph_from_data_frame(x$OT_model[, c(1, 2)])
    ot_fg <- graph_from_adjacency_matrix(x$OT_trans_mat, weighted = TRUE)

    cbn_tree <- graph_from_data_frame(x$CBN_model[, c(1, 2)])
    cbn_fg <- graph_from_adjacency_matrix(x$CBN_trans_mat, weighted = TRUE)
    cbn_td_fg <- graph_from_adjacency_matrix(x$CBN_td_trans_mat, weighted = TRUE)

    if((length(x$MCCBN_model) != 1) && !is.na(x$MCCBN_model)) {
        mccbn_tree <- graph_from_data_frame(x$MCCBN_model[, c(1, 2)])
        mccbn_fg <- graph_from_adjacency_matrix(x$MCCBN_trans_mat, weighted = TRUE)
        mccbn_td_fg <- graph_from_adjacency_matrix(x$MCCBN_td_trans_mat, weighted = TRUE)
    }
    
    ## mhn_mat <- graph_from_adjacency_matrix(out_schill$theta, weighted = TRUE,
    ##                                        mode = "directed")

    ## to rm edges with almost 0 weight
    
    mhn_tm <- as.matrix(x$MHN_trans_mat)
    mhn_td_tm <- as.matrix(x$MHN_td_trans_mat)
    
    if(prune_edges) {
        mhn_tm[mhn_tm < 0.01] <- 0
        mhn_td_tm[mhn_td_tm < 0.01] <- 0
    }
    mhn_fg <- graph_from_adjacency_matrix(mhn_tm,
                                          weighted = TRUE)
    mhn_td_fg <- graph_from_adjacency_matrix(mhn_td_tm,
                                             weighted = TRUE)
    
    if(orientation == "horizontal") {
        if((length(x$MCCBN_model) != 1) && !is.na(x$MCCBN_model))
            op1 <- par(mfcol = c(3, 4))
        else
            op1 <- par(mfcol = c(3, 3))
    } else {
        if((length(x$MCCBN_model) != 1) && !is.na(x$MCCBN_model))
            op1 <- par(mfrow = c(4, 3))
        else
            op1 <- par(mfrow = c(3, 3))
    }
    
    plot(ot_tree, layout = layout.reingold.tilford, vertex.size = 0, main = "OT")
    if (matrix){
        ot <- as.matrix(out$OT_trans_mat)
        plot(ot[rowSums(ot)>0, colSums(ot)>0], 
             digits=1, xlab="", ylab="",
             axis.col=list(side=1, las=2), axis.row = list(side=2, las=1), 
             main="OT trans matrix", cex.axis=0.7,
             mgp = c(2, 1, 0), key=NULL)
    }
    else {
        plot_fg(ot_fg)  
    }
    
    par(mar=rep(3, 4))
    plot_sampled_genots(data)
    ## plot(type ="n", c(0, 0), c(0, 0), axes = FALSE)

    plot(cbn_tree, layout = layout.reingold.tilford, vertex.size = 0, main = "CBN")
    
    if (matrix){
        cbn_trans_mat <- as.matrix(out$CBN_trans_mat)
        plot(cbn_trans_mat[rowSums(cbn_trans_mat)>0,colSums(cbn_trans_mat)>0], 
             digits=1, cex.axis=0.7,
             main="CBN: trans matrix",
             xlab="", ylab="",
             axis.col=list(side=1, las=2), axis.row = list(side=2, las=1), 
             mgp = c(2, 1, 0), key=NULL)
        cbn_td_trans_mat <- as.matrix(out$CBN_td_trans_mat)
        plot(cbn_td_trans_mat[rowSums(cbn_td_trans_mat)>0,colSums(cbn_td_trans_mat)>0], 
             digits=1, cex.axis=0.7,
             main="CBN: trans td matrix", 
             xlab="", ylab="",
             axis.col=list(side=1, las=2), axis.row = list(side=2, las=1), 
             mgp = c(2, 1, 0), key=NULL)
        }
    else {
        plot_fg(cbn_fg)  
        plot_fg(cbn_td_fg)
    }

    if((length(x$MCCBN_model) != 1) && !is.na(x$MCCBN_model)) {
        plot(mccbn_tree, layout = layout.reingold.tilford, vertex.size = 0, main = "MCCBN")
        plot_fg(mccbn_fg)
        plot_fg(mccbn_td_fg)
    }
    
    ## I am unable to use bending edges for the theta matrix.
    ## plot_fg(mhn_mat)
    op <- par(mar=c(4, 4, 5, 3), las = 1)
    plot(x$MHN_theta, cex = 1.5, digits = 2, key = NULL,
         axis.col = list(side = 3),
         xlab = "Effect of this (effector)",
         ylab = " on this (affected)",
         mgp = c(2, 1, 0))
    par(op)
    if (matrix){
        plot(mhn_tm[rowSums(mhn_tm)>0,colSums(mhn_tm)>0], 
             digits=1, 
             main="MHN: trans matrix",cex.axis=0.7, 
             axis.col=list(side=1, las=2), axis.row = list(side=2, las=1), 
             xlab="", ylab="",
             mgp = c(2, 1, 0), key=NULL)
        plot(mhn_td_tm[rowSums(mhn_td_tm)>0,colSums(mhn_td_tm)>0], 
             digits=1,  cex.axis=0.7,
             axis.col=list(side=1, las=2), axis.row = list(side=2, las=1), 
             main="MHN: trans td matrix", 
             xlab="", ylab="",
             mgp = c(2, 1, 0), key=NULL)
    }
    else {
        plot_fg(mhn_fg)
        plot_fg(mhn_td_fg)
    }
    par(op1)
}


## Remove a fraction, frac, of the WT
## from a matrix of data. Used to examine
## the effect of having very few WT
remove_WT <- function(x, frac = 0.9) {
    which_wt <- which(rowSums(x) == 0)
    if(length(which_wt) > 0) {
        rmwt <- which_wt[1:(round(frac * length(which_wt)))]
    }
    return(x[-rmwt, ])
}

## Add N WT to the data. Used to examine the effect of having many WT.
add_WT <- function(x, N = 10000) {
    ncx <- ncol(x)
    x <- rbind(x, matrix(0, nrow = N, ncol = ncx))
    return(x)
}

library(codetools)
checkUsageEnv(env = .GlobalEnv)



#####################################################################

## ## mini example
## Dat1 <- readRDS(file="./MHN/data/BreastCancer.rds") [1:50, 1:4]
## colnames(Dat1) <- LETTERS[1:ncol(Dat1)]

## Dat1 <- as.matrix(Dat1)

## simplerun <- all_methods_2_trans_mat(Dat1)

## rm(simplerun)





## Getting an idea of sizes of data if we did not use sparse matrices
if(FALSE) {
## 8 bytes per float as can be checked doing
u <- runif((2^10) * (2^10))
print(object.size(u), units = "b")
8 * length(u)

uu <-  runif((2^15) * (2^15))
length(uu)
print(object.size(uu), units = "b")
print(object.size(uu), units = "GB") ## 8 GB

length(uu) * 8 /(1024 * 1024 * 1024)

## so an object for 2^20 * 2^20 would take

((2^20)^2) * 8 / (1024 * 1024 * 1024)

## or 8192 GB which is what R complaints about if I try this
uuu <-  runif((2^20) * (2^20))

## Yes, there are no memory limits in the machines
## install.packages("devtools", dependencies = TRUE)
##devtools::install_github("krlmlr/ulimit")
library(ulimit)
ulimit::memory_limit()

}



## MHN does not run if we use mclapply after setting the threads for OMP >
## 1

## all_data_out <- mclapply(all_data, all_methods_2_trans_mat, mc.cores
## = 36)

## Yes, you can use mclapply, just make sure to set the threads of OMP to 1
## The following shows it. It just did not make any sense above.
if(FALSE) {

    lapply(all_data, dim)

    RhpcBLASctl::omp_set_num_threads(36)

    system.time(cucu3 <- lapply(all_data[c(3, 4, 20, 21:24)], do_MHN))

    RhpcBLASctl::omp_set_num_threads(1)
    system.time(cucu4 <- mclapply(all_data[c(3, 4, 20, 21:24)], do_MHN,
                                  mc.cores = 36))
    stopifnot(identical(cucu3, cucu4))
}


## If you wanted to launch a single process, you could do this
## ## Since I have a lot fewer than 36 data sets, and easier to catch issues
## ## as they show up, do not parallelize. But use multiple threads
## all_data_out <- lapply(1:length(all_data),
##                        function(i) {
##                            cat("\n #################################")
##                            cat("\n #################################")
##                            cat("\n\n   Doing data set ", names(all_data)[i])
##                            cat(" with dim = ", dim(all_data[[i]]))
##                            return(all_methods_2_trans_mat(all_data[[i]]))
##                        }
##                        )


## Sort wrapper to generate a transition matrix.
## It first generates all accesibles genotypes given a transition matrix
cpm_access_genots_paths_w_simplified_OR <- function(data, parameter_column_name){
  tmp <- try(df_2_access_genots_and_graph_OR(data$edges[, c("From", "To")]))
  
  if(inherits(tmp, "try-error")) {
    stop("how is this happening? there was edges component!")
  } else {
    accessible_genots <- tmp$accessible_genots
    fgraph <- unrestricted_fitness_graph_sparseM(accessible_genots)
  }
  
  if (is.null(data$edges[,parameter_column_name])){
    stop("No such column")
  }
  
  weights <- unique(data$edges[, c("To", parameter_column_name)])
  rownames(weights) <- weights[, "To"]
  weighted_fgraph <- transition_fg_sparseM(fgraph, weights)
  trans_mat_genots <- rowScaleMatrix(weighted_fgraph)

  data$fgraph <- fgraph
  data$weighted_fgraph <- weighted_fgraph
  data$trans_mat_genots <- trans_mat_genots
  return(data)
}