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
MCCBN_INSTALLED <- FALSE


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
source("hypertraps-process.R", echo = FALSE, max.deparse.length = 0)
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
    
    # if(caprese_capri_minimal) {
    #     ## Most of the time we only want the graphs, that's it
    #     evaleloss <- FALSE
    #     evalprederr <- FALSE
    #     evalposterr <- FALSE
    # } else {
    #     evaleloss <- TRUE
    #     evalprederr <- TRUE
    #     evalposterr <- TRUE
    # }
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

## DAG of restrictions (as data frame) -> vector of accessible genotypes and graph of DAG of restrictions
##  Under an OR model (not XOR, not AND), such as DBN
## return all the accessible genotypes
##     from a DAG of genes plus the DAG as igraph object
df_2_access_genots_and_graph_relationships <- function(x, gene_relations = NULL) {
    ## minor detail: if x is a sparse adjacency mat.
    ##    from igraph this still works. But don't do that.

    genotype_follows_relationship <- function(genotype){
        tmp_gene_relations <- gene_relations[genotype]
        
        for(relationship in c("XOR", "AND")){
            rel_genes <- names(which(tmp_gene_relations == relationship) == TRUE)
            for ( rel_idx in rel_genes){
                parents_of_rel <- adjacent_vertices(g, rel_idx, mode=c("in"))[[rel_idx]]$name
                if(relationship == "XOR" && all(parents_of_rel %in% genotype)) return(FALSE) ## XOR relationship not futfilled
                if(relationship == "AND" && !all(parents_of_rel %in% genotype)) return(FALSE) ## AND relationship not futfilled
            }
        }
        return(TRUE)
    }

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

            ## Evaluate relationships
            if(genotype_follows_relationship(tmp_gty_1)){
                all_gty[[i]] <- tmp_gty_1 
                i <- i + 1
            } 

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
                tmp_gty_2 <- tmp_gty_2[unlist(lapply(tmp_gty_2, genotype_follows_relationship))]
                if(length(all_gty) < (i + length(tmp_gty_2))) {
                    all_gty <- c(all_gty,
                                 vector(mode = "list", length = max(length(all_gty),
                                                                    2 + length(tmp_gty_2))))
                }
                if(length(tmp_gty_2)){
                    all_gty[(i):(i + length(tmp_gty_2) - 1)] <- tmp_gty_2
                    i <- i + length(tmp_gty_2)
                }
            }
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


## ### Examples
## x2 <- data.frame(From = c("Root", "Root", "A", "B", "C"),
##                  To = c("A", "B", "C", "C", "D"))

## x3 <- data.frame(From = c("Root", "Root", "A", "D", "B", "E"),
##                  To = c("A", "D", "B", "E", "C", "C"))


## df_2_access_genots_and_graph_OR(x2)
## df_2_access_genots_and_graph_OR(x3)




## From function of same name in ruggify-functions.R

## list of accessible genotypes -> adjacency matrix of genotypes (fitness graph)
## return maximally connected fitness graph for a given set of accessible genotypes

## BEWARE! This is the maximally connected fitness graph. This works with
## CPMs. But this is wrong, if, say, fitnesses are: A = 2, B = 3, AB = 2.5.
## This will place an arrow between B and AB, but there should  be no such edge.
## See below for how to have a general procedure

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




## A simple, general procedure, to obtain the fitness graph is:
##   - create matrix of all genotypes
##   - Fill up, upper diagonal, with the fitness of destination genotype
##   - Subtract fitness of original genotype
##   - Set as 0 those with value <= 0.
## We could use the list given below as starting point (to use a smaller matrix)
## The list of accessible genotypes can also be obtained from wrap_accessibleGenotypes
## in OncoSimulR (unexported function)
## This function is now available as genots_2_fgraph_and_trans_mat



## Trying sparse matrices
## https://stackoverflow.com/questions/23107837/r-sparse-matrix-from-list-of-dimension-names
## https://stackoverflow.com/questions/26207850/create-sparse-matrix-from-a-data-frame
## https://cmdlinetips.com/2019/05/introduction-to-sparse-matrices-in-r/
## https://www.gormanalysis.com/blog/sparse-matrix-construction-and-use-in-r/


## BEWARE! This is the maximally connected fitness graph. This works with
## CPMs. But this is wrong, if, say, fitnesses are: A = 2, B = 3, AB = 2.5.
## This will place an arrow between B and AB, but there should  be no such edge.
## This function is now available as genots_2_fgraph_and_trans_mat


## list of accessible genotypes -> adjacency matrix of genotypes (fitness graph)
## return maximally connected fitness graph for a given set of accessible genotypes
##   This list contains no WT (we add it)
unrestricted_fitness_graph_sparseM <- function(gacc, plot = FALSE) {
    gs <- unlist(lapply(gacc, function(g) paste0(g, collapse = ", ")))
    gs <- c("WT", gs)
    nmut <- c(0, vapply(gacc, length, 1))
    
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
##   BEWARE: this assumes a single global maximum!!!
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

## Sort wrapper to generate a transition matrix.
## It first generates all accesibles genotypes given a transition matrix
cpm_access_genots_paths_w_simplified_OR <- function(data,               
    parameter_column_name = c("Thetas", "Probabilities")){

  tmp <- try(df_2_access_genots_and_graph_OR(data$edges[, c("From", "To")]))
  
  if(inherits(tmp, "try-error")) {
    stop("how is this happening? there was edges component!")
  } else {
    accessible_genots <- tmp$accessible_genots
    fgraph <- unrestricted_fitness_graph_sparseM(accessible_genots)
  }
  
  which_col_weights <- which(colnames(data$edges) %in% parameter_column_name)
  if (is.null(data$edges[,which_col_weights])){
    stop("No such column")
  }
  weights <- unique(data$edges[, c(2, which_col_weights)])
  rownames(weights) <- weights[, "To"]
  weighted_fgraph <- transition_fg_sparseM(fgraph, weights)
  trans_mat_genots <- rowScaleMatrix(weighted_fgraph)

  return(list(
        fgraph = fgraph,
        weighted_fgraph = weighted_fgraph,
        trans_mat_genots = trans_mat_genots
  ))
}

## Sort wrapper to generate a transition matrix.
## It first generates all accesibles genotypes given a transition matrix
cpm_access_genots_paths_w_simplified_relationships <- function(data,               
    parameter_column_name = c("Thetas", "Probabilities", "Lambdas")){

  tmp <- try(df_2_access_genots_and_graph_relationships(data$edges[, c("From", "To")]
                                                        , data$parent_set))
  
  if(inherits(tmp, "try-error")) {
    stop("how is this happening? there was edges component!")
  } else {
    accessible_genots <- tmp$accessible_genots
    fgraph <- unrestricted_fitness_graph_sparseM(accessible_genots)
  }
  
  which_col_weights <- which(colnames(data$edges) %in% parameter_column_name)
  if (is.null(data$edges[,which_col_weights])){
    stop("No such column")
  }
  weights <- unique(data$edges[, c(2, which_col_weights)])
  rownames(weights) <- weights[, "To"]
  weighted_fgraph <- transition_fg_sparseM(fgraph, weights)
  trans_mat_genots <- rowScaleMatrix(weighted_fgraph)

  return(list(
        fgraph = fgraph,
        weighted_fgraph = weighted_fgraph,
        trans_mat_genots = trans_mat_genots
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

all_methods_2_trans_mat <- function(x, cores_cbn = 1, do_MCCBN = FALSE, HT_folder = NULL) {
     
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
    
    cat("\n     Doing DBN\n\n")
    time_dbn <- system.time(
      out_dbn <- do_DBN(x))["elapsed"]
    
    cat("\n  time DBN = ", time_dbn)

    cat("\n     Doing HyperTraps")
    print("By default we run it here with dry_run = TRUE.
        HyperTraPS takes a long time and I do no want to block the script.
    ")
    # time_HyperTraPS <- system.time(
    #   out_HyperTraPS <- do_HyperTraPS(x, tmp_folder = HT_folder, dry_run = TRUE, plot = FALSE))["elapsed"]

    # cat("\n  time HyperTraPS = ", time_HyperTraPS)

    # cat("\n     Doing HESBCN")
    # time_hesbcn <- system.time(
    #   out_hesbcn <- do_HESBCN(x))["elapsed"]
      out_hesbcn <- NULL
    
    # cat("\n  time HESBCN = ", time_hesbcn)
     
    cpm_out_others <- all_methods(x, cores_cbn = cores_cbn, do_MCCBN = do_MCCBN)
    pre_trans_mat_others <- lapply(cpm_out_others[methods],
        cpm_access_genots_paths_w_simplified)

    pre_trans_mat_new_CPMS <- lapply( #For the moment just for DBN (HESBCN in the future)
        list(DBN = out_dbn
        # , HyperTraPS = out_HyperTraPS
        ),
        cpm_access_genots_paths_w_simplified_OR)
    
    # pre_trans_mat_HESBCN <- lapply( #For the moment just for DBN (HESBCN in the future)
    #     list(HESBCN = out_hesbcn
    #     # , HyperTraPS = out_HyperTraPS
    #     ),
    #     cpm_access_genots_paths_w_simplified_relationships)

    pre_trans_mat_others["DBN"] <- list(pre_trans_mat_new_CPMS$DBN)
    # pre_trans_mat_others["HESBCN"] <- list(pre_trans_mat_HESBCN$HESBCN)
    # pre_trans_mat_others["HyperTraPS"] <- list(pre_trans_mat_new_CPMS$HyperTraPS)
    cat("\n    getting transition matrices for all non-mhn methods \n")

    ## ## Unweighted
    ## uw <- lapply(pre_trans_mat_others, function(x) rowScaleMatrix(x$fgraph))
    ## Weighted
    wg <- lapply(pre_trans_mat_others[c("OT", "MCCBN" , "CBN_ot", "DBN", "HyperTraPS", "HESBCN")[c(TRUE, do_MCCBN, TRUE, TRUE, FALSE, FALSE)]], 
        function(x) x$trans_mat_genots)
    ## Diagonal
    td <- lapply(pre_trans_mat_others[c("MCCBN", "CBN_ot", "DBN", "HyperTraPS", "HESBCN")[c(do_MCCBN, TRUE, TRUE, FALSE, FALSE)]],
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
        DBN_f_graph = pre_trans_mat_new_CPMS$DBN$weighted_fgraph,
        DBN_trans_mat = pre_trans_mat_new_CPMS$DBN$trans_mat_genots,
        DBN_td_trans_mat = td$DBN
        # HESBCN_model = out_hesbcn$edges,
        # HESBCN_parent_set = out_hesbcn$parent_set,
        # HESBCN_f_graph = pre_trans_mat_HESBCN$HESBCN$weighted_fgraph,
        # HESBCN_trans_mat = pre_trans_mat_HESBCN$HESBCN$trans_mat_genots,
        # HESBCN_td_trans_mat = td$HESBCN
        # HyperTraPS_model = out_HyperTraPS$edges,
        # HyperTraPS_f_graph = pre_trans_mat_new_CPMS$HyperTraPS$weighted_fgraph,
        # HyperTraPS_trans_mat = pre_trans_mat_new_CPMS$HyperTraPS$trans_mat_genots,
        # HyperTraPS_td_trans_mat = td$HyperTraPS,
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
cpm_layout <- function(graph){
    # V(graph)$num_mutations <<- 
    num_mutations <- sapply(V(graph)$name, function(x){
        ifelse(x == "WT", 0, nchar(x))
    })
    V(graph)$num_mutations <- num_mutations 
    
    lyt <- matrix(0, ncol = 2, nrow = length(V(graph)))
    lyt[, 2] <-  num_mutations

    for (i in 0:max(num_mutations)) {
        level_idx <- which(num_mutations == i)
        gnt_names <- sort(V(graph)$name[level_idx], index.return = TRUE)
        spacing <- 6 / (length(level_idx) + 1)
        level_elements <- 1:length(level_idx) * spacing
        level_elements <-  rev(level_elements - max(level_elements))
        correction_rate <- max(level_elements) - min(level_elements) 
        level_elements <- level_elements + correction_rate/2
        lyt[, 1][level_idx[gnt_names$ix]] <- level_elements
        # browser()
    }
    return(lyt)
}

#' Plot transitions between genotypes
#' 
#' @param trans_mat transitions matrix to plot: can contain row counts, probabilities... 
#' Entries will be normalized
#' @param observations Original cross sectional data used to compute the model. Optional.
#' @param freqs DataFrame with column $Genotype and $Freqs with their frequencies. Optional.
#' @examples
#' png("fluxes.png")
#' par(mfrow = c(1, 3))
#' plot_genot_fg(out$MHN_trans_mat, db2, sorted_observations)
#' title("Transition Matrix", line = -3)
#' plot_genot_fg(edge_transitions, db2, sorted_observations)
#' title("Fluxes", line = -3)
#' dev.off()
plot_genot_fg <- function(trans_mat
    , observations = NULL
    , freqs = NULL
    , simplify = TRUE
    , top_edge_per_node = 5){
    # trans_mat <- as.matrix(trans_mat)
    rownames(trans_mat) <- str_replace_all(rownames(trans_mat), ", ", "")
    colnames(trans_mat) <- str_replace_all(colnames(trans_mat), ", ", "")

    graph <- graph_from_adjacency_matrix(trans_mat, weighted = TRUE)
    # print(sprintf("->Number of nodes %s Number of edges %s", length(V(graph)), length(E(graph))))

    unique_genes_names <- sort(unique(str_split(paste(rownames(trans_mat)[-1], collapse=""), "")[[1]]))

    num_genes <- length(unique_genes_names)
    graph <- graph_from_adjacency_matrix(trans_mat, weighted = TRUE)
    # browser()
    ## TODO take into account that a genoytpe can be observed
    if(simplify){
        graph <- decompose(graph)[[1]]
        min_values <- sort(trans_mat[trans_mat > 0], decreasing = TRUE)
        thr <- num_genes * top_edge_per_node
        ifelse(length(min_values > thr)
            , min_value <- min_values[thr]
            , min_value <- min_values[-1]
        )
        trans_mat[trans_mat < min_value] = 0

        graph <- graph_from_adjacency_matrix(trans_mat, weighted = TRUE)
        
        subgraphs <- decompose(graph)
        graph <- subgraphs[[1]]
        for(i in subgraphs[-1]){
            if(length(V(i)) > 5){
                EL  = get.edgelist(graph)
                EL1 = get.edgelist(i)
                ELU = rbind(EL, EL1)
                ELU = ELU[!duplicated(ELU),]
                w <- c(get.edge.attribute(graph, "weight")
                    , get.edge.attribute(i, "weight"))
                graph <- graph_from_edgelist(ELU)
                graph <- set_edge_attr(graph, "weight", value = w)
            } 
        }
    } 
        
    if (!is.null(observations)){
        observations <- as.data.frame(sampledGenotypes(observations))
        observations$Abs_Freq <- observations$Freq / sum(observations$Freq)
        observations$Genotype <- str_replace_all(observations$Genotype, ", ", "")
    }

    if (!is.null(freqs)){
        freqs$Abs_Freq <- freqs$Counts / sum(freqs$Counts)
        freqs$Genotype <- str_replace_all(freqs$Genotype, ", ", "")
    }

    lyt <- cpm_layout(graph)
    # browser()

    observed_color <- "#ff7b00"
    not_observed_color <- "#0892d0" 
    if(is.null(observations)){
        not_observed_color <- "#ff7b00"
    } 
    colors <- sapply(V(graph)$name, 
        function(gen){
            if (sum(match(observations$Genotype, gen, nomatch = 0)) == 1){
                return(observed_color)
            } 
            return(not_observed_color)
        })
    V(graph)$color <- colors
    V(graph)$frame.color <- colors

    sizes <- vapply(V(graph)$name, 
        function(gen){
            if (sum(match(freqs$Genotype, gen, nomatch = 0)) == 1){
                return(freqs$Abs_Freq[which(freqs$Genotype == gen)] * 150)
            } else if ((sum(match(observations$Genotype, gen, nomatch = 0)) == 1)){
                return(observations$Abs_Freq[which(observations$Genotype == gen)] * 150)
            } else {
                return(7.5)
            }
        }, numeric(1.0))

    if(all(sizes == 7.5)) sizes <- rep(15, length(sizes))
    V(graph)$size <- sizes

    V(graph)$label.family <- "Helvetica"

    opx <- par(mar = c(1, 1, 1, 1))

    w <- E(graph)$weight
    w <- w / max(w) * 10
    if(all(w == 10)) w <- rep(1, length(w))
    plot(graph
        , layout = lyt[, 2:1]
        , vertex.label.color = "black"
        , vertex.label.family = "Helvetica"
        , font.best = 2
        , vertex.label.cex = 1
        , vertex.frame.width = 0
        , edge.color = rgb(0.5, 0.5, 0.5, 1)
        # , edge.color = rgb(0.5, 0.5, 0.5, E(graph)$weight/max(E(graph)$weight))
        , edge.arrow.size = 0.5
        , xlab = "Number of features acquired"
        , edge.width = w
    )

    margin <- -1.15
    lines(c(-1.2, 1.2), c(margin, margin), lwd = 2)
    node_depth <- sapply(V(graph)$name
        , function(x) distances <- distances(graph, algorithm = "unweighted", to = x)["WT",]
        )
    max_node_depth <- max(node_depth[is.finite(node_depth)])
    axis(1
        , at = seq(-1, 1, length.out = max_node_depth + 1)
        , labels = 0:(max_node_depth)
        , lwd = 2
        , cex = 2
        , pos = margin)
    if(!(is.null(observations))){
        legend("bottom", c("Observed", "Not observed") 
            , box.lwd = 0, lty=c(NA, NA), lwd = c(NA, NA)
            , pch = c(21, 21)
            , col = c(observed_color, not_observed_color)
            , pt.bg = c(observed_color, not_observed_color)
            , pt.cex = c(2, 2), horiz = TRUE
            , x.intersp = c(0, 0)
            )
    }
    par(opx)
    # title(xlab = "Number of features acquired", line = -3)
}

#' Plot results from CPMs
#' 
#' By default it create a top row with the DAG of the CPM 
#' or de transtionRateMatrix for MHN
#' The bottom row has a custom plot for transition between genotypes

#' @param x output from the cpm
#' @param data Original cross sectional data used to compute the model. Optional.
#' @param models Output of the CPMs to plot. Current support is for OT, CBN, DBN, MCCBN and MHN Optional.
#' @param orientation String. If it not "vertical" will be displayed with an horizontal layout. Optional.
#' @param plot_type String. You can choose between 4 options. Optional.
#' genotypes: DAG with genotypes transitions
#' matrix: respresents the transtion matrix as a heatmap
#' transitions: shows the transitions count between genotypes in HyperTraPS style
#'              running simulations is needed before this
#' trans_mat: HyperTraps-like representation of the transition matrix
#' 
#' @examples
#' out <- all_methods_2_trans_mat(dB_c1, do_MCCBN = TRUE)
#' png("trans_at.png", width = 1000, height = 600, units = "px")
#' plot_DAG_fg(out, dB_c1, plot_type = "trans_mat")
#' dev.off()
#' out2 <- run_all_simulations(out, 100, n_genes = 5)
#' png("graph.png", width = 1000, height = 600, units = "px")
#' plot_DAG_fg(out, dB_c1, plot_type = "transitions")
#' dev.off()
plot_DAG_fg <- function(x, data, orientation = "horizontal", 
                        models = c("OT", "CBN", "DBN", "MCCBN", "MHN", "HESBCN"),
                        plot_type = "trans_mat",
                        prune_edges = TRUE) {
    
    if (!(plot_type %in% c("matrix", "transitions", "trans_mat", "genotypes"))){
        stop(sprintf("Plot type %s is not supported", plot_type))
    }
    # plot_fg <- function(fg) {
    #     ## Ideas from: https://stackoverflow.com/a/48368540
    #     lyt <- layout.reingold.tilford(fg)
    #     opx <- par(mar=c(2, 0.5, 2, 0.5))
    #     plot(fg,
    #          ## If I plot them sideways, labels in self-transitions
    #          ## overlap. FIXME. This sucks, I want them sideways!
    #          ## layout = -lyt[, 2:1],
    #          layout = lyt,
    #          edge.label = round(E(fg)$weight, 2),
    #          vertex.color = "SkyBlue2",
    #          edge.label.color = "black")
    #     par(opx)
    # }

    process_data <- function(mod) {
        dag_tree <- NULL
        
        dag_tree <- NULL
        tryCatch (expr = {
            dag_model <- get(paste(mod, "_model", sep = ""), x)
            dag_tree <- graph_from_data_frame(dag_model[, c(1, 2)])
        }, error = function(e){})

        dag_trans_mat <- get(paste(mod, "_trans_mat", sep = ""), x)
        fg <- graph_from_adjacency_matrix(dag_trans_mat, weighted = TRUE)

        if(prune_edges) {
            dag_trans_mat[dag_trans_mat < 0.01] <- 0
        }

        if(plot_type == "matrix") {
            dag_trans_mat <- as.matrix(dag_trans_mat)
            dag_trans_mat <- dag_trans_mat[rowSums(dag_trans_mat) > 0, colSums(dag_trans_mat) > 0]
        }

        td_trans_mat <- NULL
        td_fg <- NULL
        tryCatch(expr = {
            td_trans_mat <- get(paste(mod, "_td_trans_mat", sep = ""), x)
            td_fg <- graph_from_adjacency_matrix(td_trans_mat, weighted = TRUE)
            if(prune_edges) {
                td_trans_mat[td_trans_mat < 0.01] <- 0
            }
            if (plot_type == "matrix"){
                td_trans_mat <- as.matrix(td_trans_mat)
                td_trans_mat <- td_trans_mat[rowSums(td_trans_mat) > 0, colSums(td_trans_mat) > 0]
            }
        }, error = function(e){ })

        theta <- NULL
        tryCatch(expr ={
            theta <- get(paste(mod, "_theta", sep=""), x)
        }, error = function(e) { })

        return(list(dag_tree = dag_tree
            , dag_trans_mat = dag_trans_mat
            , fg = fg
            , td_trans_mat = td_trans_mat
            , td_fg = td_fg
            , theta = theta
            , parent_set = x[[sprintf("%s_parent_set", mod)]]
            , transitions = x[[sprintf("%s_genotype_transitions", mod)]]
            ))
    }

    ## List of available models
    available_models <- models[
        vapply(models, function(mod) {
            attr_name <- ifelse(mod == "MHN", "theta", "model")
            return(any(!is.na(x[[sprintf("%s_%s", mod, attr_name)]])))
        }, logical(1))
    ]
    print(available_models)

    ## DAG relationships colors 
    standard_relationship <- "gray73"
    colors_relationships <- c(standard_relationship, "coral2", "cornflowerblue", "darkolivegreen3")
    names(colors_relationships) <- c("Single", "AND", "OR", "XOR")
    
    ## Shape of the plot
    l_models <- length(available_models)
    n_rows <- ifelse(plot_type == "matrix", 3, 2)

    if(orientation == "vertical") {
        op1 <- par(mfrow = c(l_models, n_rows))
    } else {
        op1 <- par(mfcol = c(n_rows, l_models))
    }
    par(mar = c(0.5, 1, 0.5, 0.5), mai = c(0.25, 0.25, 0.25, 0.25))

    ## Plotting models
    for(mod in available_models) {
        ## Processing data
        model_data2plot <- process_data(mod)
        ## Plotting data
        if(!is.null(model_data2plot$dag_tree)) {
            if(!is.null(model_data2plot$parent_set)){
                for(i in names(model_data2plot$parent_set)){
                    E(model_data2plot$dag_tree)[.to(i)]$color <- colors_relationships[model_data2plot$parent_set[[i]]]
                }
            } else E(model_data2plot$dag_tree)$color <- standard_relationship
            plot(model_data2plot$dag_tree
                , layout = layout.reingold.tilford
                , vertex.size = 30 
                , vertex.label.color = "black"
                , vertex.label.family = "Helvetica"
                , font.best = 2
                , vertex.frame.width = 0.5
                , vertex.color = "white"
                , vertex.frame.color = "black" 
                , vertex.label.cex = 1
                , edge.arrow.size = 0.45
                , edge.width = 5
                , main = mod)
            if(!is.null(model_data2plot$parent_set)){
                legend("topleft", legend = names(colors_relationships),
                    col = colors_relationships, lty = 1, lwd = 2)
            }
        }else if(!is.null(model_data2plot$theta)) {
            op <- par(mar=c(3, 3, 5, 3), las = 1)
            plot(model_data2plot$theta, cex = 1.5, digits = 2, key = NULL
                , axis.col = list(side = 3)
                , xlab = "Effect of this (effector)"
                , ylab = " on this (affected)"
                , main = mod
                , mgp = c(2, 1, 0))
            par(op)
        }

        if (plot_type == "matrix"){
            plot(model_data2plot$dag_trans_mat
                , digits = 1, xlab = "", ylab = ""
                , axis.col = list(side = 1, las = 2)
                , axis.row = list(side = 2, las = 1) 
                , main = paste(mod, ": trans matrix", sep = " ")
                , cex.axis = 0.7
                , mgp = c(2, 1, 0)
                , key = NULL)
            
            if (!is.null(model_data2plot$td_trans_mat)){
                plot(model_data2plot$td_trans_mat, 
                    digits = 1, cex.axis = 0.7,
                    main = paste(mod, ": trans td matrix", sep = " "), 
                    xlab = "", ylab = "",
                    axis.col = list(side = 1, las = 2), axis.row = list(side = 2, las = 1), 
                    mgp = c(2, 1, 0), key = NULL)
            }
        } else if (plot_type == "genotypes") {
            plot_genot_fg(as_adjacency_matrix(model_data2plot$fg), simplify = FALSE)
        } else if (plot_type == "transitions") {
            plot_genot_fg(model_data2plot$transitions, data)
        }else if (plot_type == "trans_mat"){
            plot_genot_fg(model_data2plot$dag_trans_mat, data, simplify = FALSE)
        }

        if ((mod %in% c("OT")) & (plot_type == "matrix")) {
            par(mar = rep(3, 4))
            plot_sampled_genots(data)
        }
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
