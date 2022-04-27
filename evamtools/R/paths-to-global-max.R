## Copyright 2016, 2017, 2018, 2020, 2022 Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify it under
## the terms of the GNU Affero General Public License (AGPLv3.0) as published by
## the Free Software Foundation, either version 3 of the License, or (at your
## option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero General Public License for more details.

## You should have received a copy of the GNU Affero General Public License along
## with this program.  If not, see <http://www.gnu.org/licenses/>.


## The code that follows is tested, and is the basis of much simplified
## code in evam-wrapper.R
## But this code itself is not actually used here.
## This was used in Diaz-Uriarte and Vasallo, 2019
## This was used to obtain paths to the global maximum from the CPMs
## all of which assumed a single global maximum.


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

    if(inherits(x, "try-error") || all(is.na(x)) || is.null(x)) {
        ## The CPM analysis produced no edges component, so
        ## nothing can be done
        if(inherits(x, "try-error")) likely_error <- "Error_in_run"
        if(all(is.na(x))) likely_error <- "ncol_x"
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
    tmp <- try(DAG_2_access_genots(x[, c("From", "To")]))
    
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
        ## weighted_fgraph need not have each row sum to 1.
        ## Obvious if they are lambdas
        ## from CBN for instance.
        ## Neither do they sum to 1 from some OT models.
        ## So make sure they are transition matrices
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
            evam_shannonI(weighted_paths_to_max[, "probability"])

    } else {
        weighted_paths_to_max <- NA
        diversity_weighted_paths_to_max <- NA
    }
    
    return(list(accessible_genots = accessible_genots,
                num_accessible_genots = length(accessible_genots),
                CPM_DAG_as_igraph = tmp$graph,
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



## From function of same name in ruggify-functions.R

## list of accessible genotypes -> adjacency matrix of genotypes (fitness graph)
## return maximally connected fitness graph for a given set of accessible genotypes

## BEWARE! This is the maximally connected fitness graph. This works with
## CPMs. But this is wrong, if, say, fitnesses are: A = 2, B = 3, AB = 2.5.
## This will place an arrow between B and AB, but there should  be no such edge.
## See unrestricted_fitness_graph_sparseM for
## similar comments and a function that only returns
## the truly accessible

unrestricted_fitness_graph <- function(gacc) {
    
    gs <- unlist(lapply(gacc, function(g) paste0(g, collapse = ", ")))
    gs <- c("WT", gs)
    nmut <- c(0, vapply(gacc, length, 1))
    gs <- gs[order(nmut, gs)]
    
    adjmat <- matrix(0L, nrow = length(gs), ncol = length(gs))
    rownames(adjmat) <- colnames(adjmat) <- gs
    
    adjmat["WT", gs[which(nmut == 1)]] <- 1L
    
    for(m in 2:max(nmut)){
        g <- gs[which(nmut == m)]
        for (gn in g) {
            parents <- gacc[which(nmut == m-1)-1]
            gns <- unlist(strsplit(gn, ", "))
            parents <-
                parents[which(unlist(lapply(parents,
                                            function(p)
                                                length(setdiff(gns, p)))) == 1)]
            for (p in parents){
                adjmat[paste0(p, collapse = ", "), gn] <- 1L
            }
        }
    }
    ## if (plot)
    ##     mccbn::plot_poset(adjmat) ## , title = "G0 (unrestricted)")

    stopifnot(all(adjmat %in% c(0L, 1L) ))
    storage.mode(adjmat) <- "integer"

    return(adjmat)
}


## fitness graph, weights of probs/lambdas descendant gene given parent gene
##                                            -> weighted fitness graph
##             the fitness graph with weights (relative to each node)
##             of jumping to each of the descendant genotypes
transition_fg <- function(x, weights) { 
    pos_do <- which(x == 1, arr.ind = TRUE)
    wfg <- x
    wfg[] <- 0
    tmp <- unlist(
        Map(function(r, c)
             get_single_lambda(x, r, c, weights),
            pos_do[, 1], pos_do[, 2]))
    wfg[pos_do] <- tmp
        
    return(wfg)
}


# library(codetools)
## checkUsageEnv(env = .GlobalEnv)
