## Copyright 2016, 2017, 2018, 2020, 2022 Ramon Diaz-Uriarte

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


## wrap calling the CPMs.
## Given a data set (patients as rows, genes as columns) return the
## transition matrices between genotypes according to all methods.

## Done by function: evam (formerly all_methods_2_trans_mat)
## almost at the bottom.


######################################################################
######################################################################
###
###
###       Cores, OMP, etc
###
###
######################################################################
######################################################################

## By default, this will use the available cores
## for MHN. For CBN we fix the OMP cores to one (setting
## that as an env. var. right before calling cbn via system)
##
## See more comments and possible issues at bottom and
## in the cbn-process.R file.
## For more control, you might want to set
## thiscores <- 36 ## or whatever
## Then
## export OPENBLAS_NUM_THREADS=thiscores
## export OMP_NUM_THREADS=thiscores
## in the script that launches analyses, unless those runs are parallelized.

## If runs are parallelized, and you do many parallel runs, it might
## be simpler and faster to set cores to 1 for everything (MHN included).
## (See bottom: MHN will not run with cores > 1 if also using forking
##  from mclapply)
## And you can also do it here for R itself
## library(RhpcBLASctl)
## RhpcBLASctl::blas_set_num_threads(thiscores)
## RhpcBLASctl::omp_set_num_threads(thiscores)


######################################################################
######################################################################
###
###
###       DAG_2_access_genots* and cpm2tm
###
###
######################################################################
######################################################################

## cpm2tm : CPM -> transition matrices
## DAG_2_access_genots*: CPM output -> accessible genotypes

## cpm2tm can call 3 different DAG_2_access_genots*

## OT and CBN: DAG_2_access_genots
##    These were the first written. Assume and AND if more than one parent

## OncoBN: if CBN model, same as OT and CBN. If DBN model:
##   DAG_2_access_genots_OR. Derived from previous, uses OR
##   if multiple parents

## HESBCN: if only AND, as for OT/CBN. If only OR, as for OncoBN with DBN.
##  o.w. (XOR or mixtures of any of AND/OR/XOR): cpm2tm_relationships
##    A vector that specifies the relationship if multiple parents
##    is used.

## We try to use the simplest possible for speed and simplicity of code.





## Some of this is coming from
## Cancer_Data_sets/CPM-weighted-paths-biol.R
## in the supplementary material for Diaz-Uriarte and Vasallo.
## We leave in there things we don't really need. Simpler.

## DAG of restrictions (as data frame) ->
##              vector of accessible genotypes and graph of DAG of restrictions
## Under an AND model, such as CBN
## return all the accessible genotypes
##     from a DAG of genes plus the DAG as igraph object
DAG_2_access_genots <- function(x) {
    
    ## minor detail: if x is a sparse adjacency mat.
    ##    from igraph this still works. But don't do that.
    g <- igraph::graph_from_data_frame(x, directed = TRUE)
    
    ## g: an igraph object, with a "Root"
    ## returns a list
    children <- sort(setdiff(V(g)$name, "Root"))
    node_depth <-
        unlist(lapply(children,
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
        tmp_gty_1 <-
            sort(setdiff(names(subcomponent(g, v = names(node_depth)[j],
                                            mode = "in")),
                                  "Root"))
        all_gty[[i]] <- tmp_gty_1

        if (i > 1) {
            to_unite <-  which(
                unlist(lapply(all_gty[1:(i-1)],
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
                             vector(mode = "list",
                                    length = max(length(all_gty),
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



## DAG of restrictions (as data frame) ->
##              vector of accessible genotypes and graph of DAG of restrictions
##  Under an OR model (not XOR, not AND), such as DBN
## return all the accessible genotypes
##     from a DAG of genes plus the DAG as igraph object
DAG_2_access_genots_OR <- function(x) {
    
    ## minor detail: if x is a sparse adjacency mat.
    ##    from igraph this still works. But don't do that.
    g <- igraph::graph_from_data_frame(x, directed = TRUE)
    
    ## g: an igraph object, with a "Root"
    ## returns a list
    children <- sort(setdiff(V(g)$name, "Root"))
    node_depth <-
        unlist(lapply(children,
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
    for (j in seq_along(node_depth)) {
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
                to_unite <-
                    which(unlist(lapply(all_gty[1:(i-1)],
                                        function(z) !all(z %in% tmp_gty_1))))
            } else {
                to_unite <- vector(length = 0)
            }

            if(length(to_unite)) {
                ## we need unique as some sets are the same
                tmp_gty_2 <-
                    unique(lapply(all_gty[to_unite],
                                  function(u) sort(union(u, tmp_gty_1))))
                ## check remaining space in preallocated list.
                ## expand if needed.
                if(length(all_gty) < (i + 1 + length(tmp_gty_2))) {
                    all_gty <- c(all_gty,
                                 vector(mode = "list",
                                        length = max(length(all_gty),
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

## DAG of restrictions (as data frame), options info about gene_relations
##            -> vector of accessible genotypes and graph of DAG of restrictions
##  Under OR, AND, XOR models, as specified by gene_relations
##  Called from cpm2tm_relationships
##   which itself is used for HESBCN
DAG_2_access_genots_relationships <- function(x,
                                                       gene_relations = NULL) {

    ## check if genotype follows the XOR/AND relationship.
    genotype_follows_relationship <- function(genotype) {
        tmp_gene_relations <- gene_relations[genotype]
        for (relationship in c("XOR", "AND")) {
            rel_genes <-
                names(which(tmp_gene_relations == relationship) == TRUE)
            for (rel_idx in rel_genes) {
                parents_of_rel <-
                    adjacent_vertices(g, rel_idx, mode = c("in"))[[rel_idx]]$name
                ## XOR relationship not futfilled
                if (relationship == "XOR" &&
                    (sum(parents_of_rel %in% genotype) > 1)) return(FALSE)
                ## AND relationship not futfilled
                if (relationship == "AND" &&
                    !all(parents_of_rel %in% genotype)) return(FALSE)
            }
        }
        return(TRUE)
    }

    ## g: an igraph object, with a "Root"
    ## returns a list
    g <- igraph::graph_from_data_frame(x, directed = TRUE)

    children <- sort(setdiff(V(g)$name, "Root"))
    node_depth <-
        unlist(lapply(children,
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
    for (j in seq_along(node_depth)) {
        ## FIXME: we did this traversing above. Reuse that
        all_tmp_gty_1 <- lapply(
            all_simple_paths(g, from = "Root", to = names(node_depth)[j],
                             mode = "out"),
            names)

        all_tmp_gty_1 <- lapply(all_tmp_gty_1,
                                function(u) sort(setdiff(u, "Root")))
        for (k in seq_along(all_tmp_gty_1)) {
            tmp_gty_1 <- all_tmp_gty_1[[k]]

            ## Evaluate relationships
            if (genotype_follows_relationship(tmp_gty_1)) {
                all_gty[[i]] <- tmp_gty_1
                i <- i + 1
            }

            if (i > 1) {
                to_unite <-
                    which(unlist(lapply(all_gty[1:(i - 1)],
                                        function(z) !all(z %in% tmp_gty_1))))
            } else {
                to_unite <- vector(length = 0)
            }

            if (length(to_unite)) {
                ## we need unique as some sets are the same
                tmp_gty_2 <-
                    unique(lapply(all_gty[to_unite],
                                  function(u) sort(union(u, tmp_gty_1))))
                ## check remaining space in preallocated list.
                ## expand if needed.
                tmp_gty_2 <-
                    tmp_gty_2[unlist(lapply(tmp_gty_2,
                                            genotype_follows_relationship))]
                if (length(all_gty) < (i + length(tmp_gty_2))) {
                    all_gty <- c(all_gty,
                                 vector(mode = "list",
                                        length = max(length(all_gty),
                                                     2 + length(tmp_gty_2))))
                }
                if (length(tmp_gty_2)){
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


## DAG_2_access_genots_OR(x2)
## DAG_2_access_genots_OR(x3)






## list of accessible genotypes -> adjacency matrix of genotypes (fitness graph)
## return maximally connected fitness graph for a given set of accessible
## genotypes.   This list contains no WT (we add it)
## BEWARE! This is the maximally connected fitness graph. This works with
## CPMs. But this is wrong, if, say, fitnesses are: A = 2, B = 3, AB = 2.5.
## This will place an arrow between B and AB, but there should  be no such edge.
## function genots_2_fgraph_and_trans_mat
## only returns the truly accessible

unrestricted_fitness_graph_sparseM <- function(gacc) {
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
            parents <-
                parents[which(unlist(lapply(parents,
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
    return(adjmat)
}




## output of CPM analysis ->
##             fitness graph (fgraph): just the adj. matrix
##             weighted fitness graph (weighted_fgraph)
##             transition matrix (trans_mat_genots)
##  Assumes an AND relationship

## To be used with output from CBN, MCCBN, OT

## Based on cpm2tm  but only with necessary output
##  for both speed and size and using sparse matrices.
##  No paths are computed, so multiple global max is a moot point

## For CBN and MCCBN the weighted fitness graph is the same
## as the transition rate matrix (with 0 in the diagonal)

## We use the same code for OT and CBN, but lambda and OT_edgeWeight are
## fundamentally different things: the first are rates, the second conditional
## probabilities of an untimed oncogenetic tree.  Computing transition
## probabilities between genotypes (as well as "weighted fitness graphs") is
## abusing what OTs strictly provide.

cpm2tm <- function(x, 
                   names_weights_paths =
                       c("rerun_lambda", ## CBN, from rerun
                         "lambda", ## MCCBN
                         "OT_edgeWeight", ## OT
                         "Lambdas", ## HESBCN
                         "Thetas" ## OncoBN
                         )) {
    if(inherits(x, "try-error") || all(is.na(x)) || is.null(x)) {
        return(list(
            fgraph = "ERROR_CPM_ANALYSIS",
            weighted_fgraph = "ERROR_CPM_ANALYSIS",
            trans_mat_genots = "ERROR_CPM_ANALYSIS"
        ))
    }

    ## Logic of obtaining transition matrices and transition rate matrix
    ##  (weighted fitness graph is trans. rate matrix for CBN and HESBCN)
    ##  - obtain all accessible genotypes and fitness graph (unrestricted
    ##      fitness graph) from CPM
    ##  - from the lambdas/probs. of genes given genes obtain
    ##       probabilities/lambdas of genotypes given genotypes. The call
    ##       to "transition_fg" (where weights is the conditional
    ##       prob./lambda) of descendant gene given parent gene
    ##  - when using CBN, might not be probabilities. Make sure they are
    ##    transition probabilities between genotypes: sweep

    
    ## #############################
    ## Only place in code where dealing with CBN and OT
    ##       can differ from DBN from HESBCN
    ##    If at all possible, always use the simplest processing possibe

    ## sanity check
    if(exists("parent_set", x) && !(exists("Relation", x$edges)))
        stop("x has parent_set but no Relation. Old format?")

    
    ## OT, CBN, and OncoBN if using CBN, HESBCN only Single and AND
    if(!exists("Relation", x$edges) ||
       isTRUE(all(x$edges$Relation %in% c("AND", "Single")))) {
        tmp <- try(DAG_2_access_genots(x$edges[, c("From", "To")]))
    } else if (all(x$edges$Relation %in% c("OR", "Single"))) { ## OncoBN with DBN model, HESBCN if only ORs
        tmp <- try(DAG_2_access_genots_OR(x$edges[, c("From", "To")]))
    } else { ## HESBCN
        ## Use the "parent_set" component, which is returned
        ## from HESBCN itself, not the "Relation" component
        ## which we create from the paernt_set. "Relation" is shown
        ## for a human to easily interpret the $edges data frame
        tmp <-
            try(DAG_2_access_genots_relationships(
                x$edges[, c("From", "To")],
                x$parent_set
            ))
    }

    if(inherits(tmp, "try-error")) {
        stop("Some error here")
    } else {
        accessible_genots <- tmp$accessible_genots
        fgraph <- unrestricted_fitness_graph_sparseM(accessible_genots)
    }

    which_col_weights <- which(colnames(x$edges) %in% names_weights_paths)
    if (length(which_col_weights) > 1) {
        stop("more than one column with weights")
    } else if (length(which_col_weights) == 1) {
        stopifnot(colnames(x$edges)[2] == "To")
        weights <- unique(x$edges[, c(2, which_col_weights)])
        if (any(duplicated(weights[, "To"]))) {
            stop("Different lambda/weight for same destination gene.",
                 " This is not allowed with DAGs")
        }
        rownames(weights) <- weights[, "To"]
        weighted_fgraph <- transition_fg_sparseM(fgraph, weights)
    } else {
        weighted_fgraph <- NA
    }

    ## weighted_fgraph need not have each row sum to 1.  Not if they are lambdas
    ## from CBN, neither do they sum to 1 from some OT models  (see example ex1
    ## in test.trans-rates-f-graphs.R).
    if (length(which_col_weights) == 1) {
        trans_mat_genots <- rowScaleMatrix(weighted_fgraph)
    } else {
        trans_mat_genots <- NA
    }
    
    return(list(
        ## fgraph = fgraph,
        ## again: the weighted_fgraph is the transition rate matrix
        ## for CBN and MCCBN. Not for OT or OncoBN.
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




## fitness graph, weights of probs/lambdas descendant gene given parent gene
##                                           -> weighted fitness graph
##             the fitness graph with weights (relative to each node)
##             of jumping to each of the descendant genotypes
## using sparse matrices
## This is a transition rate matrix (except for the diagonal) if the
## weights are rates, as in CBN
transition_fg_sparseM <- function(x, weights) { 
    ## pos_do <- which(x == 1, arr.ind = TRUE)
    pos_do <- as.matrix(summary(x)[, c("i", "j")])
    wfg <- x
    wfg[] <- 0
    tmp <- unlist(
        Map(function(r, c)
             get_single_lambda(x, r, c, weights),
            pos_do[, 1], pos_do[, 2]))
    wfg[pos_do] <- tmp
    return(wfg)
}




## Main function. data frame or matrix -> output
evam <- function(x,
                 methods = c("CBN", "OT", "HESBCN", "MHN", "OncoBN"),
                 max_cols = 15,
                 cbn_cores = 1,
                 cbn_init_poset = "OT",
                 hesbcn_steps = 100000,
                 hesbcn_seed = NULL,
                 mhn_lambda = 1/nrow(x),
                 oncobn_model = "DBN",
                 oncobn_algorithm = "DP",
                 oncobn_epsilon = min(colMeans(x)/2),
                 oncobn_silent = TRUE,
                 ot_with_errors_dist_ot = TRUE
                 ) {

    if(!(cbn_init_poset %in% c("OT", "linear")))
        stop("cbn_init_poset must be one of OT or linear. ",
             " Custom not allowed in call from evam.")
    
    if ("MCCBN" %in% methods) {
        do_MCCBN <- TRUE
    } else {
        do_MCCBN <- FALSE
    }

    if ("OncoBN" %in% methods) {
        do_OncoBN <- TRUE
    } else {
        do_OncoBN <- FALSE
    }

    do_HyperTraPS <- FALSE


    ##########      Preprocessing: common to all methods
    x <- df_2_mat_integer(x)
    xoriginal <- x
    
    x <- add_pseudosamples(x, n00 = "auto3")
    ## remove.constant makes no difference IFF we add pseudosamples, as
    ## there can be no constant column when we add pseudosamples
    x <- pre_process(x, remove.constant = FALSE,
                     min.freq = 0, max.cols = max_cols)

    if(ncol(x) < 2) {
        warning("Fewer than 2 columns in data set")
        return(NA)
    }

    
    ## ## Not using HyperTraPS for now. A placeholder
    ## if(do_HyperTraPS) {
    ##     ##     message("Doing HyperTraps")
    ##     ##     message("By default we run it here with dry_run = TRUE.
    ##     ##     HyperTraPS takes a long time and I do no want to block the script.
    ##     ## ")
    ##     ## FIXME: HT_folder undefined
    ##     HT_folder <- NULL
    ##     time_HyperTraPS <- system.time(
    ##         out_HyperTraPS <-
    ##             do_HyperTraPS(x, tmp_folder = HT_folder,
    ##                           dry_run = TRUE, plot = FALSE))["elapsed"]
    ##     message("time HyperTraPS = ", time_HyperTraPS)
    ## } 
    ## if (do_HyperTraPS) {
    ##     ## FIXME: this will break, as pre_trans_mat_new_CPMS$HyperTraPS
    ##     ## has never been created
    ##     pre_trans_mat_others["HyperTraPS"] <-
    ##         list(pre_trans_mat_new_CPMS$HyperTraPS)
    ## } 


    ## ######################################################################
    ##      Run MHN
    ##      MHN is fully self-contained: CPM, transition matrix, etc.
    ## ######################################################################
    
    message("Doing MHN")
    time_MHN <- system.time(
        MHN_out <- do_MHN2(x, lambda = mhn_lambda))["elapsed"]
    message("time MHN = ", time_MHN)


    ## ####################################################################
    ##    Run each one of the remaining CPMs
    ## ####################################################################
    
    message("Doing HESBCN")
    time_hesbcn <- system.time(
        HESBCN_out <- do_HESBCN(x,
                                n_steps = hesbcn_steps,
                                seed = hesbcn_seed))["elapsed"]
    message("time HESBCN = ", time_hesbcn)


    
    message("Doing OT")
    time_ot <-
        system.time(
            OT_out <- try(
                suppressMessages(
                    ot_proc(x,
                            nboot = 0,
                            distribution.oncotree = TRUE,
                            with_errors_dist_ot = ot_with_errors_dist_ot))))["elapsed"]
    message("time OT = ", time_ot)

    message("Doing CBN")
    time_cbn_ot <- system.time(
           CBN_out <- try(cbn_proc(x,
                                   addname = "tmpo",
                                   init.poset = cbn_init_poset,
                                   nboot = 0,
                                   parall = TRUE,
                                   cores = cbn_cores)))["elapsed"]
    message("time CBN = ", time_cbn_ot)

    if(do_OncoBN) {
        message("Doing OncoBN\n\n")
        time_dbn <- system.time(
            OncoBN_out <- do_OncoBN(x,
                                    model = oncobn_model,
                                    algorithm = oncobn_algorithm,
                                    epsilon = oncobn_epsilon,
                                    silent = oncobn_silent))["elapsed"]
        message("time OncoBN = ", time_dbn)
    } 
    
    if(do_MCCBN) {
        stop("MC-CBN disabled")
        ## if( !requireNamespace("mccbn", quietly = TRUE)) {
        ##     warning("MC-CBN (mccbn) no installed. Not running MC-CBN")
        ## } else {
        ##     stop("MC-CBN disabled now")
        ##     ## message("Doing MC-CBN")
        ##     ## time_mccbn <-
        ##     ##     system.time(MCCBN_out <- try(mccbn_proc(x)))["elapsed"]
        ##     ## message("time MC-CBN = ", time_mccbn)
        ## }
    }

    ## ######################################################################
    ##   From CPM output, obtain weighted fitness graph and transition matrix
    ##   between genotypes.
    ## ######################################################################
    
    OT_fg_tm <- cpm2tm(OT_out)
    CBN_fg_tm <- cpm2tm(CBN_out)
    HESBCN_fg_tm <- cpm2tm(HESBCN_out)
    if(do_MCCBN)
        MCCBN_fg_tm <- cpm2tm(MCCBN_out)
    if(do_OncoBN) 
        OncoBN_fg_tm <- cpm2tm(OncoBN_out)
    
    
    ## ######################################################################
    ##    For methods for which it makes sense, obtain time discretized
    ##    transition rate matrix
    ##    We use the same function for all cases.
    ## ######################################################################
    
    CBN_td_trans_mat <-  trans_rate_to_trans_mat(CBN_fg_tm$weighted_fgraph,
                                       method = "uniformization")
    if(do_MCCBN)
        MCCBN_td_trans_mat <-  trans_rate_to_trans_mat(MCCBN_fg_tm$weighted_fgraph,
                                       method = "uniformization")
    HESBCN_td_trans_mat <-  trans_rate_to_trans_mat(HESBCN_fg_tm$weighted_fgraph,
                                       method = "uniformization")
    

    ## Return list always has the same elements
    ## Set to NA those for not run models.
    if (!do_MCCBN) {
        MCCBN_out <- list(edges = NA)
        MCCBN_fg_tm <- list(weighted_fgraph = NA, trans_mat_genots = NA)
        MCCBN_td_trans_mat <- NA
    }
    if (!do_OncoBN) {
        OncoBN_out <- list(edges = NA, likelihood = NA)
        OncoBN_fg_tm <- list(weighted_fgraph = NA, trans_mat_genots = NA)
    }

    HyperTraPS_model <- NA
    HyperTraPS_trans_rate_mat <- NA
    HyperTraPS_trans_mat <- NA
    HyperTraPS_td_trans_mat <- NA

    ## Future: Getting all paths to global maximum
    ## Simply call do_weighted_paths_to_max on a list of
    ## transition matrices and a pre-created paths_to_max
    ## Make sure we are not repeating expensive operations
    ## When testing, compare against cpm2tm
    ## for OT and CBN

    ## f_graph: remember this is the transition rate matrix
    ## for CBN, MCCBN, HESBCN
    ## For OT ... well, it is something else, but not really probabilities
    ## of anything.
    return(list(
        OT_model = OT_out$edges,
        OT_f_graph = OT_fg_tm$weighted_fgraph,
        OT_trans_mat = OT_fg_tm$trans_mat_genots,
        OT_genots_predicted = OT_out$genots_predicted,

        CBN_model = CBN_out$edges,
        CBN_trans_rate_mat = CBN_fg_tm$weighted_fgraph,
        CBN_trans_mat = CBN_fg_tm$trans_mat_genots,
        CBN_td_trans_mat = CBN_td_trans_mat,

        MCCBN_model = MCCBN_out$edges,
        MCCBN_trans_rate_mat = MCCBN_fg_tm$weighted_fgraph,
        MCCBN_trans_mat = MCCBN_fg_tm$trans_mat_genots,
        MCCBN_td_trans_mat = MCCBN_td_trans_mat,

        MHN_theta = MHN_out$theta,
        MHN_exp_theta = exp(MHN_out$theta),
        MHN_trans_rate_mat = MHN_out$transitionRateMatrix,
        MHN_trans_mat = MHN_out$transitionMatrixCompExp,
        MHN_td_trans_mat = MHN_out$transitionMatrixTimeDiscretized,

        OncoBN_model = OncoBN_out$edges,
        OncoBN_likelihood = OncoBN_out$likelihood, 
        OncoBN_f_graph = OncoBN_fg_tm$weighted_fgraph, 
        OncoBN_trans_mat = OncoBN_fg_tm$trans_mat_genots,
        OncoBN_genots_predicted = OncoBN_out$genots_predicted,
        OncoBN_fitted_model = OncoBN_out$model,
        OncoBN_epsilon = OncoBN_out$epsilon,

        HESBCN_model = HESBCN_out$edges,
        HESBCN_parent_set = HESBCN_out$parent_set,
        HESBCN_trans_rate_mat = HESBCN_fg_tm$weighted_fgraph,
        HESBCN_trans_mat = HESBCN_fg_tm$trans_mat_genots,
        HESBCN_td_trans_mat = HESBCN_td_trans_mat,

        HyperTraPS_model = HyperTraPS_model,
        HyperTraPS_trans_rate_mat = HyperTraPS_trans_rate_mat,
        HyperTraPS_trans_mat = HyperTraPS_trans_mat,
        HyperTraPS_td_trans_mat = HyperTraPS_td_trans_mat,
        original_data = xoriginal,
        analyzed_data = x
         ))
}





#####################################################################
#####################################################################
##                Getting an idea of sizes of data if
##                we did not use sparse matrices
##
#####################################################################
#####################################################################

# if(FALSE) {
# ## 8 bytes per float as can be checked doing
# u <- runif((2^10) * (2^10))
# print(object.size(u), units = "b")
# 8 * length(u)

# uu <-  runif((2^15) * (2^15))
# length(uu)
# print(object.size(uu), units = "b")
# print(object.size(uu), units = "GB") ## 8 GB

# length(uu) * 8 /(1024 * 1024 * 1024)

# ## so an object for 2^20 * 2^20 would take

# ((2^20)^2) * 8 / (1024 * 1024 * 1024)

# ## or 8192 GB which is what R complaints about if I try this
# uuu <-  runif((2^20) * (2^20))

# ## Yes, there are no memory limits in the machines
# ## install.packages("devtools", dependencies = TRUE)
# ##devtools::install_github("krlmlr/ulimit")
# library(ulimit)
# ulimit::memory_limit()

# }



######################################################################
######################################################################
##
##      mapply, OMP threads, etc
##
######################################################################
######################################################################


## MHN does not run if we use mclapply after setting the threads for OMP >
## 1

## all_data_out <- mclapply(all_data, all_methods_2_trans_mat, mc.cores
## = 36)

## Yes, you can use mclapply, just make sure to set the threads of OMP to 1
## The following shows it. It just did not make any sense above.
# if(FALSE) {

#     lapply(all_data, dim)

#     RhpcBLASctl::omp_set_num_threads(36)

#     system.time(cucu3 <- lapply(all_data[c(3, 4, 20, 21:24)], do_MHN))

#     RhpcBLASctl::omp_set_num_threads(1)
#     system.time(cucu4 <- mclapply(all_data[c(3, 4, 20, 21:24)], do_MHN,
#                                   mc.cores = 36))
#     stopifnot(identical(cucu3, cucu4))
# }


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











## ## output of CPM analysis ->
## ##             fitness graph (fgraph): just the adj. matrix
## ##             weighted fitness graph (weighted_fgraph)
## ##             transition matrix (trans_mat_genots)
## ## Like cpm2tm but for OR relationships

## ## To be used with output from DBN

## ## DBN seems similar to OT: these are conditional probabilities, not transition
## ## rates. Beware of interpretations. See comments in
## ## cpm2tm
## cpm2tm_OR <- function(data,
##                                                     parameter_column_name = c("Thetas")) {
    
##     tmp <-
##         try(DAG_2_access_genots_OR(data$edges[, c("From", "To")]))

##     if (inherits(tmp, "try-error")) {
##         stop("Some error here")
##     } else {
##         accessible_genots <- tmp$accessible_genots
##         fgraph <- unrestricted_fitness_graph_sparseM(accessible_genots)
##     }

##     which_col_weights <-
##         which(colnames(data$edges) %in% parameter_column_name)
    
##     if (is.null(data$edges[,which_col_weights])){
##         stop("No such column")
##     }

##     weights <- unique(data$edges[, c(2, which_col_weights)])
##     rownames(weights) <- weights[, "To"]
##     weighted_fgraph <- transition_fg_sparseM(fgraph, weights)
##     trans_mat_genots <- rowScaleMatrix(weighted_fgraph)

##     return(list(
##         fgraph = fgraph,
##         ## Seems similar to OT. This ain't the transition rate matrix.
##         weighted_fgraph = weighted_fgraph,
##         trans_mat_genots = trans_mat_genots
##     ))
## }


## ## output of CPM analysis ->
## ##             fitness graph (fgraph): just the adj. matrix
## ##             weighted fitness graph (weighted_fgraph)
## ##             transition matrix (trans_mat_genots)
## ## To be used with output from HESBCN

## ## The weighted fitness graph is the same
## ## as the transition rate matrix (with 0 in the diagonal)

## cpm2tm_relationships <-
##     function(data,
##              parameter_column_name = c("Lambdas")) {             

##         ## Use the "parent_set" component, which is returned
##         ## from HESBCN itself, not the "Relation" component
##         ## which we create from the paernt_set. "Relation" is shown
##         ## for a human to easily interpret the $edges data frame
##         tmp <-
##             try(DAG_2_access_genots_relationships(
##                 data$edges[, c("From", "To")],
##                 data$parent_set
##             ))

##         if (inherits(tmp, "try-error")) {
##             stop("Some error here")
##         } else {
##             accessible_genots <- tmp$accessible_genots
##             fgraph <- unrestricted_fitness_graph_sparseM(accessible_genots)
##         }

##         which_col_weights <-
##             which(colnames(data$edges) %in% parameter_column_name)
##         if (is.null(data$edges[, which_col_weights])) {
##             stop("No such column")
##         }

##         weights <- unique(data$edges[, c(2, which_col_weights)])
##         rownames(weights) <- weights[, "To"]
##         weighted_fgraph <- transition_fg_sparseM(fgraph, weights)
##         trans_mat_genots <- rowScaleMatrix(weighted_fgraph)

##         return(list(
##             fgraph = fgraph,
##             ## weighted_fgraph is the transition rate matrix
##             weighted_fgraph = weighted_fgraph,
##             trans_mat_genots = trans_mat_genots
##         ))
## }
