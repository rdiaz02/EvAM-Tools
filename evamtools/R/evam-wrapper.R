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
    children <- evam_string_sort(setdiff(V(g)$name, "Root"))
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
            evam_string_sort(setdiff(names(subcomponent(g,
                                                        v = names(node_depth)[j],
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
                                       function(u) evam_string_sort(union(u, tmp_gty_1))))
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
    return(list(accessible_genots = all_gty
              ##, graph = g
                ))
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
    children <- evam_string_sort(setdiff(V(g)$name, "Root"))
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
                                function(u) evam_string_sort(setdiff(u, "Root")))
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
                                  function(u) evam_string_sort(union(u, tmp_gty_1))))
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
    return(list(accessible_genots = all_gty
              ## , graph = g
                ))
}

## DAG of restrictions (as data frame), options info about gene_relations
##            -> vector of accessible genotypes and graph of DAG of restrictions
##  Under OR, AND, XOR models, as specified by gene_relations
##  Called from cpm2tm_relationships
##   which itself is used for HESBCN

## We could just rm the "gene_relations" and create it as we do inside
## Left as a consistency check
DAG_2_access_genots_relationships <- function(x,
                                              gene_relations = NULL) {

    if (is.null(gene_relations)) {
        message("gene_relations (parent_set) NULL. ",
                "Creating it from edges component.")
        gene_relations <- parent_set_from_edges(x)
    } else {
        ## A consistency check
        ng <- evam_string_sort(names(gene_relations))
        parent_set_input <- gene_relations[ng]
        parent_set_edges <- parent_set_from_edges(x)[ng]
        stopifnot(identical(parent_set_input, parent_set_edges))
    }

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

    children <- evam_string_sort(setdiff(V(g)$name, "Root"))
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
                                function(u) evam_string_sort(setdiff(u, "Root")))
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
                                  function(u) evam_string_sort(union(u, tmp_gty_1))))
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
    return(list(accessible_genots = all_gty
              ## , graph = g
                ))
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

## In gacc, each genotype is a vector of genes, and the vector is already
## sorted. This ensures that the naming of genotypes is consistent later.
unrestricted_fitness_graph_sparseM <- function(gacc) {
    gs <- unlist(lapply(gacc, function(g) paste0(g, collapse = ", ")))
    gs <- c("WT", gs)
    nmut <- c(0, vapply(gacc, length, 1))
    gs <- gs[order(nmut, gs)]
    
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


## Based on cpm_access_genots_paths_w  but only with necessary output
##  for both speed and size and using sparse matrices,
##  and allowing OR, XOR, and mixes of OR, AND, XOR.
##  No paths are computed, so multiple global max is a moot point

## For CBN and MCCBN and HESBCN the weighted fitness graph is the same
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
                         "theta" ## OncoBN
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
    if (!exists("Relation", x$edges) ||
       isTRUE(all(x$edges$Relation %in% c("AND", "Single")))) {
        tmp <- try(DAG_2_access_genots(x$edges[, c("From", "To")]))
    } else if (all(x$edges$Relation %in% c("OR", "Single"))) {
        ## OncoBN with DBN model, HESBCN if only ORs
        tmp <- try(DAG_2_access_genots_OR(x$edges[, c("From", "To")]))
    } else { ## HESBCN
        ## Use the "parent_set" component, which is returned
        ## from HESBCN itself, not the "Relation" component
        ## which we create from the parent_set. "Relation" is shown
        ## for a human to easily interpret the $edges data frame.
        ## We check for consistency anyway.
        tmp <-
            try(DAG_2_access_genots_relationships(
                x$edges[, c("From", "To", "Relation")],
                x$parent_set
            ))
    }

    if (inherits(tmp, "try-error")) {
        stop("Some error here: ", tmp)
    } else {
        fgraph <- unrestricted_fitness_graph_sparseM(tmp$accessible_genots)
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

    ## OncoBN can return edges with theta 0. So rm the inaccessible
    ## genotypes. Do not remove by predicted genotype freq. = 0.
    ## That would be wrong: with a theta of 1, the parent genotype
    ## can have a predicted freq. of 0, and yet be part of the path.
    if ("theta" %in% colnames(x$edges)) {
        weighted_fgraph <- adjm_rm_no_access(weighted_fgraph)
    }

    ## Check
    if (sum(colSums(weighted_fgraph) == 0) > 1) ## WT always has colSum = 0
        warning("weighted_fgraph contains unreachable destinations")

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
                 methods = c("CBN", "OT", "HESBCN", "MHN", "OncoBN", "MCCBN", "HyperHMM"),
                 max_cols = 15,
                 cores = detectCores(),
                 paths_max = FALSE,
                 mhn_opts = list(lambda = 1/nrow(x),
                                 omp_threads = ifelse(cores > 1, 1, detectCores())
                                ),
                 hyperhmm_opts = list(precursors = NA, # precursors (e.g. ancestors) -- blank for cross-sectional
                                      nboot = 1, # number of boostrap resamples
                                      random.walkers = 0, # run random walkers for each resample? 0 no, 1 yes
                                      label = "label", # label for file I/O
                                      simulate = TRUE, # actually run HyperHMM? If not, try and pull precomputed output using "label"
                                      fork = FALSE # if we are running it, fork to a new process, or keep linear?
                                      ),
                 ot_opts = list(with_errors_dist_ot = TRUE),
                 cbn_opts = list(
                     omp_threads = 1,
                     init_poset = "OT"
                 ),
                 hesbcn_opts = list(
                     MCMC_iter = 100000,
                     seed = NULL,
                     reg = c("bic", "aic", "loglik"),
                     silent = TRUE                     
                 ),
                 oncobn_opts = list(
                     model = "DBN",
                     algorithm = "DP",
                     k = 3,
                     epsilon = min(colMeans(x)/2),
                     silent = TRUE
                 ),
                 mccbn_opts = list(
                     model = "OT-CBN",
                     tmp_dir = NULL,
                     addname = NULL,
                     silent = TRUE,
                     L = 100,
                     sampling = c("forward", "add-remove", "backward", "bernoulli", "pool"),
                     max.iter = 100L,
                     update.step.size = 20L,
                     tol = 0.001,
                     max.lambda.val = 1e+06,
                     T0 = 50,
                     adap.rate = 0.3,
                     acceptance.rate = NULL,
                     step.size = NULL,
                     max.iter.asa = 10000L,
                     neighborhood.dist = 1L,
                     adaptive = TRUE,
                     thrds = 1L,
                     verbose = FALSE,
                     seed = NULL)
                 ) {

    ## Sanity check of gene names
    gn_comma <- stringi::stri_count_fixed(colnames(x), ",")
    if(any(gn_comma))
        stop("At least one of your gene names has a comma. That is not allowed")

    gn_backslash <- stringi::stri_count_fixed(colnames(x), "\\") 
    if(any(gn_backslash))
        stop("At least one of your gene names has a backslash. That is not allowed")

    gn_space <- stringi::stri_count_regex(colnames(x), "[\\s]") 
    if(any(gn_space))
        stop("At least one of your gene names has a space. That is not allowed")
    
    if(any(colnames(x) == "WT"))
        stop("One of your genes is called WT. That is not allowed")

    
    ## Sanity checks of methods
    methods <- unique(methods)
    accepted_methods <- c("OT", "OncoBN", "CBN", "MCCBN", "MHN", "HESBCN", "HyperHMM")
    not_valid_methods <- which(!(methods %in% accepted_methods))
    if (length(not_valid_methods)) {
        warning("Method(s) ",
                paste(methods[not_valid_methods], sep = ", ", collapse = ", "),
                " not among the available methods.",
                " Ignoring the invalid method.")
        methods <- methods[-not_valid_methods]
    }
    
    if ("MCCBN" %in% methods) {
        MCCBN_INSTALLED <- requireNamespace("mccbn", quietly = TRUE)
        if (!MCCBN_INSTALLED) {
            warning("MCCBN method requested, but mccbn packaged not installed. ",
                    "Removing MCCBN from list of requested methods.")
            methods <- setdiff(methods, "MCCBN")
        }
    }
    
    if (length(methods) == 0) stop("No valid methods given.")

    ## ########      Preprocessing: common to all methods
    x <- df_2_mat_integer(x)
    xoriginal <- x
    
    x <- add_pseudosamples(x)
    ## remove.constant makes no difference IFF we add pseudosamples, as
    ## there can be no constant column when we add pseudosamples
    x <- pre_process(x, remove.constant = FALSE,
                     min.freq = 0, max.cols = max_cols)

    if (ncol(x) < 2) {
        stop("Fewer than 2 columns in the data. ",
             "There must be at least two genes ",
             "and two different genotypes to run evam ",
             "(and remember that any genes that are ",
             "completely aliased, i.e., indistinguishable, ",
             "because they have identical patterns ",
             "---identical columns in the data matrix--- ",
             "are regarded as a single gene).")
    }
    
    
    ## ############################################################
    ##
    ##       Dealing with default arguments
    ##
    ##  The "d_*" are the list of default arguments.
    ##  The user can pass, in the call, just a subset of the
    ##  options, and the rest will be taken from the "d_*" below.
    ##  If the user passes something not in "d_*", give warning.
    ##  
    ## ############################################################
    d_mhn_opts <- list(lambda = 1/nrow(x),
                       omp_threads = ifelse(cores > 1, 1, detectCores())
                       )
    d_cbn_opts <- list(
        omp_threads = 1,
        init_poset = "OT"
    )
    d_hyperhmm_opts = list(
      precursors = NA, 
      nboot = 1,
      random.walkers = 0, 
      label = "label", 
      simulate = TRUE, 
      fork = FALSE
    )
    d_hesbcn_opts <- list(
        MCMC_iter = 100000,
        seed = NULL,
        reg = c("bic", "aic", "loglik"),
        silent = TRUE
    )
    d_oncobn_opts <- list(
        model = "DBN",
        algorithm = "DP",
        k = 3,
        epsilon = min(colMeans(x)/2),
        silent = TRUE
    )
    d_mccbn_opts <- list(
        model = "OT-CBN",
        tmp_dir = NULL,
        addname = NULL,
        silent = TRUE,
        L = 100,
        sampling = c("forward", "add-remove", "backward", "bernoulli", "pool"),
        max.iter = 100L,
        update.step.size = 20L,
        tol = 0.001,
        max.lambda.val = 1e+06,
        T0 = 50,
        adap.rate = 0.3,
        acceptance.rate = NULL,
        step.size = NULL,
        max.iter.asa = 10000L,
        neighborhood.dist = 1L,
        adaptive = TRUE,
        thrds = 1L,
        verbose = FALSE,
        seed = NULL)

    ## Not needed, but make sure we keep names distinct
    mhn_opts_2 <- fill_args_default(mhn_opts, d_mhn_opts)
    hyperhmm_opts_2 <- fill_args_default(hyperhmm_opts, d_hyperhmm_opts)
    cbn_opts_2 <- fill_args_default(cbn_opts, d_cbn_opts)
    hesbcn_opts_2 <- fill_args_default(hesbcn_opts, d_hesbcn_opts)
    oncobn_opts_2 <- fill_args_default(oncobn_opts, d_oncobn_opts)
    mccbn_opts_2 <- fill_args_default(mccbn_opts, d_mccbn_opts)

    ## rm to avoid confusion, though not needed
    rm(cbn_opts, hesbcn_opts, oncobn_opts, mccbn_opts, mhn_opts, hyperhmm_opts)
    rm(d_cbn_opts, d_hesbcn_opts, d_oncobn_opts, d_mccbn_opts, d_mhn_opts, d_hyperhmm_opts)

    if (!(cbn_opts_2$init_poset %in% c("OT", "linear")))
        stop("CBN's init_poset must be one of OT or linear. ",
             " Custom not allowed in call from evam.")

    if ("MCCBN" %in% methods) {
        stopifnot(mccbn_opts_2$model %in% c("OT-CBN", "H-CBN2"))
    }

    if ("OncoBN" %in% methods) {
        stopifnot(oncobn_opts_2$model %in% c("DBN", "CBN"))
    }
    
    if ("HyperHMM" %in% methods) {
      # possible names of executable
      cmds <- c("hyperhmm.ce", "hyperhmm.exe", "hyperhmm")
      #see if any are here
      found <- FALSE
      commandname <- ""
      for(cmd in cmds) {
        commandpath <- Sys.which(cmd)
        if(nchar(commandpath) > 0) {
          commandname <- commandpath
          found <- TRUE
        }
      }
      if (found == FALSE) { message ("Couldn't find HyperHMM executable in the PATH") }
      if(ncol(xoriginal) > 15) { message("In HyperHMM model you are welcome to try > 15 features but this may use a lot of memory.") }
    }
    
    
    ## ######################################################################
    ##   Big function so we can parallelize the calls
    ## Run inside evam. Take all additional args from envir of evam
    do_method <- function(method) {
        if (method == "MHN") {
            RhpcBLASctl::omp_set_num_threads(mhn_opts_2$omp_threads)
            time_out <- system.time({
                out <- do_MHN2(x, lambda = mhn_opts_2$lambda)
                out <- c(out,
                     predicted_genotype_freqs =
                         list(probs_from_trm(out$transitionRateMatrix)))
                })["elapsed"]
        } else if (method == "HyperHMM") {
            time_out <- system.time({
              out <- do_HyperHMM(xoriginal, commandname, precursors = NA, nboot = 1, random.walkers = 0,
                                 label = "label", simulate = TRUE,fork = FALSE)
              
            })
        } else if (method == "HESBCN") {
            time_out <- system.time({
                out <- do_HESBCN(x,
                                 MCMC_iter = hesbcn_opts_2$MCMC_iter,
                                 seed = hesbcn_opts_2$seed,
                                 silent = hesbcn_opts_2$silent,
                                 reg = hesbcn_opts_2$reg)
                out <- c(out, cpm2tm(out))
                out <- c(out,
                         td_trans_mat =
                             trans_rate_to_trans_mat(out[["weighted_fgraph"]],
                                                     method = "uniformization"))
                out <- c(out,
                         predicted_genotype_freqs =
                             list(probs_from_trm(out$weighted_fgraph)))
            })["elapsed"]
        } else if (method == "CBN") {
            time_out <- system.time({
                out <- try(cbn_proc(x,
                                    addname = "tmpo",
                                    init.poset = cbn_opts_2$init_poset,
                                    nboot = 0,
                                    parall = TRUE,
                                    omp_threads = cbn_opts_2$omp_threads))
            out <- c(out, cpm2tm(out))
            out <- c(out,
                     td_trans_mat =
                         trans_rate_to_trans_mat(out[["weighted_fgraph"]],
                                                 method = "uniformization"))
            out <- c(out,
                     predicted_genotype_freqs =
                         list(probs_from_trm(out$weighted_fgraph)))
                })["elapsed"]
        } else if (method == "MCCBN") {
            if(mccbn_opts_2$model == "OT-CBN") 
                time_out <-
                    system.time(out <- try(do_MCCBN_OT_CBN(x)))["elapsed"]
            else if(mccbn_opts_2$model == "H-CBN2") {
                mccbn_hcbn2_opts_2 <- mccbn_opts_2
                mccbn_hcbn2_opts_2$model <- NULL
                time_out <-
                    system.time(
                        out <- try(
                            do_MCCBN_HCBN2(x,
                                           mccbn_hcbn2_opts = mccbn_hcbn2_opts_2
                                           )))["elapsed"]
            }
            time_out2 <- system.time({
            out <- c(out, cpm2tm(out))
            out <- c(out,
                     td_trans_mat =
                         trans_rate_to_trans_mat(out[["weighted_fgraph"]],
                                                 method = "uniformization"))
            out <- c(out,
                     predicted_genotype_freqs =
                         list(probs_from_trm(out$weighted_fgraph)))
            })["elapsed"]
            time_out <- time_out + time_out2
        } else if (method == "OT") {
            time_out <- system.time({
                out <- try(
                    suppressMessages(
                        ot_proc(x,
                                nboot = 0,
                                distribution.oncotree = TRUE,
                                with_errors_dist_ot = ot_opts$with_errors_dist_ot)))
                out <- c(out, cpm2tm(out))
            })["elapsed"]
        } else if (method == "OncoBN") {
            time_out <- system.time({
                out <- do_OncoBN(x,
                                 model = oncobn_opts_2$model,
                                 algorithm = oncobn_opts_2$algorithm,
                                 k = oncobn_opts_2$k,
                                 epsilon = oncobn_opts_2$epsilon,
                                 silent = oncobn_opts_2$silent)
                out <- c(out, cpm2tm(out))
                })["elapsed"]
        }
        message(paste0("time ", method, ": ", round(time_out, 3)))
        out
    }

    all_out <- mclapply(methods, do_method, mc.cores = cores)
    names(all_out) <- methods

    
    ## f_graph: remember this is the transition rate matrix
    ## for CBN, MCCBN, HESBCN, HyperHMM
    ## For OT ... well, it is something else, but not really probabilities
    ## of anything.

    ## To avoid repeated code
    get_output <- function(method, component) {
        if (!exists(method, all_out)) return(NA)
        if (!exists(component, all_out[[method]])) return(NA)
        return(all_out[[method]][[component]])
    }
    
    ## To avoid repeating code
    get_paths_max <- function(method) {
        if (paths_max) {
            trans_mat_name <- ifelse(method == "MHN",
                                     "transitionMatrixCompExp",
                                     "trans_mat_genots")
            trans_mat <- get_output(method, trans_mat_name)
            if ((length(trans_mat) == 1) && is.na(trans_mat)) return(NA)
            return(paths_probs_2_df(trans_mat_2_paths_probs(trans_mat),
                                    order = "prob"))
        } else {
            return(NA)
        }
    }
    return(list(
        OT_model = get_output("OT", "edges"),
        OT_f_graph = get_output("OT", "weighted_fgraph"),
        OT_trans_mat = get_output("OT", "trans_mat_genots"),
        OT_predicted_genotype_freqs = get_output("OT",
                                                 "predicted_genotype_freqs"),
        OT_eps = get_output("OT", "eps"),
        OT_fit = get_output("OT", "ot_fit"),
        OT_paths_max = get_paths_max("OT"),
       
        CBN_model = get_output("CBN", "edges"),
        CBN_trans_rate_mat = get_output("CBN", "weighted_fgraph"),
        CBN_trans_mat = get_output("CBN", "trans_mat_genots"),
        CBN_td_trans_mat = get_output("CBN", "td_trans_mat"),
        CBN_predicted_genotype_freqs = get_output("CBN",
                                                  "predicted_genotype_freqs"),
        CBN_paths_max = get_paths_max("CBN"),
        
        MCCBN_model = get_output("MCCBN", "edges"),
        MCCBN_trans_rate_mat = get_output("MCCBN", "weighted_fgraph"),
        MCCBN_trans_mat = get_output("MCCBN", "trans_mat_genots"),
        MCCBN_td_trans_mat = get_output("CBN", "td_trans_mat"),
        MCCBN_predicted_genotype_freqs = get_output("MCCBN",
                                                    "predicted_genotype_freqs"),
        MCCBN_paths_max = get_paths_max("MCCBN"),
        
        MHN_theta = get_output("MHN", "theta"),
        MHN_trans_rate_mat = get_output("MHN", "transitionRateMatrix"),
        MHN_trans_mat = get_output("MHN", "transitionMatrixCompExp"),
        MHN_td_trans_mat = get_output("MHN", "transitionMatrixTimeDiscretized"),
        MHN_exp_theta = exp(get_output("MHN", "theta")),
        MHN_predicted_genotype_freqs = get_output("MHN",
                                                  "predicted_genotype_freqs"),
        MHN_paths_max = get_paths_max("MHN"),

        
        HyperHMM_stats.df = get_output('HyperHMM', 'stats.df'),
        HyperHMM_model = get_output('HyperHMM', 'transitions'),
        HyperHMM_features = get_output('HyperHMM', 'features'),
        HyperHMM_viz.tl = get_output('HyperHMM', 'viz.tl'),
        HyperHMM_trans_mat = get_output("HyperHMM", "trans_mat"),

        #HyperHMM = get_paths_max('HyperHMM'),

        OncoBN_model = get_output("OncoBN", "edges"),
        OncoBN_likelihood = get_output("OncoBN", "likelihood"),
        OncoBN_f_graph = get_output("OncoBN", "weighted_fgraph"), 
        OncoBN_trans_mat = get_output("OncoBN", "trans_mat_genots"),
        OncoBN_predicted_genotype_freqs = get_output("OncoBN",
                                                     "predicted_genotype_freqs"),
        OncoBN_fitted_model = get_output("OncoBN", "model"),
        OncoBN_epsilon = get_output("OncoBN", "epsilon"),
        OncoBN_parent_set = get_output("OncoBN", "parent_set"),
        OncoBN_fit = get_output("OncoBN", "fit"),
        OncoBN_paths_max = get_paths_max("OncoBN"),
        
        HESBCN_model = get_output("HESBCN", "edges"),
        HESBCN_parent_set = get_output("HESBCN", "parent_set"),
        HESBCN_trans_rate_mat = get_output("HESBCN", "weighted_fgraph"),
        HESBCN_trans_mat = get_output("HESBCN", "trans_mat_genots"),
        HESBCN_td_trans_mat = get_output("HESBCN", "td_trans_mat"),
        HESBCN_predicted_genotype_freqs = get_output("HESBCN",
                                                     "predicted_genotype_freqs"),
        HESBCN_command = get_output("HESBCN", "command"),
        HESBCN_paths_max = get_paths_max("HESBCN"),
        
        original_data = xoriginal,
        analyzed_data = x,
        genotype_id_ordered =
            stats::setNames(1:(2^ncol(x)),
                            genotypes_standard_order(colnames(x))),
        all_options = list(
            mhn_opts = mhn_opts_2,
            hyperhmm_opts = hyperhmm_opts_2,
            ot_opts = ot_opts,
            cbn_opts = cbn_opts_2,
            hesbcn_opts = hesbcn_opts_2,
            oncobn_opts = oncobn_opts_2,
            mccbn_opts = mccbn_opts_2)
    
    )
    )
}

## transition matrix to paths to maximum
## as list with two components: igraph.vs
## and log probs
trans_mat_2_paths_probs <- function(trans_mat) {
    g <- igraph::graph_from_adjacency_matrix(trans_mat,
                                             weighted = TRUE,
                                             mode = "directed")
    rp <- rank_paths(g, log_weights = TRUE)
    stopifnot(all.equal(sum(exp(rp$weights)), 1,
                        tolerance = 100 * sqrt(.Machine$double.eps)))
    return(rp)
}
    
## Take an object returned by trans_mat_2_paths_probs or rank_paths
## (list of paths and their weights)
## and return a data frame with two columns
## path and prob, where path has "genotype -> genotype ..."
## and probability is probability
paths_probs_2_df <- function(x, order = c("prob", "path")) {
    order <- match.arg(order)
    df <- data.frame(Path = vapply(x$paths,
                                   function(x) paste(igraph::as_ids(x),
                                                     collapse = " -> "), ""),
                     Prob = exp(x$weights))
    ## Always same ordering
    if (order == "prob")
        oi <- order(-df$Prob, df$Path, decreasing = FALSE)
    else if (order == "path")
        oi <- order(df$Path, decreasing = FALSE)
    return(df[oi, ])
}






## With OncoBN at least, edges of theta 0 lead to
## genotypes that cannot exist. Clean up the adjacency matrix
adjm_rm_no_access <- function(x) {
    ## A single loop should be enough, I think.
    ## But continue until done
    while (TRUE) {
        nacc <- which(colSums(x) == 0)
        wwt <- which(colnames(x) == "WT")
        nacc <- setdiff(nacc, wwt)
        if (length(nacc)) {
            x <- x[-nacc, -nacc]
        } else {
            break
        }
    }
    return(x)
}

## Return parent set from our usual edges structure
parent_set_from_edges <- function(edges) {
    ## Simpler with a table? Oh well
    ps <- aggregate(Relation ~ To,
                    data = edges,
                    FUN = function(x) {
                        if (length(x) == 1) return("Single")
                        else return(x[1])
                    }
                    )
    ps_v <- ps[, "Relation"]
    names(ps_v) <- ps[, "To"]
    return(ps_v)
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
## 1. This is known: should not use OMP on processes under mclapply
## https://luminwin.github.io/randomForestSRC/articles/parallel.html

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
