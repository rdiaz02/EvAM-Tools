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


## No longer using CAPRESE or CAPRI. It gives more output of OT, MHN, CBN than
## other former versions, such as
## code-all-methods-2-trans-matrix-max-genes-15.R. It can use MCCBN if
## MCCBN_INSTALLED is set to TRUE And we can produce plots (CPM-plotting.R)



## You need to have
## - Schill's MHN code
## - MCCBN if you use it.
## - CBN: use my version with fixes, in the repo:
##     cd to ct-cbn-0.1.04b-with-rdu-bug-fix-write-lambda-likelihood
## - OT  (the corresponding package Oncotree)
## - HESBCN
## Full details in the README.md on the upper level directory.




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
###       df_2_access_genots_and_graph* and cpm_access_genots_paths*
###
###
######################################################################
######################################################################

## cpm_access_genots_paths_w* : CPM -> transition matrices
## df_2_access_genots_and_graph*: CPM output -> accessible genotypes

## cpm_access_* calls df_2_access_*

## OT and CBN: df_2_access_genots_and_graph, cpm_access_genots_paths_w_simplified
##    These were the first written. Assume and AND if more than one parent

## DBN: df_2_access_genots_and_graph_OR, cpm_access_genots_paths_w_simplified_OR
##     Derived from previous, assume OR if multiple parents

## HESBCN: df_2_access_genots_and_graph_relationships,
##          cpm_access_genots_paths_w_simplified_relationships
##     Derived from the OT/CBN. Since HESBCN can have AND and XOR and OR, a vector
##     that specifies the relationship if multiple parents is needed.








## Some of this is coming from
## Cancer_Data_sets/CPM-weighted-paths-biol.R
## in the supplementary material for Diaz-Uriarte and Vasallo.
## We leave in there things we don't really need. Simpler.


## DAG of restrictions (as data frame) ->
##              vector of accessible genotypes and graph of DAG of restrictions
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
    ## FIXME: Could be smarter as a function of dim(x)?
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
df_2_access_genots_and_graph_OR <- function(x) {
    
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
##  Called from cpm_access_genots_paths_w_simplified_relationships
##   which itself is used for HESBCN
df_2_access_genots_and_graph_relationships <- function(x,
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
                    all(parents_of_rel %in% genotype)) return(FALSE)
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


## df_2_access_genots_and_graph_OR(x2)
## df_2_access_genots_and_graph_OR(x3)






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

## Based on cpm_access_genots_paths_w  but only with necessary output
##  for both speed and size and using sparse matrices.
##  No paths are computed, so multiple global max is a moot point

## For CBN and MCCBN the weighted fitness graph is the same
## as the transition rate matrix (with 0 in the diagonal)

## We use the same code for OT and CBN, but lambda and OT_edgeWeight are
## fundamentally different things: the first are rates, the second conditional
## probabilities of an untimed oncogenetic tree.  Computing transition
## probabilities between genotypes (as well as "weighted fitness graphs") is
## abusing what OTs strictly provide.
cpm_access_genots_paths_w_simplified <- function(x, 
                                                 names_weights_paths =
                                                     c("rerun_lambda",
                                                       "lambda",
                                                       "OT_edgeWeight")) {
    if(inherits(x, "try-error") || all(is.na(x)) || is.null(x)) {
        return(list(
            fgraph = "ERROR_CPM_ANALYSIS",
            weighted_fgraph = "ERROR_CPM_ANALYSIS",
            trans_mat_genots = "ERROR_CPM_ANALYSIS"
        ))
    }
    
    x <- x$edges
    tmp <- try(df_2_access_genots_and_graph(x[, c("From", "To")]))
    
    if(inherits(tmp, "try-error")) {
        stop("Some error here")
    } else {
        accessible_genots <- tmp$accessible_genots
        fgraph <- unrestricted_fitness_graph_sparseM(accessible_genots)
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
        weighted_fgraph <- NA
    }

    if (length(which_col_weights) == 1) {
        ## weighted_fgraph need not have each row sum to 1.
        ## Obvious if they are lambdas
        ## from CBN for instance.
        ## Neither do they sum to 1 from some OT models.
        ## (see example ex1 in test.trans-rates-f-graphs.R)
        ## So make sure they are transition matrices
        ## between genotypes.
        trans_mat_genots <- rowScaleMatrix(weighted_fgraph)
    } else {
        trans_mat_genots <- NA
    }
    
    return(list(
        fgraph = fgraph,
        ## again: the weighted_fgraph is the transition rate matrix
        ## for CBN and MCCBN. Not for OT.
        weighted_fgraph = weighted_fgraph,
        trans_mat_genots = trans_mat_genots
    ))
}


## output of CPM analysis ->
##             fitness graph (fgraph): just the adj. matrix
##             weighted fitness graph (weighted_fgraph)
##             transition matrix (trans_mat_genots)
## Like cpm_access_genots_paths_w_simplified but for OR relationships

## To be used with output from DBN

## DBN seems similar to OT: these are conditional probabilities, not transition
## rates. Beware of interpretations. See comments in
## cpm_access_genots_paths_w_simplified
cpm_access_genots_paths_w_simplified_OR <-
    function(data,
             parameter_column_name = c("Thetas")) {
##             parameter_column_name = c("Thetas", "Probabilities")) {

        tmp <-
            try(df_2_access_genots_and_graph_OR(data$edges[, c("From", "To")]))

        if (inherits(tmp, "try-error")) {
            stop("Some error here")
        } else {
            accessible_genots <- tmp$accessible_genots
            fgraph <- unrestricted_fitness_graph_sparseM(accessible_genots)
        }

        which_col_weights <-
            which(colnames(data$edges) %in% parameter_column_name)
        
        if (is.null(data$edges[,which_col_weights])){
            stop("No such column")
        }

        weights <- unique(data$edges[, c(2, which_col_weights)])
        rownames(weights) <- weights[, "To"]
        weighted_fgraph <- transition_fg_sparseM(fgraph, weights)
        trans_mat_genots <- rowScaleMatrix(weighted_fgraph)

        return(list(
            fgraph = fgraph,
            ## Seems similar to OT. This ain't the transition rate matrix.
            weighted_fgraph = weighted_fgraph,
            trans_mat_genots = trans_mat_genots
        ))
    }


## output of CPM analysis ->
##             fitness graph (fgraph): just the adj. matrix
##             weighted fitness graph (weighted_fgraph)
##             transition matrix (trans_mat_genots)
## To be used with output from HESBCN

## The weighted fitness graph is the same
## as the transition rate matrix (with 0 in the diagonal)

cpm_access_genots_paths_w_simplified_relationships <-
    function(data,
             parameter_column_name = c("Lambdas")) {             
        ## parameter_column_name = c("Thetas", "Probabilities", "Lambdas")) {

        ## Use the "parent_set" component, which is returned
        ## from HESBCN itself, not the "Relation" component
        ## which we create from the paernt_set. "Relation" is shown
        ## for a human to easily interpret the $edges data frame
        tmp <-
            try(df_2_access_genots_and_graph_relationships(
                data$edges[, c("From", "To")],
                data$parent_set
            ))

        if (inherits(tmp, "try-error")) {
            stop("Some error here")
        } else {
            accessible_genots <- tmp$accessible_genots
            fgraph <- unrestricted_fitness_graph_sparseM(accessible_genots)
        }

        which_col_weights <-
            which(colnames(data$edges) %in% parameter_column_name)
        if (is.null(data$edges[, which_col_weights])) {
            stop("No such column")
        }

        weights <- unique(data$edges[, c(2, which_col_weights)])
        rownames(weights) <- weights[, "To"]
        weighted_fgraph <- transition_fg_sparseM(fgraph, weights)
        trans_mat_genots <- rowScaleMatrix(weighted_fgraph)

        return(list(
            fgraph = fgraph,
            ## weighted_fgraph is the transition rate matrix
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









#' Runs the CPMs
#'
#' Executes all CPMS given a cross sectional data set
#'  
#' @param x cross sectional data
#' @param cores_cbn How many cores to use for CBN
#' @param methods The methods to use. For now, the only thing that matters is
#'     whether the vector includes the strings "MCCBN" and "DBN".
#' @param max.cols Maximum number of columns to use in the analysis. If x has >
#'     max.cols, selected columns are those with the largest number of events.
#' @examples
#'\dontrun{
#' dB_c1 <- matrix(
#'  c(
#'      rep(c(1, 0, 0, 0, 0), 300) #A
#'    , rep(c(0, 0, 1, 0, 0), 300) #C
#'    , rep(c(1, 1, 0, 0, 0), 200) #AB
#'    , rep(c(0, 0, 1, 1, 0), 200) #CD
#'    , rep(c(1, 1, 1, 0, 0), 100) #ABC
#'    , rep(c(1, 0, 1, 1, 0), 100) #ACD
#'    , rep(c(1, 1, 0, 0, 1), 100) #ABE
#'    , rep(c(0, 0, 1, 1, 1), 100) #CDE
#'    , rep(c(1, 1, 1, 0, 1), 100) #ABCE
#'    , rep(c(1, 0, 1, 1, 1), 100) #ACDE
#'    , rep(c(1, 1, 1, 1, 0), 50) # ABCD
#'    , rep(c(0, 0, 0, 0, 0), 10) # WT
#'  ), ncol = 5, byrow = TRUE
#' )
#' colnames(dB_c1) <- LETTERS[1:5]
#' out <- all_methods_2_trans_mat(dB_c1, do_MCCBN = FALSE)
evam <- function(x, cores_cbn = 1,
                 methods = c("CBN", "OT", "HESBCN", "MHN",
                             "DBN"),
                 max.cols = 15) {


    if ("MCCBN" %in% methods) {
        do_MCCBN <- TRUE
    } else {
        do_MCCBN <- FALSE
    }

    if ("DBN" %in% methods) {
        do_DBN <- TRUE
    } else {
        do_DBN <- FALSE
    }

    do_HyperTraPS <- FALSE


    ##########      Preprocessing: common to all methods
    x <- df_2_mat_integer(x)
    xoriginal <- x
    
    x <- add_pseudosamples(x, n00 = "auto3")
    ## remove.constant makes no difference IFF we add pseudosamples, as
    ## there can be no constant column when we add pseudosamples
    x <- pre_process(x, remove.constant = FALSE,
                     min.freq = 0, max.cols = max.cols)

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
        MHN_out <- do_MHN2(x, lambda = 1/nrow(x)))["elapsed"]
    message("time MHN = ", time_MHN)


    ## ####################################################################
    ##    Run each one of the remaining CPMs
    ## ####################################################################
    
    message("Doing HESBCN")
    time_hesbcn <- system.time(
      HESBCN_out <- do_HESBCN(x))["elapsed"]
    message("time HESBCN = ", time_hesbcn)


    
    message("Doing OT")
    time_ot <-
        system.time(
            OT_out <- try(
                suppressMessages(
                    ot_proc(x,
                            nboot = 0,
                            distribution.oncotree = TRUE))))["elapsed"]
    message("time OT = ", time_ot)

    message("Doing CBN")
    time_cbn_ot <- system.time(
           CBN_out <- try(cbn_proc(x,
                                   addname = "tmpo",
                                   init.poset = "OT",
                                   nboot = 0,
                                   parall = TRUE,
                                   cores = 1)))["elapsed"]
    message("time CBN = ", time_cbn_ot)

    if(do_DBN) {
        message("Doing DBN\n\n")
        time_dbn <- system.time(
            DBN_out <- do_DBN(x))["elapsed"]
        message("time DBN = ", time_dbn)
    } 
    
    if(do_MCCBN) {
        if( !requireNamespace("mccbn", quietly = TRUE)) {
            warning("MC-CBN (mccbn) no installed. Not running MC-CBN")
        } else {
            message("Doing MC-CBN")
            time_mccbn <-
                system.time(MCCBN_out <- try(mccbn_proc(x)))["elapsed"]
            message("time MC-CBN = ", time_mccbn)
        }
    }

    ## ######################################################################
    ##   From CPM output, obtain weighted fitness graph and transition matrix
    ##   between genotypes.
    ##
    ##   Beware that HESBCN and DBN use a different function.
    ##   This does not use lapply on purpose.
    ## ######################################################################
    
    OT_fg_tm <- cpm_access_genots_paths_w_simplified(OT_out)
    CBN_fg_tm <- cpm_access_genots_paths_w_simplified(CBN_out)
    if(do_MCCBN)
        MCCBN_fg_tm <- cpm_access_genots_paths_w_simplified(MCCBN_out)
    HESBCN_fg_tm <- cpm_access_genots_paths_w_simplified_relationships(HESBCN_out)
    if(do_DBN)
        DBN_fg_tm <- cpm_access_genots_paths_w_simplified_OR(DBN_out)

    
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
    if (!do_DBN) {
        DBN_out <- list(edges = NA, likelihood = NA)
        DBN_fg_tm <- list(weighted_fgraph = NA, trans_mat_genots = NA)
    }

    HyperTraPS_model <- NA
    HyperTraPS_trans_rate_mat <- NA
    HyperTraPS_trans_mat <- NA
    HyperTraPS_td_trans_mat <- NA

    ## Future: Getting all paths to global maximum
    ## Simply call do_weighted_paths_to_max on a list of
    ## transition matrices and a pre-created paths_to_max
    ## Make sure we are not repeating expensive operations
    ## When testing, compare against cpm_access_genots_paths_w
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

        DBN_model = DBN_out$edges,
        DBN_likelihood = DBN_out$likelihood, 
        DBN_f_graph = DBN_fg_tm$weighted_fgraph, 
        DBN_trans_mat = DBN_fg_tm$trans_mat_genots,

        HESBCN_model = HESBCN_out$edges,
        HESBCN_parent_set = HESBCN_out$parent_set,
        HESBCN_trans_rate_mat = HESBCN_fg_tm$weighted_fgraph,
        HESBCN_trans_mat = HESBCN_fg_tm$trans_mat_genots,
        HESBCN_td_trans_mat = HESBCN_td_trans_mat,

        HyperTraPS_model = HyperTraPS_model,
        HyperTraPS_trans_rate_mat = HyperTraPS_trans_rate_mat,
        HyperTraPS_trans_mat = HyperTraPS_trans_mat,
        HyperTraPS_td_trans_mat = HyperTraPS_td_trans_mat,
        csd_data = xoriginal
         ))
}





## #' Runs the CPMs
## #'
## #' Executes all CPMS given a cross sectional data set
## #'  
## #' @param x cross sectional data
## #' @param cores_cbn How many cores to use for CBN
## #' @param do_MCCBN Wether to run MCCBN. Default is FALSE.
## #' @param max.cols Maximum number of columns to use in the analysis. If x has >
## #'     max.cols, selected columns are those with the largest number of events.
## #' @examples
## #'\dontrun{
## #' dB_c1 <- matrix(
## #'  c(
## #'      rep(c(1, 0, 0, 0, 0), 300) #A
## #'    , rep(c(0, 0, 1, 0, 0), 300) #C
## #'    , rep(c(1, 1, 0, 0, 0), 200) #AB
## #'    , rep(c(0, 0, 1, 1, 0), 200) #CD
## #'    , rep(c(1, 1, 1, 0, 0), 100) #ABC
## #'    , rep(c(1, 0, 1, 1, 0), 100) #ACD
## #'    , rep(c(1, 1, 0, 0, 1), 100) #ABE
## #'    , rep(c(0, 0, 1, 1, 1), 100) #CDE
## #'    , rep(c(1, 1, 1, 0, 1), 100) #ABCE
## #'    , rep(c(1, 0, 1, 1, 1), 100) #ACDE
## #'    , rep(c(1, 1, 1, 1, 0), 50) # ABCD
## #'    , rep(c(0, 0, 0, 0, 0), 10) # WT
## #'  ), ncol = 5, byrow = TRUE
## #' )
## #' colnames(dB_c1) <- LETTERS[1:5]
## #' out <- all_methods_2_trans_mat(dB_c1, do_MCCBN = FALSE)
## all_methods_2_trans_mat <- function(x, cores_cbn = 1,
##                                     methods = c("CBN", "OT", "HESBCN", "MHN",
##                                                 "DBN"),
##                                     max.cols = 15) {


##     if ("MCCBN" %in% methods) {
##         do_MCCBN <- TRUE
##     } else {
##         do_MCCBN <- FALSE
##     }

##     if ("DBN" %in% methods) {
##         do_DBN <- TRUE
##     } else {
##         do_DBN <- FALSE
##     }

##     do_HyperTraPS <- FALSE


##     ## FIXME: the same information is available from some of the
##     ## pre_trans_mat_others and pre_trans_mat_new_CPMS and pre_trans_mat_HESBCN
##     ## Clarify and avoid duplication: costly with large number of genotypes
##     ## Seems to affect DBN and HyperTraPS, and HESBCN 
    
##     ## This function has grown over time and would benefit from rewriting.
##     ## We run MHN, HESBCN, DBN, each in its own call. OT and CBN in the same call.
##     ## Then, as appropriate, we might need to get transition matrices
##     ## on other calls.
##     ## This could be simplified/unified, but inherits from working code.
    
    
##     if(do_MCCBN && !requireNamespace("mccbn", quietly = TRUE))
##         stop("MC-CBN (mccbn) no installed")

##     ## Preprocessing: common to all methods
##     x <- df_2_mat_integer(x)
##     xoriginal <- x
    
##     x <- add_pseudosamples(x, n00 = "auto3")
##     ## remove.constant makes no difference IFF we add pseudosamples, as
##     ## there can be no constant column when we add pseudosamples
##     x <- pre_process(x, remove.constant = FALSE,
##                      min.freq = 0, max.cols = max.cols)

##     ## In the OT-CBN family: which ones we use?
##     if(do_MCCBN)
##         methods <- c("OT",
##                      "CBN_ot"
##                    , "MCCBN")
##     else
##         methods <- c("OT", "CBN_ot")

  
##     cat("\n     Doing MHN")
##     time_schill <- system.time(
##         out_schill <- do_MHN2(x, lambda = 1/nrow(x)))["elapsed"]

##     cat("\n  time MHN = ", time_schill)

##     if(do_DBN) {
##         cat("\n     Doing DBN\n\n")
##         time_dbn <- system.time(
##             out_dbn <- do_DBN(x))["elapsed"]
##         cat("\n  time DBN = ", time_dbn)
##     } else {
##         out_dbn <- NA
##         time_dbn <- NA
##     }

##     if(do_HyperTraPS) {
##     ##     cat("\n     Doing HyperTraps")
##     ##     print("By default we run it here with dry_run = TRUE.
##     ##     HyperTraPS takes a long time and I do no want to block the script.
##         ## ")
##         ## FIXME: HT_folder undefined
##         HT_folder <- NULL
##         time_HyperTraPS <- system.time(
##             out_HyperTraPS <- do_HyperTraPS(x, tmp_folder = HT_folder,
##                                             dry_run = TRUE, plot = FALSE))["elapsed"]
        
##         cat("\n  time HyperTraPS = ", time_HyperTraPS)
##     } else {
##         out_HyperTarPS <- NA
##         time_HyperTraPS <- NA
##     }

    
##     cat("\n     Doing HESBCN")
##     time_hesbcn <- system.time(
##       out_hesbcn <- do_HESBCN(x))["elapsed"]
##     cat("\n  time HESBCN = ", time_hesbcn)

##     cpm_out_others <- ot_cbn_methods(x, cores_cbn = cores_cbn,
##                                      do_MCCBN = do_MCCBN)
##     ## For CBN, OT, get transition matrices and fitness graphs
##     ## Transition matrix only used later for CBN/MCCBN
##     pre_trans_mat_others <- lapply(cpm_out_others[methods],
##         cpm_access_genots_paths_w_simplified)


##     pre_trans_mat_HESBCN <- lapply(
##         list(HESBCN = out_hesbcn
##              ),
##         cpm_access_genots_paths_w_simplified_relationships)
##     pre_trans_mat_others["HESBCN"] <- list(pre_trans_mat_HESBCN$HESBCN)

##     if (do_HyperTraPS) {
##         ## FIXME: this will break, as pre_trans_mat_new_CPMS$HyperTraPS
##         ## has never been created
##         pre_trans_mat_others["HyperTraPS"] <-
##             list(pre_trans_mat_new_CPMS$HyperTraPS)
##     } 
        
##     if (do_DBN) {
##         pre_trans_mat_new_CPMS <- lapply( 
##             list(DBN = out_dbn
##                  ),
##             cpm_access_genots_paths_w_simplified_OR)
##         pre_trans_mat_others["DBN"] <- list(pre_trans_mat_new_CPMS$DBN)
##     } 
        
##     cat("\n    getting transition matrices for all non-mhn methods \n")

##     ## Weighted: the transition matrices between genotypes
##     ## FIXME: why are we calling this on HESBCN? Or on all, for that matter?
##     ## Do not call this at all.
##     wg <-
##         lapply(pre_trans_mat_others[c("OT", "MCCBN" ,
##                                       "CBN_ot", "DBN",
##                                       "HyperTraPS", "HESBCN")[c(TRUE, do_MCCBN,
##                                                                 TRUE, do_DBN,
##                                                                 do_HyperTraPS,
##                                                                 TRUE)]],
##                function(x) x$trans_mat_genots)
##     ## Time discretized, via uniformization
##     ## Only makes sense for methods that return rates
##     td <- lapply(pre_trans_mat_others[c("MCCBN", "CBN_ot",
##                                         "HyperTraPS",
##                                         "HESBCN")[c(do_MCCBN, TRUE,
##                                                     do_HyperTraPS, TRUE)]],
##                  function(x) trans_rate_to_trans_mat(x$weighted_fgraph,
##                                                      method = "uniformization"))

##     ## Return list always has the same elements
##     ## Set to NA those for not run models.
##     if (do_MCCBN) {
##         MCCBN_model <- cpm_out_others$MCCBN$edges
##         MCCBN_f_graph <- pre_trans_mat_others$MCCBN$weighted_fgraph
##         MCCBN_trans_mat <- wg$MCCBN
##         MCCBN_td_trans_mat <- td$MCCBN
##     } else {
##         MCCBN_model <- NA
##         MCCBN_f_graph <- NA
##         MCCBN_trans_mat <- NA
##         MCCBN_td_trans_mat <- NA
##     }
    
##     if (do_DBN) {
##         DBN_model <- out_dbn$edges
##         DBN_likelihood <- out_dbn$likelihood
##         DBN_f_graph <- pre_trans_mat_new_CPMS$DBN$weighted_fgraph
##         DBN_trans_mat <- pre_trans_mat_new_CPMS$DBN$trans_mat_genots
##     } else {
##         DBN_model <- NA
##         DBN_likelihood <- NA
##         DBN_f_graph <- NA
##         DBN_trans_mat <- NA
##     }

##     if (do_HyperTraPS) {
##         HyperTraPS_model <-  out_HyperTraPS$edges
##         HyperTraPS_f_graph <- pre_trans_mat_new_CPMS$HyperTraPS$weighted_fgraph
##         HyperTraPS_trans_mat <- pre_trans_mat_new_CPMS$HyperTraPS$trans_mat_genots
##         HyperTraPS_td_trans_mat <- td$HyperTraPS
##     } else {
##         HyperTraPS_model <- NA
##         HyperTraPS_f_graph <- NA
##         HyperTraPS_trans_mat <- NA
##         HyperTraPS_td_trans_mat <- NA
##     }

##     ## Getting all paths to global maximum
##     ## Simply call do_weighted_paths_to_max on a list of
##     ## transition matrices and a pre-created paths_to_max
##     ## Make sure we are not repeating expensive operations
    
##     return(list(
##         OT_model = cpm_out_others$OT$edges,
##         OT_f_graph = pre_trans_mat_others$OT$weighted_fgraph,
##         OT_trans_mat = wg$OT,
##         OT_genots_predicted = cpm_out_others$OT$genots_predicted,
##         CBN_model = cpm_out_others$CBN_ot$edges,
##         CBN_f_graph = pre_trans_mat_others$CBN_ot$weighted_fgraph,
##         CBN_trans_mat = wg$CBN_ot,
##         CBN_td_trans_mat = td$CBN_ot,
##         MCCBN_model = MCCBN_model,
##         MCCBN_f_graph = MCCBN_f_graph,
##         MCCBN_trans_mat = MCCBN_trans_mat,
##         MCCBN_td_trans_mat = MCCBN_td_trans_mat,
##         MHN_theta = out_schill$theta,
##         MHN_exp_theta = exp(out_schill$theta),
##         MHN_transitionRateMatrix = out_schill$transitionRateMatrix,
##         MHN_trans_mat = out_schill$transitionMatrixCompExp,
##         MHN_td_trans_mat = out_schill$transitionMatrixTimeDiscretized,
##         DBN_model = DBN_model,
##         DBN_likelihood = DBN_likelihood,
##         DBN_f_graph = DBN_f_graph,
##         DBN_trans_mat = DBN_trans_mat,
##         HESBCN_model = out_hesbcn$edges,
##         HESBCN_parent_set = out_hesbcn$parent_set,
##         HESBCN_f_graph = pre_trans_mat_HESBCN$HESBCN$weighted_fgraph,
##         ## FIXME: is this correct? Don't you want wg$HESBCN?
##         HESBCN_trans_mat = pre_trans_mat_HESBCN$HESBCN$trans_mat_genots,
##         HESBCN_td_trans_mat = td$HESBCN,
##         HyperTraPS_model = HyperTraPS_model,
##         HyperTraPS_f_graph = HyperTraPS_f_graph,
##         HyperTraPS_trans_mat = HyperTraPS_trans_mat,
##         HyperTraPS_td_trans_mat = HyperTraPS_td_trans_mat,
##         csd_data = xoriginal
##          ))
## }






#####################################################################

## ## mini example
## Dat1 <- readRDS(file="./MHN/data/BreastCancer.rds") [1:50, 1:4]
## colnames(Dat1) <- LETTERS[1:ncol(Dat1)]

## Dat1 <- as.matrix(Dat1)

## simplerun <- all_methods_2_trans_mat(Dat1)

## rm(simplerun)





## Getting an idea of sizes of data if we did not use sparse matrices
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
