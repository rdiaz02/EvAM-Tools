## Copyright 2021 Ramon Diaz-Uriarte

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


######################################################################
######################################################################
######################################################################

## Use OncoSimulR to obtain accessible genotypes and transition matrices between
## genotypes for CBN, OT, PMCE and DBN.

## See comments and examples in 
## inst/miscell/Using_OncoSimulR_to_get_accessible_genotypes_trans_mats.tex

## Additionally, since we compute and output fitness, can be used for simulations
## (and shows the relationship between lambdas and selection coefficients).

## For CBN, MCCBN, should allow to use time discretization.
##   But this requires proper definition of this process.


## Run as cpm2F2tm(cpm_output)
## See example  below

######################################################################
######################################################################
######################################################################


## fitness, target max fitness. WT fitness always 1.
scale_fitness_2 <- function(x, max_f) {
    max_x <- max(x)
    if(max_x > 1e10) {
        warning("Maximum fitness > 1e10. Expect numerical problems.")
    }
    return(1.0 +  (x - 1) * ((max_f - 1) / (max_x - 1)))
}


## output from CPMs, sh if restrictions not satisfied ->
##                            input for OncoSimulR evalAllGenotypes
cpm_out_to_oncosimul <- function(x, sh = -Inf) {
    sh <- sh
    
    if ("rerun_lambda" %in% names(x)) { ## CBN
        s <- x$rerun_lambda
        typeDep <- "AND"
    } else if ("lambda" %in% names(x)) { ## MCCBN-HCBN
        s <- x$lambda
        typeDep <- "AND"
    } else if ("Relation" %in% names(x)) { ## HESBCN (same thing as PMCE)
        ## Also using this for DBN, as it could return an AND
        if("Lambdas" %in% names(x) ) ## HESBCN
            s <- x$Lambdas
        if("theta" %in% names(x) ) ## DBN
            s <- x$theta
        typeDep <- x$Relation
        typeDep[typeDep == "Single"] <- "AND"
    } else if ("OT_edgeWeight" %in% names(x)) { ## OT
        s <- x$OT_edgeWeight
        typeDep <- "AND"
    } else {
        stop("Input not recognized")
    }
    ## To get the colors right
    ## Fix this in OncoSimulR
    typeDep[typeDep == "AND"] <- "MN"
    typeDep[typeDep == "OR"] <- "SM"
    typeDep[typeDep == "XOR"] <- "XMPN"    
    x1 <- data.frame(parent = x$From,
                     child  = x$To,
                     s = s,
                     sh = sh,
                     typeDep = typeDep
                     )
    return(x1)
}


## output from CPMs, max final fitness, sh when restrictions not stasified,
##         max num genots -> fitness of all genotypes
##    if max_f is NULL: no rescaling of fitness
##       max_f should have no effect in probs transition
##    max_genots: argument max of evalAllGenotypes
cpm_to_fitness_genots <- function(x, max_f = NULL, sh = -Inf, max_genots = 2^15) {
    x1 <- cpm_out_to_oncosimul(x, sh)
    x1 <- evalAllGenotypes(fitnessEffects = allFitnessEffects(rT = x1),
                           addwt = TRUE, max = max_genots)

    ## In newer OncoSimulR, column names for Fitness can now be called Birth
    fitness_birth_column <- ifelse("Fitness" %in% colnames(x1),
                                   "Fitness", "Birth")
    if (!is.null(max_f)) {
        if (max_f < 1) stop("max_f must be larger than min_f")

        if (fitness_birth_column == "Fitness") {
            x1$Fitness[x1$Fitness > 0.0] <-
                scale_fitness_2(x1$Fitness[x1$Fitness > 0.0], max_f)
        } else if (fitness_birth_column == "Birth") {
            x1$Birth[x1$Birth > 0.0] <-
                scale_fitness_2(x1$Birth[x1$Birth > 0.0], max_f)
        } else {
            stop("The column should be called Birth or Fitness")
        }
    }
    return(x1)
}



## named vector of genotype fitness -> fitness graph and transition matrix and
##                                     accessible genotypes
##                                   Only accessible genotypes shown in output.
##   We assume WT is 1 if not given explicitly
genots_2_fgraph_and_trans_mat <- function(x) {
    ## Logic:
    ##  - Construct a fully connected fitness graph
    ##    between the given genotypes. So any genotype connected to
    ##    genotypes with one extra mutation.
    ##    This matrix might contain genotypes that are not truly accessible.
    ##  - Construct matrix of fitness differences between ancestor and immediate
    ##    descendants. Likely slow if many genotypes.
    ##  - Set to non accessible if fitness difference <= 0

    ##  Could be done faster, by not creating the unrestricted fitness graph
    ##  and instead maybe using OncoSimulR's wrap_accessibleGenotypes
    ##  But would need to check that works with partial lists of genotypes
    ##  and use allGenotypes_to_matrix. And would need to change
    ##  OncoSimulR's wrap_accessibleGenotype and how it uses the th.

    ## THINK: but the call to genots_2_fgraph_and_trans_mat inside
    ## cpm_to_trans_mat_oncosimul has already used wrap_accessibleGenotypes
    ## and only accessible genotypes are passed to this function.
    ## Note, though, that this function has been extenssively tested.

    ## We use this approach, to minimize the number of genotypes we call
    ## unrestricted_fitness_graph_sparseM in cpm_to_trans_mat_oncosimul

    which_wt <- which(names(x) == "WT")
    if (length(which_wt) == 1) {
        fit_wt <- x["WT"]
        if (fit_wt != 1.0) message("Your WT has fitness different from 1.",
                                  " Using WT with the fitness you provided.")
        x <- x[-which_wt]
    } else {
        fit_wt <- c("WT" = 1.0)
    }

    ## Silly? Inside the next, we now put them together. FIXME?
    access_genots_as_list <- lapply(names(x),
                                    function(v) strsplit(v, ", ")[[1]])
    ## For ordered output
    no <-  order(unlist(lapply(access_genots_as_list, length)), names(x))
    access_genots_as_list <- access_genots_as_list[no]

    fgraph <- unrestricted_fitness_graph_sparseM(access_genots_as_list)

    genots_fitness <- c(fit_wt, x)[colnames(fgraph)]
    mf <- matrix(rep(genots_fitness, nrow(fgraph)),
                 nrow = nrow(fgraph), byrow = TRUE)
    stopifnot(identical(dim(mf), dim(fgraph)))

    fdiff <- mf - genots_fitness
    fdiff <- fgraph * fdiff
    fgraph[fdiff <= 0] <- 0
    tm <- fdiff
    tm[tm < 0] <- 0
    tm <- tm / ifelse(rowSums(tm) != 0, rowSums(tm), 1)
    ## This we know is always 0
    tm[nrow(tm), ] <- 0
    ## This we can set
    tm[rowSums(fgraph) == 0, ] <- 0

    ## First simple filtering for potentially expensive call to igraph
    accessible_genotypes_candidates <-
        genots_fitness[colnames(fgraph)[colSums(fgraph) >= 1]]

    ig_fgraph <- igraph::graph_from_adjacency_matrix(fgraph)

    num_paths_from_WT <-
        vapply(names(accessible_genotypes_candidates),
               function(x)
                   length(igraph::all_simple_paths(ig_fgraph,
                                                   from = "WT",
                                                   to = x, mode = "out")),
               FUN.VALUE = 0)

    accessible_genotypes <-
        accessible_genotypes_candidates[num_paths_from_WT > 0]

    ## Remove unaccessibe genotypes from matrices before returning
    name_ret <- c("WT", names(accessible_genotypes))

    return(
        list(fitness_graph = fgraph[name_ret, name_ret],
             transition_matrix = tm[name_ret, name_ret],
             fitness_differences = fdiff[name_ret, name_ret],
             accessible_genotypes = accessible_genotypes))

}





## output from CPMs,  max final fitness, sh if restrictions not satisfied,
##               max num genots ->
##                            list(accessible genotypes ,
##                                 fitness graph,
##                                 transition matrix between genotypes)
##    max_f  should have no effect in probs transition
##    accessible_genotypes: gives the genotypes and their fitness,
##                          computed from the lambdas (thetas) and possibly scaled
##   By default, WT has fitness 1.
cpm_to_trans_mat_oncosimul <- function(x, max_f = NULL, sh = -Inf,
                                       max_genots = 2^15,
                                       WT_fitness = 1.0) {
    fitness_all <- cpm_to_fitness_genots(x, max_f = max_f, sh = sh,
                                         max_genots = max_genots)

    ## If we were to use wrap_accessibleGenotypes we would not
    ## do this filtering, so that we get all the genotypes.
    ## O.w. we can filter because from OncoSimulR, WT always has
    ## fitness 1. WT_fitness = 1.

    ## access_genots_fitness <-
    ##     fitness_all$Fitness[fitness_all$Fitness > WT_fitness]
    ## names(access_genots_fitness) <-
    ##     fitness_all$Genotype[fitness_all$Fitness > WT_fitness]

    
    ## Using wrap_accessibleGenotypes likely to be faster
    ## and this will only give truly accessible from WT, except
    ## for the borderline cases of equal fitness ancestor-descendant
    
    fitness_mat <- evam_allGenotypes_to_matrix(fitness_all)
    access_genots <- evam_wrap_accessibleGenotypes(fitness_mat, th = 0)

    ## In newer OncoSimulR, column names for Fitness can now be called Birth

    fitness_birth_column <- ifelse("Fitness" %in% colnames(fitness_all),
                                   "Fitness", "Birth") 
    
    access_genots_fitness <- fitness_all[access_genots, fitness_birth_column]
    names(access_genots_fitness) <- fitness_all[access_genots, "Genotype"]

    
    gfo <- genots_2_fgraph_and_trans_mat(access_genots_fitness)

    return(list(accessible_genotypes = gfo$accessible_genotypes,
                fitness_graph = gfo$fitness_graph,
                ## For checking; these are not necessarily the original lambdas
                ## but the original multiplied by a factor, like changing the rate of all.
                ## lambdas = gfo$fitness_differences, 
                transition_matrix = gfo$transition_matrix
                ))
}

## shorter
cpm2F2tm <- cpm_to_trans_mat_oncosimul

