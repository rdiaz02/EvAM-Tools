## Copyright 2021 Ramon Diaz-Uriarte

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


######################################################################
######################################################################
######################################################################

## Use OncoSimulR to obtain accessible genotypes and transition matrices between
## genotypes for CBN, OT, PMCE and DBN.

## See comments and examples in ../OncoSimul-for-accessible-and-more-on-PMCE.org 


## Additionally, since we compute and output fitness, can be used for simulations
## (and shows the equivalence between lambdas and selection coefficients).

## Code for CBN and OT tested. For DBN not tested (nor fully implemented). For
## PMCE not tested.

## FIXME:
## For CBN, MCCBN, should allow to use time discretization.
##   But this requires proper definition of this process.


## Run as cpm2tm(cpm_output)
## See example  below

######################################################################
######################################################################
######################################################################



## TODO lots of things could be removed from this file

## fitness, target max fitness. WT fitness always 1.
scale_fitness_2 <- function(x, max_f) {
    max_x <- max(x)
    return(1.0 +  (x - 1) * ((max_f - 1) / (max_x - 1)))
}


## output from CPMs, sh if restrictions not satisfied -> input for OncoSimulR evalAllGenotypes

cpm_out_to_oncosimul <- function(x, sh = -Inf) {
    sh <- sh
    
    if("rerun_lambda" %in% names(x)) { ## CBN
        s <- x$rerun_lambda
        typeDep <- "AND"
    } else if("Relation" %in% names(x)) { ## PMCE
        if(exists("Lambdas", where = x))
            s <- x$Lambdas
        # else if(exists("lambda", where = x))
        #     s <- x$lambda
        typeDep <- x$Relation
        typeDep[typeDep == "Single"] <- "AND"
    } else if("OT_edgeWeight" %in% names(x) ) { ## OT
        s <- x$OT_edgeWeight
        typeDep <- "AND"
    } else if("Thetas" %in% names(x)) { ## Something for DB
        s <- x$Thetas
        typeDep <- "OR"
        
    } else if("otro" %in% names(x) ) { ## MCCBN
        
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
    if(!is.null(max_f)) {
        if(max_f < 1) stop("max_f must be larger than min_f")
        x1$Fitness[x1$Fitness > 0.0] <-
            scale_fitness_2(x1$Fitness[x1$Fitness > 0.0], max_f)
    }
    return(x1)
}


## named vector of genotype fitness -> fitness graph and transition matrix and
##                                     accessible genotypes
##   We assume WT is 1 if not given explicitly
##   BEWARE: only tested with one minimal test, and by comparing against CBN and OT!! 
genots_2_fgraph_and_trans_mat <- function(x) {
    ## Logic:
    ##  - Construct a fully connected fitness graph
    ##    between the given genotypes. So any genotype connected to
    ##    genotypes with one extra mutation.
    ##  - Construct matrix of fitness differences between ancestor and immediate
    ##    descendants. Likely slow if many genotypes.
    ##  - Set to non accessible if fitness difference <= 0

    ##  Could be done faster, by not creating the unrestricted fitness graph
    ##  and instead maybe using OncoSimulR's wrap_accessibleGenotypes
    ##  But would need to check that works with partial lists of genotypes
    ##  and use allGenotypes_to_matrix.

    ## FIXME: use igraph's function to all genots as destinations?
    ## NOPE! use wrap_accessibleGenotypes for accessible
    ##      what for fitness graph?
    ##      and this code for transition
    ##      And then, check for consistency.

    
    which_wt <- which(names(x) == "WT")
    if(length(which_wt) == 1) {
        fit_wt <- x["WT"]
        if(fit_wt != 1.0) message("Your WT has fitness different from 1.",
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
    accessible_genotypes <- genots_fitness[colnames(fgraph)[colSums(fgraph) >= 1]]

    return(list(fitness_graph = fgraph,
           transition_matrix = tm,
           fitness_differences = fdiff,
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
cpm_to_trans_mat_oncosimul <- function(x, max_f = NULL, sh = -Inf,
                                       max_genots = 2^15) {
    fitness_all <- cpm_to_fitness_genots(x, max_f = max_f, sh = sh,
                                         max_genots = max_genots)

    access_genots_fitness <- fitness_all$Fitness[fitness_all$Fitness >= 1.0]
    names(access_genots_fitness) <- fitness_all$Genotype[fitness_all$Fitness >= 1.0]

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
cpm2tm <- cpm_to_trans_mat_oncosimul

# library(codetools)
# checkUsageEnv(env = .GlobalEnv)


# source(file = "tests-access_genots_from_oncosimul.R", echo = TRUE)

######################################################################
######################################################################

#########      Examples

###### CBN
# ex_cbn_out2 <- structure(list(From = c("Root", "A", "Root", "C"),
#                               To = c("A", "B", "C", "D"),
#                               edge = c("Root -> A", "A -> B", "Root -> C", "C -> D"),
#                               init_lambda = c(0.689845, 2.734304, 0.6988, 2.463583),
#                               final_lambda = c(0.689845, 2.734304, 0.6988, 2.463583),
#                               rerun_lambda = c(0.689845, 2.734304, 0.6988, 2.463583),
#                               CBN_edgeBootFreq = c(NA, NA, NA, NA)),
#                          class = "data.frame", row.names = c("A", "B", "C", "D"))

# cpm2tm(ex_cbn_out2)

# ######    PMCE
# ## From email 20-May-2021, at 17:32:51
# ex_pmce_out1 <- read.table("ex_pmce_out1.txt", header = TRUE)

# outp1 <- cpm2tm(ex_pmce_out1)

# ## plot

# plot(allFitnessEffects(cpm_out_to_oncosimul(ex_pmce_out1)))

# plot(allFitnessEffects(cpm_out_to_oncosimul(ex_pmce_out1)),
#      "igraph", layout = layout.reingold.tilford)

## 
# ex_pmce_email <- read.table("ex_pmce_email.txt", header = TRUE)

# out_em <- cpm2tm(ex_pmce_email)

# pdf(file = "DAG_pmce_example.pdf")
# plot(allFitnessEffects(cpm_out_to_oncosimul(ex_pmce_email)))
# dev.off()

## 
# ex_pmce_email <- read.table("ex_pmce_email.txt", header = TRUE)

# out_em <- cpm2tm(ex_pmce_email)

# pdf(file = "DAG_pmce_example.pdf")
# plot(allFitnessEffects(cpm_out_to_oncosimul(ex_pmce_email)))
# dev.off()



# #### Tiny OR and XOR

# ex_or <- read.table("ex_pmce_or.txt", header = TRUE, sep = "\t")
# out_or <- cpm2tm(ex_or)


# ex_xor <- read.table("ex_pmce_xor.txt", header = TRUE, sep = "\t")
# out_xor <- cpm2tm(ex_xor)


# ## Load the stomach output
# load("stomach_pmce.RData")

# stomach_out <- cpm2tm(stomach_pmce)

# plot(allFitnessEffects(
#     cpm_out_to_oncosimul(stomach_pmce)),
#      "igraph", layout = layout.reingold.tilford)


# plot(allFitnessEffects(
#     cpm_out_to_oncosimul(stomach_pmce)))
     


# stomach_out$accessible_genotypes


