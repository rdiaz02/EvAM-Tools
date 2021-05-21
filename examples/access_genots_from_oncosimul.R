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

## Additionally, since we compute and output fitness, can be used for simulations
## (and shows the equivalence between lambdas and selection coefficients).

## Code for CBN and OT tested. For DBN not tested (nor fully implemented). For
## PMCE not tested.

## For CBN, MCCBN, should allow to use time discretization. Not done yet. FIXME?


## Run as cpm2tm(cpm_output)
## See example  below

######################################################################
######################################################################
######################################################################


## Required deps
## Yes, slow. Will need to separate testing.
pwd0 <- getwd()
setwd("../code_from_what_genotype_next/")
source("code-all-methods-minimal.R")
setwd(pwd0)
rm(pwd0)


## fitness, target max fitness. WT fitness always 1.
scale_fitness_2 <- function(x, max_f) {
    max_x <- max(x)
    return(1.0 +  (x - 1) * ( (max_f - 1)/(max_x - 1) ))
}


## output from CPMs, fitness scale -> input for OncoSimulR evalAllGenotypes
##   min_f: must pass this argument: minimum fitness
cpm_out_to_oncosimul <- function(x) {
    sh <- -Inf
    
    if("rerun_lambda" %in% names(x)) { ## CBN
        s <- exp(x$rerun_lambda) - 1
        typeDep <- "AND"
    } else if("Relation" %in% names(x)) { ## PMCE
        s <- exp(x$Lambda) - 1        
        typeDep <- x$Relation
        typeDep[typeDep == "Single"] <- "AND"
    } else if("OT_edgeWeight" %in% names(x) ) { ## OT
        s <- exp(x$OT_edgeWeight) - 1                
        typeDep <- "AND"
    } else if("whatever" %in% names(x)) { ## Something for DB
        
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


## output from CPMs, min and max final fitness, max num genots -> fitness of all genotypes
##    if max_f is NULL: no rescaling of fitness
##       max_f and min_f should have no effect in probs transition
##    max_genots: argument max of evalAllGenotypes
cpm_to_fitness_genots <- function(x, max_f = 3, max_genots = 2^15) {
    x1 <- cpm_out_to_oncosimul(x)
    x1 <- evalAllGenotypes(fitnessEffects = allFitnessEffects(rT = x1),
                           addwt = TRUE, max = max_genots)
    x1$Fitness[x1$Fitness > 0.0] <- log(x1$Fitness[x1$Fitness > 0.0]) + 1
    if(!is.null(max_f)) {
        if(max_f < 1) stop("max_f must be larger than min_f")
        x1$Fitness[x1$Fitness > 0.0] <-
            scale_fitness_2(x1$Fitness[x1$Fitness > 0.0], max_f)
    }
    return(x1)
}


## output from CPMs, min and max final fitness, max num genots ->
##                            list(accessible genotypes ,
##                                 fitness graph,
##                                 transition matrix between genotypes)
##    max_f and min_f should have no effect in probs transition
##    accessible_genotypes: gives the genotypes and their fitness,
##                          computed from the lambdas (thetas) and possibly scaled
cpm_to_trans_mat_oncosimul <- function(x, max_f = 3,
                                       max_genots = 2^15) {
    fitness_all <- cpm_to_fitness_genots(x, max_f = max_f,
                                         max_genots = max_genots)

    access_genots_fitness <- fitness_all$Fitness[fitness_all$Fitness > 1.0]
    names(access_genots_fitness) <- fitness_all$Genotype[fitness_all$Fitness > 1.0]
    
    ## Silly? Inside the next, we now put them together. FIXME?
    access_genots_as_list <- lapply(names(access_genots_fitness),
                                    function(v) strsplit(v, ", ")[[1]])

    fgraph <- unrestricted_fitness_graph_sparseM(access_genots_as_list)
    
    ## Not very efficient, but works
    o_access_genots_fitness <- c("WT" = 1, access_genots_fitness)[colnames(fgraph)]
    mf <- matrix(rep(o_access_genots_fitness, nrow(fgraph)),
                 nrow = nrow(fgraph), byrow = TRUE)
    ## As difference w.r.t. original genotype
    mf <- mf - o_access_genots_fitness
    tm <- fgraph * mf
    tm_s <- tm / rowSums(tm)
    tm_s[nrow(tm_s), ] <- 0
    return(list(accessible_genotypes = access_genots_fitness,
                fitness_graph = fgraph,
                ## For checking; these are not necessarily the original lambdas
                ## but a the original multiplied by a factor, like changing the rate of all.
                ## lambdas = tm, 
                transition_matrix = tm_s
                ))
}

## shorter
cpm2tm <- cpm_to_trans_mat_oncosimul

library(codetools)
checkUsageEnv(env = .GlobalEnv)


source(file = "tests-access_genots_from_oncosimul.R", echo = TRUE)

######################################################################
######################################################################

#########      Examples

###### CBN
ex_cbn_out2 <- structure(list(From = c("Root", "A", "Root", "C"),
                              To = c("A", "B", "C", "D"),
                              edge = c("Root -> A", "A -> B", "Root -> C", "C -> D"),
                              init_lambda = c(0.689845, 2.734304, 0.6988, 2.463583),
                              final_lambda = c(0.689845, 2.734304, 0.6988, 2.463583),
                              rerun_lambda = c(0.689845, 2.734304, 0.6988, 2.463583),
                              CBN_edgeBootFreq = c(NA, NA, NA, NA)),
                         class = "data.frame", row.names = c("A", "B", "C", "D"))

cpm2tm(ex_cbn_out2)

######    PMCE
## From email 20-May-2021, at 17:32:51
ex_pmce_out1 <- read.table("ex_pmce_out1.txt", header = TRUE)

outp1 <- cpm2tm(ex_pmce_out1)

## plot

plot(allFitnessEffects(cpm_out_to_oncosimul(ex_pmce_out1)))

plot(allFitnessEffects(cpm_out_to_oncosimul(ex_pmce_out1)),
     "igraph", layout = layout.reingold.tilford)



#### Tiny OR and XOR

ex_or <- read.table("ex_pmce_or.txt", header = TRUE, sep = "\t")
out_or <- cpm2tm(ex_or)


ex_xor <- read.table("ex_pmce_xor.txt", header = TRUE, sep = "\t")
out_xor <- cpm2tm(ex_xor)









