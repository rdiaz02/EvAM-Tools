

## First, I create a few transition rate matrices from MHN

## setwd("./MHN")
## source("UtilityFunctions.R")
## source("ModelConstruction.R")
## source("Likelihood.R")
## source("RegularizedOptimization.R")

## set.seed(1)

## theta8 <- Random.Theta(n=8, sparsity=0.50)
## pTh8 <- Generate.pTh(theta8)

## set.seed(1)
## theta4 <- Random.Theta(n=4, sparsity=0.50)
## pTh4 <- Generate.pTh(theta4)

## save(file = "two_thetas.RData", theta8, theta4)


source("schill-trans-mat.R")

colnames(theta4) <- rownames(theta4) <- LETTERS[1:4]
colnames(theta8) <- rownames(theta8) <- LETTERS[1:8]

trm4 <- theta_to_trans_rate_3_SM(theta4)
trm8 <- theta_to_trans_rate_3_SM(theta8)

## indiv_sample_from_trm is equivalent to simulate_sample_2
## population_sample_from_trm is equivalent to simulate_population_2



## Note differences with "simulate_sample_2"
##  I do not compute diagonals unless needed

##  I pass the trm directly, not objects from the trm that will use a lot of RAM
##  and possibly CPU on creation


## transition rate matrix, time of sampling of a case/individual ->
##                       sampled genotype, trajectory, and accumulated time
indiv_sample_from_trm <- function(trm, T_sampling, ngenots = NULL,
                            genot_names = NULL) {
    if(is.null(ngenots)) ngenots <- ncol(trm)
    if(is.null(genot_names)) genot_names <- colnames(trm)
    row <- 1
    t_accum <- 0
    genotype <- "WT" 
    trajectory <- "WT"
    
    while(TRUE) {
        qii <- sum(trm[row, ]) ## Or qs in Gotovos et al terminology
        ## Special case of genotype that does not transition to anything
        if(qii == 0.0) break 
        t_transition <- rexp(n = 1, rate = qii)
        t_accum <- t_accum + t_transition
        if(t_accum >= T_sampling ) break
        ## We transition. To what genotype?
        pj <- trm[row, ]/qii
        j <- sample(x = seq_len(ngenots), size = 1, prob = pj)
        genotype <- genot_names[j]
        trajectory <- c(trajectory, genotype)
        ## For next iteration
        row <- j
    }
    return(list(genotype = genotype,
                trajectory = trajectory,
                t_accum = t_accum))
}


## Like indiv_sample_from_trm, but for multiple times
population_sample_from_trm <- function(trm, n_samples = 10, T_sampling = NULL) {
    if(is.null(T_sampling) && !is.null(n_samples)) {
        T_sampling <- rexp(n = n_samples, rate = 1)
    }
    if(!is.null(T_sampling) && !is.null(n_samples)) {
        message("Ignoring n_samples as passing T_sampling")
    }
    if(is.null(T_sampling) && is.null(n_samples)) {
        stop("Pass either n_samples or T_sampling vector")
    }

    
    ngenots <- ncol(trm)
    genot_names <- colnames(trm)
    lapply(T_sampling,
           function(x) indiv_sample_from_trm(trm = trm,
                                             T_sampling = x,
                                             ngenots = ngenots,
                                             genot_names = genot_names))
}



## Examples

population_sample_from_trm(trm4, 12)

population_sample_from_trm(trm8, T_sampling = rexp(23))
