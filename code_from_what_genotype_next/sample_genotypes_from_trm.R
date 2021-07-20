

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


## Note differences
##  I do not compute diagonals unless needed

##  I pass the trm directly, not objects from the trm that will use a lot of RAM
##  and possibly CPU on creation


## transition rate matrix, time of sampling of a case/individual ->
##                       sampled genotype, trajectory, and accumulated time
trm_2_sample <- function(trm, T_sampling) {
    ngenots <- ncol(trm)
    genot_names <- colnames(trm)
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


