library(parallel)

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

## theta10 <- Random.Theta(n=10, sparsity=0.50)

## save(file = "three_thetas.RData", theta8, theta4, theta10)


source("schill-trans-mat.R")

colnames(theta4) <- rownames(theta4) <- LETTERS[1:4]
colnames(theta8) <- rownames(theta8) <- LETTERS[1:8]
colnames(theta10) <- rownames(theta10) <- LETTERS[1:10]

trm4 <- theta_to_trans_rate_3_SM(theta4)
trm8 <- theta_to_trans_rate_3_SM(theta8)
trm10 <- theta_to_trans_rate_3_SM(theta10)

save(file = "three_trm.Data", trm4, trm8, trm10)

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





######################################

## We repeat the sums to compute the diagonal and the division
## If we sample a large number of times, possibly worth it to
## have those precomputed

## This is what I call "transition matrix standardized":
##    Diagonal is passed separately, entries in matrix are probabilities

## This will only be used in called after standardization
## therefore the "is.null" can be removed

## transition rate matrix "standardized",
##  diagonal of transition rate matrix, time of sampling of a case/individual ->
##                       sampled genotype, trajectory, and accumulated time
indiv_sample_from_trm_pre <- function(trmstd,
                                      diag,
                                      T_sampling,
                                      ngenots,
                                      genot_names) {
    ## if(is.null(ngenots)) ngenots <- ncol(trm)
    ## if(is.null(genot_names)) genot_names <- colnames(trm)
    row <- 1
    t_accum <- 0
    genotype <- "WT" 
    trajectory <- "WT"
    
    while(TRUE) {
        ## qii <- diag[row] ## Or qs in Gotovos et al terminology
        ## Special case of genotype that does not transition to anything
        if(diag[row] == 0.0) break 
        t_transition <- rexp(n = 1, rate = diag[row])
        t_accum <- t_accum + t_transition
        if(t_accum >= T_sampling ) break
        ## We transition. To what genotype?
        ## pj <- trm[row, ]
        j <- sample(x = seq_len(ngenots), size = 1, prob = trmstd[row, ])
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

## transition rate matrix, number of samples or times of samples,
##  whether or not to precompute entries of the trans rate matrix (for speed)
##  number of cores (pass 1 is you do not want to parallelize)

population_sample_from_trm <- function(trm, n_samples = 10,
                                       T_sampling = NULL,
                                       pre_compute = TRUE,
                                       cores = detectCores()) {
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

    if(pre_compute) {
        ## Like code in trans_rate_to_trans_mat
        sx <- rowSums(trm)
        ii <- which(sx > 0)
        for(i in ii) {
            trm[i, ] <- trm[i, ]/sx[i]
        }
        
        out <- mclapply(T_sampling,
               function(x)
                   indiv_sample_from_trm_pre(trmstd = trm,
                                             diag = sx,
                                             T_sampling = x,
                                             ngenots = ngenots,
                                             genot_names = genot_names),
               mc.cores = cores)   
        
    } else {    
        out <- mclapply(T_sampling,
               function(x)
                   indiv_sample_from_trm(trm = trm,
                                         T_sampling = x,
                                         ngenots = ngenots,
                                         genot_names = genot_names),
               mc.cores = cores)
    }
    ## Structure output as Pablo's  simulate_population_2
    ## Otherwise, we could just exist from the above
    ## This will add time and increase RAM usage

    return(list(
        T_sampling = T_sampling
      , T_sum_events = unlist(lapply(out, function(x) x$t_accum))
      , trans_table = NA ## I do not know what this is for
      , trajectory = lapply(out, function(x) x$trajectory)
      , obs_events = unlist(lapply(out, function(x) x$genotype))
        ))
}



## Examples

population_sample_from_trm(trm4, 12)

population_sample_from_trm(trm8, T_sampling = rexp(23))

population_sample_from_trm(trm8, T_sampling = rexp(10000))

uu <- population_sample_from_trm(trm10, T_sampling = rexp(10000))


## 2.6 secs
system.time(null <- population_sample_from_trm(trm8, T_sampling = rexp(10000)) )

## 46 secs
system.time(null <- population_sample_from_trm(trm8, T_sampling = rexp(200000)) )


## 3.3 secs
system.time(null <- population_sample_from_trm(trm10, T_sampling = rexp(10000)) )


## 47 secs
system.time(null <- population_sample_from_trm(trm10, T_sampling = rexp(200000)) )



## Take a sample (a vector), with genotypes as "A, B", etc
## and return a vector of frequencies (counts) in the exact same
## order as used by MHN
## A simple implementation that can be slow when the sample is large

## vector of genotypes, total number of genes ->
##             counts of all genotypes in same order as used by MHN
sample_to_pD_order0 <- function(x, ngenes) {
    x <- gsub("WT", "", x)
    xs <- strsplit(x, ", ")

    Data <-  do.call(rbind,
                   lapply(xs, function(z) as.integer(LETTERS[1:ngenes] %in% z)))

    ## What follows is from Data.to.pD
    ## except we do not divide
    n <- ncol(Data)
    N <- 2^n
    Data <- apply(Data, 1, State.to.Int)
    pD <- tabulate(Data, nbins = N)
    return(pD)
}

## Take a sample (a vector), with genotypes as "A, B", etc
## and return a vector of frequencies (counts) in the exact same
## order as used by MHN
## A much faster implementation

## vector of genotypes, total number of genes ->
##             counts of all genotypes in same order as used by MHN
sample_to_pD_order <- function(x, ngenes) {
    x <- as.data.frame(table(x), stringsAsFactors = FALSE)
    
    genot_int <- x[, 1]
    genot_int <- gsub("WT", "", genot_int)
    genot_int <- vapply(genot_int,
                        function(z)
                            State.to.Int(as.integer(LETTERS[1:ngenes] %in%
                                                    strsplit(z, ", ")[[1]])),
                        numeric(1))
    ## all_genots <- rep(unname(genot_int), x[, 2])
    pD <- tabulate(rep(unname(genot_int), x[, 2]),
                   nbins = 2^ngenes)
}


## xa <- sample_to_pD_order(uu$obs_events, 4)
## xb <- sample_to_pD_order0(uu$obs_events, 4)
## xc <- sample_to_pD_order(null$obs_events, 8)
