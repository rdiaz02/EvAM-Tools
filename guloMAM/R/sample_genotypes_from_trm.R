
#' @title Sample an indivial based on a transition rate matrix
#' 
#' @param trm transition rate matrix
#' @param T_sampling Time to compute 
#' @param ngenots Number of genotypes
#' @param genot_names String array with genotype names
#' 
#' @return sampled genotype, trajectory, and accumulated time
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


#' @title Sample an indivial based on a transition rate matrix
#' 
#' We repeat the sums to compute the diagonal and the division
#' If we sample a large number of times, possibly worth it to
#' have those precomputed

#' This is what I call "transition matrix standardized":
#'    Diagonal is passed separately, entries in matrix are probabilities

#' This will only be used in called after standardization
#' therefore the "is.null" can be removed

#' @param trmstd transition rate matrix "standardized",
#' @param diag diagonal of transition rate matrix, time of sampling of a case/individual
#' @param T_sampling Time to compute 
#' @param ngenots Number of genotypes
#' @param genot_names String array with genotype names
#' 
#' @return sampled genotype, trajectory, and accumulated time
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


#' @title Sample a population
#' 
#' @description Like indiv_sample_from_trm, but for multiple times
#' 
#' @param trm transition rate matrix, number of samples or times of samples,
#' @param n_samples Int with the number of samples to be computed
#' @param T_sampling Time at wich each individual in sample. By default they 
#' are randomly generated
#' @param pre_compute whether or not to precompute entries of the trans rate matrix (for speed)
#' @param cores number of cores (pass 1 is you do not want to parallelize)
#' 
#' @return List with precompunted sampling time of sampling, the actual time sampled
#' observed for each sample (s), the complete trajectory of acquired mutations
#' and the observed genotype

population_sample_from_trm <- function(trm, n_samples = 10,
                                       T_sampling = NULL,
                                       pre_compute = TRUE,
                                       cores = detectCores()) {
    if(is.null(T_sampling) && is.null(n_samples)) {
        stop("Pass either n_samples or T_sampling vector")
    }
    if(!is.null(T_sampling) && !is.null(n_samples)) {
        message("Ignoring n_samples as passing T_sampling")
    }
    if(is.null(T_sampling) && !is.null(n_samples)) {
        T_sampling <- rexp(n = n_samples, rate = 1)
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
    #   , trans_table = NA ## I do not know what this is for
      , trajectory = lapply(out, function(x) x$trajectory)
      , obs_events = unlist(lapply(out, function(x) x$genotype))
        ))
}

#' @title Process samples
#' 
#' @description Generate trajectories from simulated data
#' 
#' @param sim list generated with mccbn::sample_genotypes. Relevant
#' fields are described below
#' $T_sum_events time of events for the mutations of each gene
#' $obs_events data.frame with mutated before the end of the sampling time
#' @param n_genes number of genes observed
#' @param output type of output that we want
#' 
#' @return List with a list of trajectories (the order in which gene mutations
#' are acquired), genotype frequencies and genotypes transition matrix (with
#' counts of how many transitions between each genotype have been observed) 
process_samples <- function(sim, n_genes, output = c("frequencies", "state_counts", "transitions")){

    #Checking input
    params <- c("trajectory", "obs_events")
    for (i in params){
        if (!(i %in% names(sim))) 
            stop(sprintf("%s is missing from your samples", i))
    }

    #Checking output variables
    valid_output <- c("frequencies", "state_counts", "transitions")
    out_params <- valid_output %in% output
    names(out_params) <- valid_output

    if (sum(out_params) == 0) stop("Specify valid output")

    not_valid_params <- output[which(!(output %in% valid_output))]
    if (length(not_valid_params) > 0) 
        warning(sprintf("The following parameters cannot be returned: %s"
            , paste(not_valid_params, collapse = ", " )))

    #Set up
    output <- list()
    n_states <- 2**n_genes
    sorted_genotypes <- vapply(0:(n_states - 1), int2str, character(1))
    trajectories <- sim$trajectory

    #Calculate frequencies
    if(out_params["frequencies"]){
        frequencies <- sample_to_pD_order(sim$obs_events, n_genes)
        frequencies <- data.frame(
            Genotype = sorted_genotypes,
            Counts = frequencies
        )
        rownames(frequencies) <- NULL

        output$frequencies <- frequencies
    }

    #Calculate transitions
    if(out_params["transitions"]){
        t <- matrix(0L, nrow = n_states, ncol = n_states)
        colnames(t) <- rownames(t) <- sorted_genotypes
        for(traj in trajectories){
            steps <- length(traj) - 1 
            if(steps > 0){
                for(i in 1:steps){
                    t[traj[i], traj[i + 1]] <-
                        t[traj[i], traj[i + 1]] + 1
                }
            }
        } 

        output$transitions <- t
    }

    #Calculate state_counts
    if(out_params["state_counts"]){
        state_counts <- sample_to_pD_order(unlist(sim$trajectory), n_genes)
        state_counts <- data.frame(
            Genotype = sorted_genotypes,
            Counts = state_counts
        )

        output$state_counts <- state_counts
    }

    return(output)
}

#' @title Run samples for all outputs of CPMs
#' 
#' TODO some methods do not make sense to be here, like OT
#' that do not raise a transition rate matrix
#' 
#' @param cpm_output Output from calling all_methods2trans_mat
#' @param N_samples Number of samples to generate
#' @param n_genes Number of samples that are in the sample
#' @param methods List of methods that we want to sample
#' 
#' @return modified cpm_outputd including a matrix with genotype transitions
sample_all_CPMS <- function(cpm_output
    , N_samples
    , n_genes
    , methods = c("OT", "CBN", "MCCBN", "DBN", "MHN", "HESBCN")){
    output <- cpm_output
    

    for(method in methods){
        if(method == "MHN") trm <- output$MHN_transitionRateMatrix
        else trm <- output[[sprintf("%s_f_graph", method)]]

        if(any(!is.na(trm))){
            sims <- population_sample_from_trm(trm, n_samples = N_samples)
            output[[sprintf("%s_genotype_transitions", method)]] <- process_samples(sims, 
                n_genes, output = c("transitions"))$transitions
        } 
        else output[[sprintf("%s_genotype_transitions", method)]] <- NA
        
    }

    return(output)
}

## ## Take a sample (a vector), with genotypes as "A, B", etc
## ## and return a vector of frequencies (counts) in the exact same
## ## order as used by MHN
## ## A simple implementation that can be slow when the sample is large

## ## vector of genotypes, total number of genes ->
## ##             counts of all genotypes in same order as used by MHN
## sample_to_pD_order0 <- function(x, ngenes) {
##     x <- gsub("WT", "", x)
##     xs <- strsplit(x, ", ")

##     Data <-  do.call(rbind,
##                    lapply(xs, function(z) as.integer(LETTERS[1:ngenes] %in% z)))

##     ## What follows is from Data.to.pD
##     ## except we do not divide
##     n <- ncol(Data)
##     N <- 2^n
##     Data <- apply(Data, 1, State.to.Int)
##     pD <- tabulate(Data, nbins = N)
##     return(pD)
## }


#' @title Count genotypes 
#' 
#' Take a sample (a vector), with genotypes as "A, B", etc
#' and return a vector of frequencies (counts) in the exact same
#' order as used by MHN
#' A much faster implementation

#' @param x vector of genotypes
#' @param ngenes total number of genes
#' 
#' @return counts of all genotypes in same order as used by MHN
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
    return(tabulate(rep(unname(genot_int), x[, 2]),
                   nbins = 2^ngenes))
}

#' @title Compute p-values for MHN sample vs our sampling
#' 
#' For testing purposes 
#' 
#' Generates a random transitions rate matrix
#' and run samples using the MHN engine and ours
#' Returns the p-value of a chis-sq test to compare if
#' both genotype distribution come from the same initial
#' distribution
#' 
#' @param ngenes Number of genes to sample
#' @param n_samples How many random samples do we want
#' @param B an integer specifying the number of replicates used in the
#'          Monte Carlo test.
#' 
#' @return list of p-values (one for each sample)
pv_one_comp <- function(ngenes, n_samples, B = 2000) {

    theta <- Random.Theta(n = ngenes, sparsity = runif(1, 0.2, 0.8))
    pTh <- Generate.pTh(theta)

    trmx <- cpm_access_genots_paths_w_simplified(theta)

    spop <- suppressMessages(
        population_sample_from_trm(trmx, n_samples = n_samples, cores = 1))
    spop_o <- sample_to_pD_order(spop$obs_events, ngenes)
    pv <- chisq.test(x = spop_o, p = pTh,
                     simulate.p.value = TRUE, B = B)$p.value
    ## this allows to catch cases with small values
    ## if(pv < 0.01)
    ##     browser()
    return(pv)
}

# print(sum(p_values5 < 0.01)/M) ## 0.0108
# print(sum(p_values5 < 0.05)/M) ## 0.0488
# print(sum(p_values5 < 0.005)/M) ## 0.0053, though questionable this can be estimated well with B = 2000
# print(sum(p_values5 < 0.001)/M) ## 0.0012: again, expect a lot of noise here with M and B used.

# save(file = "p_values5_mccbn.RData", p_values5)
# stop()
## This is very slow, because what is slow is simulating the p value
## in the chi-square test
## system.time(
## for(i in 1:M) {
##     cat("\n Doing i = ", i, "\n")
##     p_values[i] <- pv_one_comp(4, 10000)
## }
## )


## M <- 144 ## 213 seconds; so about 52 secs per set.
# M <- 10000
# Ngenes <- 8
# Nsampl <- 50000
# system.time(
#     p_values8 <- unlist(mclapply(1:M,
#                          function(x) pv_one_comp(Ngenes, Nsampl),
#                          mc.cores = detectCores()
#                          ))
# )


# sum(p_values8 < 0.05)/M ## 0.0488
# sum(p_values8 < 0.01)/M ## 0.0108
# sum(p_values8 < 0.005)/M ## 0.0053, though questionable this can be estimated well with B = 2000
# sum(p_values7 < 0.001)/M ## 0.0012: again, expect a lot of noise here with M and B used.

# save(file = "p_values8_test.RData", p_values8)


# M <- 10000 ## 144 in 230 seconds; so about 57 seconds per set. 
# Ngenes <- 7
# Nsampl <- 50000
# system.time(
#     p_values7 <- unlist(mclapply(1:M,
#                          function(x) pv_one_comp(Ngenes, Nsampl, B = 4000),
#                          mc.cores = detectCores()
#                          ))
# )


# sum(p_values7 < 0.05)/M ## 0.0497
# sum(p_values7 < 0.01)/M ## 0.0103
# sum(p_values7 < 0.005)/M ## 0.0054
# sum(p_values7 < 0.001)/M ## 7e-4, but this will not be well estimated?

# save(file = "p_values7_test.RData", p_values7)


# ## both look OK
# hist(p_values8)
# hist(p_values7)
# ## For ks test, recall we are using permutation test, so min p.value is not 0
# ## but 1/(B + 1). Though this minor thing makes no difference
# ks.test(p_values8, "punif", 1/2001, 1) ## p-value = 0.5
# ks.test(p_values7, "punif", 1/4001, 1) ## p-value = 1

# ## From https://stats.stackexchange.com/a/406717
# plot(ecdf(p_values8))
# curve(punif(x, 1/2001, 1), add = TRUE, col = "blue")

# plot(ecdf(p_values7))
# curve(punif(x, 1/4001, 1), add = TRUE, col = "blue")


# if(FALSE) {
#     ## For the hell of it, if we want, run this later
#     ## This will be very slow!!
#     M <- 30000 
#     Ngenes <- 10
#     Nsampl <- 200000
#     system.time(
#         p_values10 <- unlist(mclapply(1:M,
#                                       function(x) pv_one_comp(Ngenes, Nsampl,
#                                                               B = 10000),
#                                       mc.cores = detectCores()
#                                       ))
#     )


#     sum(p_values10 < 0.05)/M ## 
#     sum(p_values10 < 0.01)/M ## 
#     sum(p_values10 < 0.005)/M ## 
#     sum(p_values10 < 0.001)/M ## 

#     save(file = "p_values10_test.RData", p_values10)

#     hist(p_values10)
#     ks.test(p_values10, "punif", 1/10001, 1) 
#     plot(ecdf(p_values10))
#     curve(punif(x, 1/10001, 1), add = TRUE, col = "blue")
# }





## Can run it as
## nohup R --vanilla -f sample_genotypes_from_trm.R &> sample_genotypes_from_trm.Rout & 
## or, with changes in systemd
## loginctl enable-linger
## systemd-run --scope --user R --vanilla --slave -f sample_genotypes_from_trm.R &> sample_genotypes_from_trm.Rout &

################################################################

###  Examples and older stuff

## colnames(theta4) <- rownames(theta4) <- LETTERS[1:4]
## colnames(theta8) <- rownames(theta8) <- LETTERS[1:8]
## colnames(theta10) <- rownames(theta10) <- LETTERS[1:10]

## trm4 <- theta_to_trans_rate_3_SM(theta4)
## trm8 <- theta_to_trans_rate_3_SM(theta8)
## trm10 <- theta_to_trans_rate_3_SM(theta10)
## trm9 <- theta_to_trans_rate_3_SM(theta9)

## save(file = "three_trm.Data", trm4, trm8, trm10, trm9)


## ## Examples

## population_sample_from_trm(trm4, 12)

## population_sample_from_trm(trm8, T_sampling = rexp(23))

## population_sample_from_trm(trm8, T_sampling = rexp(10000))

## uu <- population_sample_from_trm(trm10, T_sampling = rexp(10000))


## ## 2.6 secs
## system.time(null <- population_sample_from_trm(trm8, T_sampling = rexp(10000)) )

## ## 46 secs; 6 seconds in Draco
## system.time(null <- population_sample_from_trm(trm8, T_sampling = rexp(200000)) )


## ## 3.3 secs
## system.time(null <- population_sample_from_trm(trm10, T_sampling = rexp(10000)) )


## ## 47 secs; 9 seconds in Draco
## system.time(null <- population_sample_from_trm(trm10, T_sampling = rexp(200000)))

## ## 7 seconds in Draco
## system.time(null <- population_sample_from_trm(trm9, T_sampling = rexp(200000)))
