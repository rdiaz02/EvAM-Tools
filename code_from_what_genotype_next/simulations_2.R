source("./utils.R")

#' @title Simulate sample
#' 
#' @description Create the tumor development process for a patient
#' 
#' @param T_events Vector with transition times to all possible genotypes
#' @param T_sampling Numeric > 0. Time at which the sampled is observed.
#' @param genotype Integer >=0, <(n_genes **2 - 1). Starting genotype 
#' @param transitions Dataframe. Lookup table with all posible transitions between genotypes created from a transitions rate matrix
#' @param n_genes Int. number of genes simulates
#' Mutation time at which each gene has mutated. 0 if it not mutated.
#' 
#' @return T_sampling Vector with time of observed mutations
simulate_sample <- function(
    # T_events
    state_transitions
    , transitions
    , n_genes
    , T_sampling  = NULL
    , genotype  = 0){
    if(is.null(T_sampling)) T_sampling <- rexp(1, 1)
    else if (!(is.numeric(T_sampling))) stop("Time should be a number")
    else if(T_sampling <= 0) stop("Sampling time should be > 0")

    if ((genotype < 0) || genotype >= (n_genes ** 2))
        stop("Genotype out of bounds")

    if (n_genes < 0) stop("n_genes should be a positive integer")

    tr <- transitions 
    trajectory <- c(genotype)
    T_sum_events <- rep(-1, n_genes)
    accessible_genotypes <- tr$TO[tr$FROM == genotype]
    accessible_probabilities <- tr$PROBS[tr$FROM == genotype]
    accessible_gene_mutated <- tr$GENE_MUTATED[tr$FROM == genotype]
    T_cum <- rexp(1, state_transitions[genotype + 1])
    while (T_cum < T_sampling 
        && length(accessible_genotypes) > 0
        ){
        prob_new_genotype <- runif(1)
        new_genotype_idx <- which(prob_new_genotype < accessible_probabilities)[1]
        new_genotype <- accessible_genotypes[new_genotype_idx]
        gene_mutated <- accessible_gene_mutated[new_genotype_idx]

        genotype <- new_genotype
        T_sum_events[gene_mutated] <- T_cum
        trajectory <- c(trajectory, new_genotype)

        # Set up for the new round 
        accessible_genotypes <- tr$TO[tr$FROM == genotype]
        accessible_probabilities <- tr$PROBS[tr$FROM == genotype]
        accessible_gene_mutated <- tr$GENE_MUTATED[tr$FROM == genotype]
        time2mutation <- rexp(1, state_transitions[genotype + 1])
        # if(is.na(time2mutation)) browser()
        T_cum <- T_cum + time2mutation
    }

    obs_events <-  as.integer((T_sum_events <= T_sampling) & (T_sum_events > 0))
    trajectory <- trajectory[0:(sum(obs_events) + 1)]
    return(list(
        T_sampling = T_sampling, 
        T_sum_events = T_sum_events,
        trajectory = trajectory,
        obs_events = trajectory[length(trajectory)]) ## Last element
    )
}

#' @title Simulates population from a transtion rate matrix
#' 
#' @description Create the tumor development process for a population of patients
#' 
#' @param trm Transition rate matrix n_genes**2 x n_genes**2 float matrix with rates transitions between all genotype
#' @param n_samples Int > 0
#' @param T_events Vector with transition times to all possible genotypes
#' @param T_sampling Numeric > 0. Time at which the sampled is observed.
#' 
#' @inheritParams simulation_sample
simulate_population <- function(trm
    , n_samples = 10, T_sampling = NULL){


    if (ncol(trm) != nrow(trm)) 
        stop("Transition matrix should be squared")
    n_genes <- log2(ncol(trm)) 
    
    ## Build data.frame
    trans_table <- as.data.frame(which(trm > 0, arr.ind = TRUE))
    colnames(trans_table) <- c("FROM", "TO")
    rownames(trans_table) <- NULL
    RATES <- mapply(function(x, y) trm[x, y]
        , x = trans_table$FROM
        , y = trans_table$TO)
    trans_table$RATES <- RATES
    
    FROM <- vapply(rownames(trm)[trans_table$FROM]
        , str2int, numeric(1))
    trans_table$FROM <- FROM
    TO <- vapply(rownames(trm)[trans_table$TO]
        , str2int, numeric(1))
    trans_table$TO <- TO

    trans_table <- trans_table[order(FROM), ]

    probs <- c()
    for (i in unique(trans_table$FROM)){
        tmp_rates <- trans_table$RATES[trans_table$FROM == i]
        probs <- c(probs, cumsum((tmp_rates / sum(tmp_rates))))
    } 

    trans_table$PROBS <- probs
   
    trans_table$GENE_MUTATED <- mapply(
        function(to, from) log2(to - from) + 1 ## Difference gives the int genotype of the single gene mutated
        , trans_table$TO, trans_table$FROM)

    rownames(trans_table) <- NULL
    
    # trans_table$NUM_GENES_MUTATED <- sapply(trans_table$FROM
    #     , function(x) sum(int2binary(x, n = n_genes)))

    # browser()
    
    # number_transitions <- length(RATES)
    # T_events <- matrix(0, n_samples, number_transitions)
    # for (i in 1:number_transitions){
    #     T_events[ , i] <- rexp(n_samples, rate = trans_table$RATES[i])
    # }
    # T_events_2 <- list(rep(NA, n_samples))
    # for (i in 1:n_samples){
    #     T_events_2[i] <- list(T_events[i, ])
    # }

    state_transitions <- apply(trm, 1, sum)
    ordering <- order(sapply(names(state_transitions), str2int)) 
    state_transitions <- state_transitions[ordering]

    T_sampling <- rexp(n_samples, rate = 1)

    genotypes <- rep(0, n_samples)
    # browser()
    output <- mapply(simulate_sample
        # , T_events = T_events_2
        , T_sampling = T_sampling
        , genotype = genotypes
        , MoreArgs = list(transitions = trans_table
            , state_transitions = state_transitions
            , n_genes = n_genes)
        , SIMPLIFY = FALSE
        )

    return(
        list(
            T_sampling = T_sampling 
            , T_sum_events = t(sapply(output, function(x) x$T_sum_events))
            , trans_table = trans_table
            , trajectory = sapply(output, function(x) x$trajectory)
            , obs_events = t(vapply(output, function(x) x$obs_events, numeric(1)))
            # , T_events = T_events
            )
        )
    }
