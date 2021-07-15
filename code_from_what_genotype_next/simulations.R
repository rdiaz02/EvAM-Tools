source("./utils.R")

#' @title Process simulations
#' 
#' @description Generate trajectories from simulated data
#' 
#' @param sim list generated with mccbn::sample_genotypes. Relevant
#' fields are described below
#' @param sim$T_sum_events time of events for the mutations of each gene
#' @param sim$obs_events data.frame with mutated before the end of the sampling time
process_simulations <- function(sim, output = c("frequencies", "trajectories", "state_counts", "transitions")){

    #Checking input
    params <- c("T_sum_events", "obs_events")
    for (i in params){
        if (!(i %in% names(sim))) 
            stop(sprintf("%s is missing from your simulations", i))
    }

    #Checking output variables
    valid_output <- c("frequencies", "trajectories", "state_counts", "transitions")
    out_params <- valid_output %in% output
    names(out_params) <- valid_output

    if (sum(out_params) == 0) stop("Specify valid output")

    not_valid_params <- output[which(!(output %in% valid_output))]
    if (length(not_valid_params) > 0) 
        warning(sprintf("The following parameters cannot be returned: %s"
            , paste(not_valid_params, collapse = ", " )))
    not_valid_params <- not_valid_params[not_valid_params]

    #Set up
    output <- list()
    n_genes <- ncol(sim$T_sum_events)
    n_states <- 2**n_genes

    sorted_genotypes <- generate_sorted_genotypes(n_genes, index.return = TRUE)
    int_sorted_genotypes <- sorted_genotypes$ix
    sorted_genotypes <- sorted_genotypes$x

    trans_table <- as.data.frame(cbind(
        "INT" = int_sorted_genotypes
        , "STR" = sorted_genotypes
        , "BIN" = sapply(int_sorted_genotypes, function(x) I(list(int2binary(x, n = n_genes))))
    ))

    if (length(sim$obs_events)/n_genes == nrow(sim$T_sum_events)) {
        str_obs_events <- apply(sim$obs_events
            , 1
            , function(x) paste(x, collapse = ""))
        str_bin_genotypes <- vapply(trans_table$BIN
            , function(x) paste(x, collapse = "")
            , character(1))
        sim$obs_events <- as.vector(unlist(sapply(str_obs_events
            , function(x) trans_table$INT[which(x == str_bin_genotypes)]
        )))
    }

    #Calculate frequencies
    if(out_params["frequencies"]){
        frequencies <- table(sim$obs_events)
        frequencies <- data.frame(
            Genotype = sorted_genotypes,
            Counts = as.vector(vapply(int_sorted_genotypes
                , function(x) {frequencies[as.character(x)]}
                , numeric(1)))
        )
        rownames(frequencies) <- NULL
        frequencies[is.na(frequencies)] <- 0

        output$frequencies <- frequencies
    }

    #Calculate trajectories 
    if(out_params["trajectories"] ||
        out_params["state_counts"] ||
        out_params["transitions"]){

        trajectories <- list(rep(NA, nrow(sim$T_sum_events)))
        for(i in 1:nrow(sim$T_sum_events)){
            trajectories[i] <- 
                list(
                    sample2trajectory(
                        sim$T_sum_events[i, ], 
                        unlist(trans_table$BIN[trans_table$INT == sim$obs_events[i]])
                    )
                )
        }
        
        if(out_params["trajectories"]) 
            output$trajectories <- trajectories
    }

    #Calculate transitions
    if(out_params["transitions"]){
        t <- matrix(0L, nrow = n_states, ncol = n_states)
        for(traj in trajectories){
            traj <- traj + 1
            if(length(traj) > 1){
                for(i in 1:(length(traj) - 1)){
                    t[traj[i], traj[i + 1]] <-
                        t[traj[i], traj[i + 1]] + 1
                }
            }
        } 

        t <- t[int_sorted_genotypes + 1, int_sorted_genotypes + 1]
        rownames(t) <- sorted_genotypes
        colnames(t) <- sorted_genotypes

        output$transitions <- t
    }

    #Calculate state_counts
    if(out_params["state_counts"]){
        state_counts <- table(unlist(trajectories))
        
        state_counts <- data.frame(
            Genotype = sorted_genotypes,
            Counts = as.vector(sapply(int_sorted_genotypes, function(x){state_counts[as.character(x)]}))
        )
        state_counts[is.na(state_counts)] <- 0

        output$state_counts <- state_counts
    }

    return(output)
}

#' @title Sample to trajectory
#' 
#' @description Transforms a single sample of observed event 
#' in a trajectory of genotypes
#' 
#' @param sum_time_event Positive vector. Times of appeareance of every mutation
#' @param obs_event Vector of 0 and 1. Determines which mutation are present at the time of sampling. 
#' @return List of ordered genotypes in binary format starting from the WT
sample2trajectory <- function(sum_time_events, obs_events){
    if(sum(which(sum_time_events < -1) > 0)){
        stop("Negative sampling times are not allowed")
    }

    if(!all(obs_events %in% c(0, 1))){
        stop("Observations should be defined with 0 (not present) or 1 (present)")
    }

    if(length(sum_time_events) != length(obs_events)){
        stop("Mismatching sizes of observations and sampling times")
    }

    mutated_genes <- (1:length(obs_events))[obs_events == 1] 
    gene_order <- mutated_genes[order(sum_time_events[mutated_genes])]
    # by_time_gene_order <- order(sum_time_events)[1:sum(obs_events)]

    trajectory <- c(0)
    if(length(gene_order) == 0) return(trajectory)
    
    # if(any(gene_order != by_time_gene_order)){
    #     warning("Observations do not respect the sampling time")
    # }

    base_genotype <- rep(0, length(sum_time_events))
    for(i in gene_order){
        base_genotype[i] <- 1
        trajectory <- c(trajectory, binary2int(base_genotype))
    }
    return(trajectory)
}



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
simulate_sample <- function(T_events
    , transitions
    , n_genes
    , T_sampling  = NULL
    , genotype  = 0){
    # browser()
    if(is.null(T_sampling)) T_sampling <- rexp(1, 1)
    else if (!(is.numeric(T_sampling))) stop("Time should be a number")
    else if(T_sampling <= 0) stop("Sampling time should be > 0")

    if ((genotype < 0) || genotype >= (n_genes ** 2))
        stop("Genotype out of bounds")

    if (nrow(transitions) != length(T_events))
        stop("Transitions and T_events must have the same length")

    if (n_genes < 0) stop("n_genes should be a positive integer")

    tr <- transitions 
    trajectory <- c(genotype)
    T_sum_events <- rep(0, n_genes)
    accessible_genotypes <- tr$TO[tr$FROM == genotype]
    accessible_gene_mutated <- tr$GENE_MUTATED[tr$FROM == genotype]
    T_cum <- 0
    while (length(accessible_genotypes) > 0){
        accessible_genotypes_idx <- which(tr$FROM == genotype)
        accessible_rates <- T_events[accessible_genotypes_idx]
        next_genotype_idx <- which.min(accessible_rates)
        time2mutation <- accessible_rates[next_genotype_idx]
        T_cum <- T_cum + time2mutation
        if (T_cum > T_sampling) break
        new_genotype <- accessible_genotypes[next_genotype_idx]
        gene_mutated <- accessible_gene_mutated[next_genotype_idx]
        # gene_mutated <- log2(new_genotype - genotype) + 1 ## Difference gives the int genotype of the single gene mutated
        genotype <- new_genotype
        T_sum_events[gene_mutated] <- T_cum
        trajectory <- c(trajectory, new_genotype)
        # browser()
        accessible_genotypes <- tr$TO[tr$FROM == genotype]
        accessible_gene_mutated <- tr$GENE_MUTATED[tr$FROM == genotype]
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
    trans_table$NUM_GENES_MUTATED <- vapply(trans_table$FROM
        , function(x) sum(int2binary(x, n = n_genes)), numeric(1))
    ##
    trans_table <- trans_table[order(FROM), ]
    trans_table$GENE_MUTATED <- mapply(
        function(to, from) log2(to - from) + 1 ## Difference gives the int genotype of the single gene mutated
        , trans_table$TO, trans_table$FROM)

    rownames(trans_table) <- NULL
   
    
    number_transitions <- length(RATES)
    T_events <- matrix(0, n_samples, number_transitions)
    for (i in 1:number_transitions){
        T_events[ , i] <- rexp(n_samples, rate = trans_table$RATES[i])
    }
    T_events_2 <- list(rep(NA, n_samples))
    for (i in 1:n_samples){
        T_events_2[i] <- list(T_events[i, ])
    }

    T_sampling <- rexp(n_samples, rate = 1)

    output <- mapply(simulate_sample
        , T_events = T_events_2
        , T_sampling = T_sampling
        , MoreArgs = list(
            transitions = trans_table
            , n_genes = n_genes
        )
        , SIMPLIFY = FALSE
        # 
    )

    return(
        list(
            T_sampling = T_sampling, 
            T_sum_events = t(sapply(output, function(x) x$T_sum_events)),
            trans_table = trans_table,
            trajectory = sapply(output, function(x) x$trajectory),
            obs_events = t(vapply(output, function(x) x$obs_events, numeric(1))),
            T_events = T_events)
        )
    }


  #' @title Simulate sample
#' 
#' @description Create the tumor development process for a patient
#' 
#' @param state_transitions Vector with transition rate from all possible genotypes. -sum of the diagonal of the trm
#' @param transitions Dataframe. Lookup table with all posible transitions between genotypes created from a transitions rate matrix
#' @param n_genes Int. number of genes simulates
#' @param T_sampling Numeric > 0. Time at which the sampled is observed.
#' @param genotype Integer >=0, <(n_genes **2 - 1). Starting genotype 
#' Mutation time at which each gene has mutated. 0 if it not mutated.
#' 
#' @return T_sampling Vector with time of observed mutations
#' @return T_sum_events: times at which each gene is mutated
#' @return trans_table: dataframe with FROM, TO, RATE, GENE_MUTATED, PROBS
#' @return trajectory Int vector with with genotypes 
#' @return obs_events. Int. Last genotype observed.
simulate_sample_2 <- function(
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
    T_sum_events <- rep(0, n_genes)
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
#' @param T_sampling Numeric > 0. Time at which the sampled is observed.
#' 
#' @return List with
#' @return T_sampling: vector of number 
#' @return T_sum_events: times at which each gene is mutated
#' @return trans_table: dataframe with FROM, TO, RATE, GENE_MUTATED, PROBS
#' @return trajectory Int vector with with genotypes 
#' @return obs_events. Int. Last genotype observed.
simulate_population_2 <- function(trm
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
    
    state_transitions <- apply(trm, 1, sum)
    ordering <- order(sapply(names(state_transitions), str2int)) 
    state_transitions <- state_transitions[ordering]

    T_sampling <- rexp(n_samples, rate = 1)

    output <- mapply(simulate_sample_2
        , T_sampling = T_sampling
        , MoreArgs = list(
            transitions = trans_table
            , state_transitions = state_transitions
            , n_genes = n_genes
        )
        , SIMPLIFY = FALSE
        )

    return(
        list(
            T_sampling = T_sampling 
            , T_sum_events = t(sapply(output, function(x) x$T_sum_events))
            , trans_table = trans_table
            , trajectory = sapply(output, function(x) x$trajectory)
            , obs_events = t(vapply(output, function(x) x$obs_events, numeric(1)))
            )
        )
    }
  
