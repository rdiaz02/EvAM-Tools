source("./utils.R")

#' @title  Process simulations
#' 
#' @description Generate trajectories from simulated data
#' 
#' @param sim list generated with mccbn::sample_genotypes. Relevant
#' fields are described below
#' @param sim$T_sampling list with the time of sampling of each sample
#' @param sim$T_events time of events for the mutations of each gene 
#' given that all parents are satisfied
#' @param sim$T_sum_events time of events for the mutations of each gene
#' @param sim$n number of simulated events
#' @param sim$obs_events data.frame with mutated before the end of the sampling time
process_simulations <- function(sim){
    #Set up
    n_genes <- ncol(sim$T_sum_events)
    n_states <- 2**n_genes

    states <- sapply(0:(n_states - 1), int2str)
    all_genotypes <- sapply(1:(n_genes**2 - 1), int2str)
    sorted_genotypes <- c("WT", all_genotypes[order(sapply(all_genotypes, nchar))])
    int_sorted_genotypes <- as.vector(sapply(sorted_genotypes, str2int))

    trans_table <- as.data.frame(cbind(
        "INT" = 0:(n_states - 1)
        , "STR" = states
        , "BIN" = sapply(0:(n_states - 1), function(x) I(list(int2binary(x, n = n_genes))))
    ))

    #Calculate frequencies
    # browser()
    frequencies <- table(sim$obs_events)
    frequencies <- data.frame(
        Genotype = sorted_genotypes,
        Counts = as.vector(sapply(int_sorted_genotypes, function(x){frequencies[as.character(x)]}))
    )
    rownames(frequencies) <- NULL
    frequencies[is.na(frequencies)] <- 0

    #Calculate trajectories 
    # browser()
    trajectories <- list(rep(NA, length(sim$T_sampling)))
    for(i in 1:length(sim$T_sampling)){
        trajectories[i] <- 
            list(
                sample2trajectory(
                    sim$T_sum_events[i, ], 
                    unlist(trans_table$BIN[trans_table$INT == sim$obs_events[i]])
                    # sim$obs_events[i, ]
                )
        )
    }
    #Calculate transitions

    t <- matrix(0L, nrow = n_states, ncol = n_states)
    for(traj in trajectories){
        traj <- traj + 1
        # str_trajs <- sapply(traj, binary2str)
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

    #Calculate state_counts
    
    # state_counts <- state_count_from_trajectories(trajectories)
    # all_trajs <- matrix(unlist(trajectories), ncol = n_genes, byrow = TRUE)
    state_counts <- table(unlist(trajectories))
    
    state_counts <- data.frame(
        Genotype = sorted_genotypes,
        Counts = as.vector(sapply(int_sorted_genotypes, function(x){state_counts[as.character(x)]}))
    )
    state_counts[is.na(state_counts)] <- 0

    return(list(
        trajectories = trajectories,
        transitions = t,
        state_counts = state_counts,
        frequencies = frequencies
    ))
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
    if(sum(which(sum_time_events < 0) > 0)){
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
#' @param transition_rate_matrix
#' @param T_sampling Numeric > 0. Time at which the sampled is observed.
#' @param sampled_time Numeric > 0. Time at which the observations start
#' @param genotype Vector 0 and 1. Starting genotype 
#' @param T_sum_events Vector numeric. Same length as genotype.
#' Mutation time at which each gene has mutated. 0 if it not mutated.
#' 
#' @return T_sampling Vector with time of observed mutations
simulation_sample <- function(transition_rate_matrix, T_sampling = NA,
    sampled_time = 0, genotype = NULL, T_sum_events = NULL, checks = TRUE){
    
    n_genes <- (ncol(transition_rate_matrix))**0.5

    if (checks){
        if(is.na(T_sampling)) T_sampling <- rexp(1, 1)
        else if (!(is.numeric(T_sampling))) stop("Time should be a number")
        else if(T_sampling <= 0) stop("Sampling time should be > 0")

        if (is.null(genotype)) genotype <- 0
        else if (is.vector(genotype)){
            if (length(genotype) != n_genes) stop("Shape Mismatch between genotype and number of genes")
            else if (!(all(genotype %in% c(0, 1)))) stop("Genotype should only contain 0 or 1")
            else genotype <- binary2int(genotype)
        } else if (genotype > (ncol(transition_rate_matrix) - 1)
            | genotype < 0) stop("This genotype is not supported with the current number of genes")

        if (sampled_time < 0) stop("Sampled time should be > 0") 

        if (is.null(T_sum_events)) T_sum_events <- rep(0, n_genes)
        else if (length(T_sum_events) != n_genes) stop("Shape mismatch between T_sum_events and number of genes")
        else if (any(T_sum_events < 0)) stop("All times should be positive")
        
        if (length(genotype) != length(T_sum_events)) stop("Shape mismatch between genotype and T_sum_events")
        else if (any(T_sum_events[which(genotype == 0)] > 0)){
            stop("Not-mutated genes can not have sampled time > 0")
        }
    }
    # browser()
    tmp_rates <- transition_rate_matrix[as.character(genotype), ]
    tmp_rates <- tmp_rates[tmp_rates > 0]
    trajectory <- c(0)
    while( sampled_time < T_sampling 
        & any(tmp_rates > 0) ){
        time_mutations <- sort.int(vapply(tmp_rates, function(x) {
            rexp(1, x)
        }, numeric(1)))[1]
        # browser()

        # z <- sapply(tmp_rates, function(x) {
        #     rexp(1, x)
        # })

        min_time2mutation <- time_mutations[1]
        new_genotype <- as.integer(names(min_time2mutation)) 
        new_gene_mutated <- which(
            int2binary(new_genotype - genotype) == 1)
        sampled_time <- sampled_time + as.numeric(min_time2mutation)
        if(sampled_time <= T_sampling){
            T_sum_events[new_gene_mutated] <- sampled_time
            genotype <- new_genotype
            trajectory <- c(trajectory, genotype)
            tmp_rates <- transition_rate_matrix[as.character(genotype), ]
            tmp_rates <- tmp_rates[tmp_rates > 0]
        }
    }

    # browser()
    return(list(
        T_sampling = T_sampling, 
        T_sum_events = T_sum_events,
        trajectory = trajectory,
        obs_events = genotype))
}

#' @title Simulates population from a transtion rate matrix
#' 
#' @description Create the tumor development process for a population of patients
#' 
#' @inheritParams simulation_sample
simulate_population <- function(transition_rate_matrix, n_samples = 10, T_sampling = NULL,
    sampled_time = NULL, genotype = NULL, T_sum_events = NULL){
    transition_rate_matrix <- as.matrix(transition_rate_matrix)
    rownames(transition_rate_matrix) <- as.vector(sapply(
        rownames(transition_rate_matrix), function(x) str2int(x)))
    colnames(transition_rate_matrix) <- as.vector(sapply(
        colnames(transition_rate_matrix), function(x) str2int(x)))

    n_genes <- ncol(transition_rate_matrix) ** 0.5

    all_params <- list(T_sampling = T_sampling, sampled_time = sampled_time,
        genotype = genotype, T_sum_events = T_sum_events)
    
    is_param_null <- sapply(all_params, is.null)

    length_params <- sapply(all_params[which(is_param_null == FALSE)], length)
    if (length(length_params) > 0) {
        if (length(unique(length_params)) != 1) stop("All parameters must have the same length")
        else n_samples <- length_params[1]
    }

    for (param in names(all_params[which(is_param_null == TRUE)])){
        if (param %in% c("T_sum_events"))
            all_params[[param]] <- matrix(0, ncol = n_genes, nrow = n_samples)
        else if (param == "T_sampling") 
            all_params[[param]] <- rexp(n_samples, 1)
        else all_params[[param]] <- rep(0, n_samples)
    }

    for (param in c("T_sum_events")) 
        all_params[[param]] <- split(all_params[[param]], row(all_params[[param]]))

    output <- mapply(simulation_sample
        , T_sampling = all_params$T_sampling
        , sampled_time = all_params$sampled_time
        , genotype = all_params$genotype
        , T_sum_events = all_params$T_sum_events
        , MoreArgs = list(transition_rate_matrix = transition_rate_matrix, checks = FALSE)
        , SIMPLIFY = FALSE)

    return(
        list(
            T_sampling = vapply(output, function(x) x$T_sampling, numeric(1)), 
            T_sum_events = t(sapply(output, function(x) x$T_sum_events)),
            trajectory = sapply(output, function(x) x$trajectory),
            obs_events = vapply(output, function(x) x$obs_events, numeric(1)))
        )
}


simulate_sample_2 <- function(T_events, T_sampling
    , genotype
    , transitions
    , n_genes){
    tr <- transitions 
    trajectory <- c(genotype)
    T_sum_events <- rep(0, n_genes)
    accessible_genotypes <- tr$TO[tr$FROM == genotype]
    while (length(accessible_genotypes) > 0){
        accessible_genotypes_idx <- which(tr$FROM == genotype)
        accessible_rates <- T_events[accessible_genotypes_idx]
        # print(accessible_rates)
        # browser()
        next_genotype_idx <- which.min(accessible_rates)
        time2mutation <- accessible_rates[next_genotype_idx]
        new_genotype <- accessible_genotypes[next_genotype_idx]
        gene_mutated <- log2(new_genotype - genotype) + 1 ## Difference gives the int genotype of the single gene mutated
        genotype <- new_genotype
        T_sum_events[gene_mutated] <- sum(T_sum_events) + time2mutation
        trajectory <- c(trajectory, new_genotype)
        accessible_genotypes <- tr$TO[tr$FROM == genotype]
    }

    obs_events <-  as.integer(T_sum_events <= T_sampling)
    trajectory <- trajectory[0:(sum(obs_events) + 1)]

    return(list(
        T_sampling = T_sampling, 
        T_sum_events = T_sum_events,
        trajectory = trajectory,
        obs_events = trajectory[length(trajectory)])
    )
}

simulate_population_2 <- function(transition_rate_matrix
    , n_samples = 10, T_sampling = NULL){

    n_genes <- ncol(transition_rate_matrix) ** 0.5

    ## Build data.frame
    non_empty <- as.data.frame(which(transition_rate_matrix > 0, arr.ind = TRUE))
    colnames(non_empty) <- c("FROM", "TO")
    rownames(non_empty) <- NULL
    RATES <- mapply(function(x, y) transition_rate_matrix[x, y]
        , x = non_empty$FROM
        , y = non_empty$TO)
    non_empty$RATES <- RATES
    
    FROM <- sapply(rownames(transition_rate_matrix)[non_empty$FROM]
        , str2int)
    non_empty$FROM <- FROM
    TO <- sapply(rownames(transition_rate_matrix)[non_empty$TO]
        , str2int)
    non_empty$TO <- TO
    non_empty$NUM_GENES_MUTATED <- sapply(non_empty$FROM
        , function(x) sum(int2binary(x, n = n_genes)))

    number_transitions <- length(RATES)
    T_events <- matrix(0, n_samples, number_transitions)
    for (i in 1:number_transitions){
        T_events[ , i] <- rexp(n_samples, RATES[i])
    }
    T_events_2 <- list(rep(NA, n_samples))
    for (i in 1:n_samples){
        T_events_2[i] <- list(T_events[i, ])
    }
    # T_events_2 <- apply(T_events, 1, function(x)x)
    # browser()
    T_sampling <- rexp(n_samples)

    genotypes <- rep(0, n_samples)
    # genotypes <- c(0, 1, 2, 15)
    # accesible_genotypes <- sapply(genotypes, function(x) as.vector(TO[x == FROM]))
    # rates_idx <- sapply(genotypes, function(x) as.vector(which(FROM == x)))

    output <- mapply(simulate_sample_2
        , T_events = T_events_2
        , T_sampling = T_sampling
        # , sampled_time = all_params$sampled_time
        , genotype = genotypes
        , MoreArgs = list(transitions = non_empty, n_genes = n_genes)
        , SIMPLIFY = FALSE
        )
    # browser()

    return(
        list(
            T_sampling = T_sampling, 
            T_sum_events = t(sapply(output, function(x) x$T_sum_events)),
            trajectory = sapply(output, function(x) x$trajectory),
            obs_events = t(sapply(output, function(x) x$obs_events)),
            T_events = T_events)
        )
    
    }


# freqs <- rep(0, binary2int(c(1,1,1,1)) + 1)
# t <- data.frame(
#     From = integer(), 
#     To = integer(), 
#     Counts = integer(),
#     stringsAsFactors = FALSE)

# for(i in c(1:1000)){
#     sim <- sample(out$MHN_transitionRateMatrix)
#     int_genotype <- binary2int((sim$genotype))
#     freqs[int_genotype + 1] <- freqs[int_genotype + 1] + 1
#     if(length(sim$trajectory) > 1){

#         for(idx in 1:(length(sim$trajectory) - 1)){
#             start_gene <- sim$trajectory[idx]
#             end_gene <- sim$trajectory[idx + 1]
#             x <- t[(t$From == start_gene 
#                     & t$To == end_gene),]$Counts 

#             if((start_gene == 1) & (end_gene == 1)){ browser()}
#             if(length(x) == 0){
#                 t[nrow(t) + 1, ] <- c(start_gene, end_gene, 1)
#             } else {
#                 t[(t$From == start_gene 
#                     & t$To == end_gene),]$Counts <- x + 1
#             }        
#         }
#     }
# }

# # Make plots
# ## Compare Frequencies
# all_freqs <- data.frame(
#     Analytical_freqs = sorted_observations$Abs_Freqs,
#     Simualated_freqs = freqs/sum(freqs))
# rownames(all_freqs) <- sorted_observations$Genotype

# barplot(t(as.matrix(all_freqs)), 
#     beside=TRUE , 
#     legend.text=T,col=c("blue" , "skyblue") ,
#     las = 2, 
#     ylim=c(0, 1) , 
#     ylab="Absolute Frequencies")

# ## HyperTraps
# library(igraph)
# colnames(t) <- c("from", "to", "weight")
# G <- graph_from_data_frame(t
#     , directed = TRUE
#     # , weighted=TRUE
#     )
# A <- as_adjacency_matrix(G, attr = "weight")

# state_as_strings <- sapply(rownames(A), function(x)int2str(x))
# rownames(A) <- state_as_strings
# colnames(A) <- state_as_strings

# freqs2 <- data.frame(
#     Genotype = sapply(c(0:(length(freqs) -1)), function(x)int2str(x)),
#     Freqs = freqs)

# plot_genot_fg(A, db2, freqs2)