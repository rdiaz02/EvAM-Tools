## Copyright 2022 Ramon Diaz-Uriarte, Pablo Herrera Nieto.

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



#' @title Sample an indivial from a transition rate matrix
#' 
#' @param trm transition rate matrix
#' @param T_sampling Time at which sampling happens.
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


#' @title Sample an individual from a transition rate matrix
#' 
#' We repeat the sums to compute the diagonal and the division
#' If we sample a large number of times, possibly worth it to
#' have those precomputed

#' This is what I call "transition matrix standardized":
#'    Diagonal is passed separately, entries in matrix are probabilities

#' This will only be called after standardization

#' @param trmstd transition rate matrix "standardized",
#' @param diag diagonal of transition rate matrix, time of sampling of a case/individual
#' @param T_sampling Time at which sampling happens.
#' @param ngenots Number of genotypes
#' @param genot_names String array with genotype names
#' 
#' @return sampled genotype, trajectory, and accumulated time
indiv_sample_from_trm_pre <- function(trmstd,
                                      diag,
                                      T_sampling,
                                      ngenots,
                                      genot_names) {
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


#' @title Sample a population from a transition rate matrix
#' 
#' @description Like indiv_sample_from_trm, but for multiple times
#' 
#' @param trm transition rate matrix, number of samples or times of samples,
#' @param n_samples Int with the number of samples to be computed
#' @param T_sampling Time at wich each individual in sample. By default they 
#' are randomly generated from an exponential distribution of rate 1.
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
#' @description Generate trajectories from data simulated from a given model.
#' 
#' @param sim list generated with population_sample_from_trm. Relevant
#' fields are described below
#' $T_sum_events time of events for the mutations of each gene
#' $obs_events data.frame with mutated before the end of the sampling time
#' @param n_genes number of genes observed
#' @param output type of output that we want
#' 
#' @return List with a list of trajectories (the order in which gene mutations
#' are acquired), genotype frequencies and genotypes transition matrix (with
#' counts of how many transitions between each genotype have been observed) 
process_samples <- function(sim, n_genes, gene_names = NULL,
                            output = c("frequencies",
                                       "state_counts",
                                       "transitions")){
    if(is.null(gene_names)) gene_names <- LETTERS[1:n_genes]
    
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
    sorted_genotypes <- vapply(0:(n_states - 1), function(x){
        tmp_genotype <- paste(gene_names[int2binary(x, n_genes) == 1]
            , collapse = ", ")
        tmp_genotype <- ifelse(tmp_genotype == "", "WT", tmp_genotype)
        return(tmp_genotype)
    }, character(1))
    trajectories <- sim$trajectory

    #Calculate frequencies
    if(out_params["frequencies"]){
        frequencies <- sample_to_pD_order(sim$obs_events, n_genes, gene_names)
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
        state_counts <- sample_to_pD_order(unlist(sim$trajectory),
                                           n_genes, gene_names)
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
#' OT already provides the predicted genotypes
#' though not from a transition rate matrix
#' 
#' @param cpm_output Output from calling all_methods2trans_mat
#' @param n_samples Number of samples to generate
#' @param n_genes Number of samples that are in the sample
#' @param methods List of methods that we want to sample
#' 
#' @return modified cpm_outputd including a matrix with genotype transitions
sample_all_CPMs <- function(cpm_output
    , n_samples
    , n_genes
    , gene_names = NULL
    , methods = c("Source", "CBN", "MCCBN", "DBN", "MHN", "HESBCN", "OT")) {
    output <- cpm_output
    ## I have removed OT from the list of CPM to sample
    ## And I have "Source" for a source data type for the web server
    if (is.null(gene_names)) gene_names <- LETTERS[1:n_genes]

    for (method in methods) {
        if (method == "OT") {
            tmp_data <- output$OT_genots_predicted
            genots <- tmp_data[2:(ncol(tmp_data) - 1)]

            genots_2 <- unlist(apply(genots, 1, 
                function(x) paste(names(genots)[x == 1], collapse = ", ")))
            names(genots_2) <- NULL
            genots_2[genots_2 == ""] <- "WT"

            tmp_genotypes_sampled <- sample_to_pD_order(
                sample(genots_2, n_samples, 
                    prob = tmp_data$Prob, replace = TRUE),
                n_genes, gene_names)
            
            output[[sprintf("%s_genotype_freqs", method)]] <-
                data.frame(
                    Genotype = generate_sorted_genotypes(n_genes, gene_names),
                    Counts = tmp_genotypes_sampled
                )
            ## The next one is NOT implicitly available.
            ##   see OT_transition_matrices.org
            output[[sprintf("%s_genotype_transitions", method)]] <- NULL
        } else {
            if (method == "MHN") {
                trm <- output$MHN_transitionRateMatrix
            } else {
                trm <- output[[sprintf("%s_f_graph", method)]]
            }
            if (any(!is.na(trm))) {
                print(sprintf("Running %s", method))
                sims <- population_sample_from_trm(trm, n_samples = n_samples)
                psamples <- process_samples(sims,
                                            n_genes, gene_names,
                                            output = c("transitions",
                                                       "frequencies")
                                            )
                output[[sprintf("%s_genotype_transitions", method)]] <-
                    psamples$transitions
                output[[sprintf("%s_genotype_freqs", method)]] <-
                    psamples$frequencies
            } else {
                output[[sprintf("%s_genotype_transitions", method)]] <- NULL
            }
        }
        }
    return(output)
}

# evamtools_pipeline <- function(data){
#     n_genes <- ncol(data)
#     cpm_output <- all_methods_2_trans_mat(data) 
#     output <- sample_all_CPMs(cpm_output, 10000, n_genes)
#     return(output)
# }


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
sample_to_pD_order <- function(x, ngenes, gene_names = NULL) {
    if(is.null(gene_names)) gene_names <- LETTERS[1:ngenes]
    x <- as.data.frame(table(x), stringsAsFactors = FALSE)
    
    genot_int <- x[, 1]
    genot_int <- gsub("WT", "", genot_int)
    genot_int <- vapply(genot_int,
                        function(z)
                            State.to.Int(as.integer(gene_names %in%
                                                    strsplit(z, ", ")[[1]])),
                        numeric(1))
    ## all_genots <- rep(unname(genot_int), x[, 2])
    return(tabulate(rep(unname(genot_int), x[, 2]),
                   nbins = 2^ngenes))
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
