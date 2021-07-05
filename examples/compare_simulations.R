library(mccbn)

pwd0 <- getwd()
setwd("../code_from_what_genotype_next")
source("simulations.R")
source("code-all-methods-minimal.R")
setwd("MHN")
source("UtilityFunctions.R")
source("ModelConstruction.R")
source("Likelihood.R")
source("RegularizedOptimization.R")

setwd("..")
source("../data/toy_datasets.R")
# source("cbn-process.R")
setwd(pwd0)
rm(pwd0)

# I have to do this with several tests, not just with one, choose several of them
# Data set up

out <- readRDS("../data/out_cpms.rds")

for (i in c("AND", "c2", "c4")){
    out <- all_methods_2_trans_mat(all_examples[[i]])
    # out <- out

    ## MHN
    p_distribution <- Generate.pTh(out$MHN_theta)
    n_simulation_samples <- 50000
    observations <- Finite.Sample(p_distribution, n_simulation_samples) * n_simulation_samples

    edge_transitions <- t(t(out$MHN_trans_mat) %*% diag(observations))
    genotypes <- sapply(0:(ncol(edge_transitions) - 1), int2str)

    observations <- observations[match(rownames(out$MHN_trans_mat), genotypes)]
    observations[is.na(observations)] <- 0
    sorted_observations <- data.frame(Genotype = rownames(out$MHN_trans_mat), Freq = observations)
    colnames(edge_transitions) <- sorted_observations$Genotype
    rownames(edge_transitions) <- sorted_observations$Genotype
    sorted_observations$Abs_Freqs <- sorted_observations$Freq/sum(sorted_observations$Freq)

    ## Custom simulations

    sim <- simulate_population(out$MHN_transitionRateMatrix, n_samples = 50000)
    trajs <- process_simulations(sim)

    ## CBN simulations
    # create poset from model
    # genes <- unique(c(out$CBN_model$From, out$CBN_model$To))[-1]
    # poset <- matrix(0L, ncol = length(genes), nrow = length(genes))
    # rownames(poset) <- genes
    # colnames(poset) <- genes
    # for (i in 1:nrow(out$CBN_model)){
    #   tmp_row <- out$CBN_model[i, ]
    #   try(poset[tmp_row$From, tmp_row$To] <- 1
    #     , silent = TRUE)
    # }

    # cbn_sims <- mccbn::sample_genotypes(50000, poset,
    #                                         sampling_param = 1,
    #                                         lambdas = out$CBN_model$final_lambda)
    # cbn_sims$obs_events <- apply(cbn_sims$obs_events, 1, binary2int)
    # cbn_trajs <- process_simulations(cbn_sims)

    # Compare genotype frequencies between MHN and custom simulations
    all_freqs <- data.frame(
        Analytical_freqs = sorted_observations$Abs_Freqs,
        Simulated_MHN_freqs = (trajs$frequencies$Counts/sum(trajs$frequencies$Counts))[match(rownames(out$MHN_trans_mat), trajs$frequencies$Genotype)]
        # Simulated_CBN_freqs = NULL
        )
    rownames(all_freqs) <- sorted_observations$Genotype

    png(sprintf("compare_freqs_%s.png", i))
    barplot(t(as.matrix(all_freqs)), 
        beside=TRUE , 
        legend.text=T,col=c("blue" , "skyblue") ,
        las = 2, 
        # ylim=c(0, 1) , 
        ylab="Absolute Frequencies")
    dev.off()
}
# # Make HyperTraps plots with the resutls from my simulations and from CBN


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
