pwd0 <- getwd()
setwd("../../code_from_what_genotype_next/")
# source("code-all-methods-minimal.R")
source("schill-trans-mat.R")
source("simulations.R")
setwd(pwd0)
rm(pwd0)

out <- readRDS("../../data/out_cpms.rds")
p_distribution <- Generate.pTh(out$MHN_theta)
n_simulate_samples <- 100000
observations <- Finite.Sample(p_distribution, n_simulate_samples) * n_simulate_samples

edge_transitions <- t(t(out$MHN_trans_mat) %*% diag(observations))
genotypes <- c("WT")
for (i in 1:(ncol(edge_transitions) - 1)){
    genotypes[i + 1] <- paste(LETTERS[which(as.integer(base::intToBits(i)) == 1)], collapse = ", ")
}

observations <- observations[match(rownames(out$MHN_trans_mat), genotypes)]
sorted_observations <- data.frame(Genotype = rownames(out$MHN_trans_mat), Freq = observations)
colnames(edge_transitions) <- sorted_observations$Genotype
rownames(edge_transitions) <- sorted_observations$Genotype
sorted_observations$Abs_Freqs <- sorted_observations$Freq/sum(sorted_observations$Freq)


rates <- as.vector(out$MHN_transitionRateMatrix[out$MHN_transitionRateMatrix > 0])
T_events <- vapply(rates, function(x) rexp(1, x), numeric(1.0)) 
transitions <-  as.data.frame(which(out$MHN_transitionRateMatrix > 0, arr.ind = TRUE))
colnames(transitions) <- c("FROM", "TO")
rownames(transitions) <- NULL
RATES <- mapply(function(x, y) out$MHN_transitionRateMatrix[x, y]
    , x = transitions$FROM
    , y = transitions$TO)


test_that("Simulations properly handles parameters", {
    expect_error(
        simulate_sample(T_events, transitions, T_sampling = -5)
        , "Sampling time should be > 0")
    expect_error(simulate_sample(T_events, transitions, T_sampling = "abc")
        , "Time should be a number")
    expect_error(simulate_sample(T_events, transitions, n_genes = 4, genotype = -1) 
        , "Genotype out of bounds")
    expect_error(simulate_sample(T_events, transitions, n_genes = 4, genotype = 16) 
        , "Genotype out of bounds")
    expect_error(simulate_sample(T_events, transitions, n_genes = -4) 
        , "n_genes should be a positive integer")
    expect_error(simulate_sample(T_events[1:13], transitions, n_genes = -4) 
        , "Transitions and T_events must have the same length")
    # expect_error(simulate_sample(out$MHN_transitionRateMatrix, genotype = c(1, 1, 1)) 
    #     , "Shape Mismatch between genotype and number of genes")
    # expect_error(simulate_sample(out$MHN_transitionRateMatrix, sampled_time = -4)
    #     , "Sampled time should be > 0")
    # expect_error(simulate_sample(out$MHN_transitionRateMatrix, T_sum_events = c(1, 2))
    #     , "Shape mismatch between T_sum_events and number of genes")
    # expect_error(simulate_sample(out$MHN_transitionRateMatrix, T_sum_events = c(-4, 1, 1, 3)) 
    #     , "All times should be positive")
    # expect_error(simulate_sample(out$MHN_transitionRateMatrix, T_sum_events = c(2, 1, 4), genotype = c(1, 1, 0, 1),) 
    #     , "Shape mismatch between T_sum_events and number of genes")
    # expect_error(simulate_sample(out$MHN_transitionRateMatrix, T_sum_events = c(2, 1, 4, 3), genotype = c(1, 0, 1, 1)) 
    #     , "Not-mutated genes can not have sampled time > 0")
    #handles genotyes as integers 
    
})

test_that("Output of simulations has the proper attributes", {
    sim <- simulate_sample(T_events, transitions, 4)
    expect_type(sim, "list")
    expect_named(sim, c("T_sampling"
        , "T_sum_events"
        , "trajectory"
        , "obs_events"))
    expect_gte(min(sim$T_sum_events), 0)
    expect_gte(sim$T_sampling, 0)
    expect_equal(sim$trajectory[1], 0)
    expect_equal(sim$trajectory[length(sim$trajectory)], sim$obs_events)
})

test_that("Output behaves_properly", {
    #Expect similar frequencies
    sim <- simulate_population(out$MHN_transitionRateMatrix, n_samples = 100)
    expect_type(sim, "list")
    expect_named(sim, c("T_sampling"
        , "T_sum_events"
        , "trans_table"
        , "trajectory"
        , "obs_events"
        , "T_events"))
    expect_gte(min(sim$T_sum_events), 0)
    expect_gte(min(sim$T_sampling), 0)
    expect_length(sim$T_sampling, 100)
    expect_length(sim$obs_events, 100)
    expect_length(sim$T_sum_events, 100 * 4)

    for (i in 1:nrow(sim$trans_table)){
        gen1 <- sim$trans_table$FROM[i]
        gen2 <- sim$trans_table$TO[i]
        rate <- sim$trans_table$RATE[i]
        expect_equal(rate, out$MHN_transitionRateMatrix[int2str(gen1), int2str(gen2)])
    }
    
    trajs <- process_simulations(sim)
    expect_equal(trajs$trajectories, sim$trajectory)
    # sum(as.integer(unlist(mapply(function(x, y) x != y, sim$trajectory, trajectories))))
    ## COmpare frequencies
    ## How much should I be? 
    # browser()
})

## Do or adapt test to run with simulate_sample_2 and simulate_population_2