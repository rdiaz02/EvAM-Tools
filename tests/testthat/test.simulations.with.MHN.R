pwd0 <- getwd()
setwd("../../code_from_what_genotype_next/")
# source("code-all-methods-minimal.R")
source("schill-trans-mat.R")
source("simulations.R")
setwd(pwd0)
rm(pwd0)

out <- readRDS("../../data/out_cpms.rds")
p_distribution <- Generate.pTh(out$MHN_theta)
n_simulation_samples <- 100000
observations <- Finite.Sample(p_distribution, n_simulation_samples) * n_simulation_samples

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


test_that("Simulations properly handles parameters", {
    expect_error(
        simulation_sample(out$MHN_transitionRateMatrix, T_sampling = -5)
        , "Sampling time should be > 0")
    expect_error(simulation_sample(out$MHN_transitionRateMatrix, T_sampling = "abc")
        , "Time should be a number")
    expect_error(simulation_sample(out$MHN_transitionRateMatrix, genotype = c(1, 1, 1, 2)) 
        , "Genotype should only contain 0 or 1")
    expect_error(simulation_sample(out$MHN_transitionRateMatrix, genotype = c(1, 1, 1)) 
        , "Shape Mismatch between genotype and number of genes")
    expect_error(simulation_sample(out$MHN_transitionRateMatrix, sampled_time = -4)
        , "Sampled time should be > 0")
    expect_error(simulation_sample(out$MHN_transitionRateMatrix, T_sum_events = c(1, 2))
        , "Shape mismatch between T_sum_events and number of genes")
    expect_error(simulation_sample(out$MHN_transitionRateMatrix, T_sum_events = c(-4, 1, 1, 3)) 
        , "All times should be positive")
    expect_error(simulation_sample(out$MHN_transitionRateMatrix, T_sum_events = c(2, 1, 4), genotype = c(1, 1, 0, 1),) 
        , "Shape mismatch between T_sum_events and number of genes")
    expect_error(simulation_sample(out$MHN_transitionRateMatrix, T_sum_events = c(2, 1, 4, 3), genotype = c(1, 0, 1, 1),) 
        , "Not-mutated genes can not have sampled time > 0")
    #handles genotyes as integers 
    
})

test_that("Output of simulations has the proper attributes", {
    sim <- simulation_sample(out$MHN_transitionRateMatrix)
    expect_type(sim, "list")
    expect_named(sim, c("T_sampling", "T_sum_events", "obs_events"))
    expect_gte(min(sim$T_sum_events), 0)
    expect_gte(min(sim$T_sampling), 0)
})

test_that("Output is comparable to observations from MHN code", {
    #Expect similar frequencies
    sim <- simulate_population(out$MHN_transitionRateMatrix, n_samples = 10000)
    expect_type(sim, "list")
    expect_named(sim, c("T_sampling", "T_sum_events", "obs_events"))
    expect_gte(min(sim$T_sum_events), 0)
    expect_gte(min(sim$T_sampling), 0)
    expect_length(sim$T_sampling, 100)
    expect_length(sim$obs_events, 100 * 4)
    expect_length(sim$T_sum_events, 100 * 4)

    expect_error( simulate_population(out$MHN_transitionRateMatrix
        , T_sampling = c(0, 0, 0)
        , sampled_time = c(0, 0, 0)
        , T_sum_events = matrix(0, ncol = 4, nrow = 3)
        , genotype = matrix(0, ncol = 3, nrow = 3)
         ), "All parameters must have the same length")
    
    trajs <- process_simulations(sim)
    expect_equal(trajs, sim$trajectory)
    sum(as.integer(unlist(mapply(function(x, y) x != y, sim$trajectory, trajectories))))
    ## COmpare frequencies
    ## How much should I be? 
    browser()
})

## Do or adapt test to run with simulation_sample_2 and simulate_population_2