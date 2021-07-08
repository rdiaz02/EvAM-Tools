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

test_that("Simulate sample gives the right results", {
    FROM <- sapply(c("WT", "WT", "A", "A, C", "B", "B, E")
        , str2int)
    TO <- sapply(c("A", "B", "A, C", "A, C, D", "B, E", "B, E, F")
        , str2int)
    tt <- data.frame(FROM = FROM, TO = TO)

    T_sampling <- 2
    
    T_events_1 <- c(0.1, 0.2, 0.1, 0.2, 0.1, 0.2)
    sim_1 <- simulate_sample(T_events_1, tt, 6, T_sampling)
    expect_equal(sim_1$T_sampling, T_sampling)
    expect_equal(sim_1$T_sum_events, c(0.1, 0, 0.2, 0.4, 0, 0))
    expect_equal(sim_1$obs_events, str2int("A, C, D"))
    expect_equal(sim_1$trajectory, as.vector(sapply(c("WT", "A", "A, C", "A, C, D"), str2int)))
    
    T_events_2 <- c(0.2, 0.1, 0.1, 0.2, 0.1, 0.2)
    sim_2 <- simulate_sample(T_events_2, tt, 6, T_sampling)
    expect_equal(sim_2$T_sampling, T_sampling)
    expect_equal(sim_2$T_sum_events, c(0, 0.1, 0, 0, 0.2, 0.4))
    expect_equal(sim_2$obs_events, str2int("B, E, F"))
    expect_equal(sim_2$trajectory, as.vector(sapply(c("WT", "B", "B, E", "B, E, F"), str2int)))

    T_sampling <- 0.3
    sim_3 <- simulate_sample(T_events_2, tt, 6, T_sampling)
    expect_equal(sim_3$T_sampling, T_sampling)
    expect_equal(sim_3$T_sum_events, c(0, 0.1, 0, 0, 0.2, 0.4))
    expect_equal(sim_3$obs_events, str2int("B, E"))
    expect_equal(sim_3$trajectory, as.vector(sapply(c("WT", "B", "B, E"), str2int)))


    T_events_1 <- c(0.1, 0.2, 0.1, 0.2, 0.1, 0.2)
    sim_4 <- simulate_sample(T_events_1, tt, 6, T_sampling)
    expect_equal(sim_4$T_sampling, T_sampling)
    expect_equal(sim_4$T_sum_events, c(0.1, 0, 0.2, 0.4, 0, 0))
    expect_equal(sim_4$obs_events, str2int("A, C"))
    expect_equal(sim_4$trajectory, as.vector(sapply(c("WT", "A", "A, C"), str2int)))

})

test_that("Output behaves_properly", {
    #Expect similar frequencies
    expect_error(simulate_population(out$MHN_transitionRateMatrix[-1, ], n_samples = 100), "Transition matrix should be squared")
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

    trm <- out$MHN_transitionRateMatrix
    tt <- simulate_population(trm, n_samples = 100)$trans_table
    expect_equal(nrow(tt), sum(trm > 0))
    for (i in 1:nrow(tt)){
        gen1 <- tt$FROM[i]
        gen2 <- tt$TO[i]
        rate <- tt$RATE[i]
        expect_equal(rate, trm[int2str(gen1), int2str(gen2)])
        expect_equal(sum(int2binary(gen2 - gen1)), 1)
    }
    expect_lte(abs(1 - mean(sim$T_sampling)), 0.05)
    trajs <- process_simulations(sim)
    expect_equal(trajs$trajectories, sim$trajectory)
    
})
