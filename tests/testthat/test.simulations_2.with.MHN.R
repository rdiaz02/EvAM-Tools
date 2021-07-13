pwd0 <- getwd()
setwd("../../code_from_what_genotype_next/")
# source("code-all-methods-minimal.R")
source("schill-trans-mat.R")
source("simulations.R")
source("simulations_2.R")
setwd(pwd0)
rm(pwd0)

set.seed(1)

out <- readRDS("../../data/out_cpms.rds")
trm <- out$MHN_transitionRateMatrix

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


rates <- as.vector(trm[trm > 0])
# T_events <- vapply(rates, function(x) rexp(1, x), numeric(1.0)) 
transitions <-  as.data.frame(which(trm > 0, arr.ind = TRUE))
colnames(transitions) <- c("FROM", "TO")
rownames(transitions) <- NULL
RATES <- mapply(function(x, y) trm[x, y]
    , x = transitions$FROM
    , y = transitions$TO)

FROM <- sapply(rownames(trm)[transitions$FROM]
    , str2int)
transitions$FROM <- FROM
TO <- sapply(rownames(trm)[transitions$TO]
    , str2int)
transitions$TO <- TO

probs <- c()
for (i in unique(transitions$FROM)){
    tmp_rates <- RATES[transitions$FROM == i]
    probs <- c(probs, cumsum((tmp_rates / sum(tmp_rates))))
} 

transitions$PROBS <- probs


transitions$GENE_MUTATED <- mapply(
    function(to, from) log2(to - from) + 1 ## Difference gives the int genotype of the single gene mutated
    , transitions$TO, transitions$FROM)

state_transitions <- apply(trm, 1, sum)
ordering <- order(sapply(names(state_transitions), str2int)) 
state_transitions <- state_transitions[ordering]


test_that("Simulations properly handles parameters", {
    expect_error(
        simulate_sample(state_transitions, transitions, T_sampling = -5)
        , "Sampling time should be > 0")
    expect_error(simulate_sample(state_transitions, transitions, T_sampling = "abc")
        , "Time should be a number")
    expect_error(simulate_sample(state_transitions, transitions, n_genes = 4, genotype = -1) 
        , "Genotype out of bounds")
    expect_error(simulate_sample(state_transitions, transitions, n_genes = 4, genotype = 16) 
        , "Genotype out of bounds")
    expect_error(simulate_sample(state_transitions, transitions, n_genes = -4) 
        , "n_genes should be a positive integer")
    ##expect transitions to have the right names in the columns
})

test_that("Output of simulations has the proper attributes", {
    sim <- simulate_sample(state_transitions, transitions, 4)
    expect_type(sim, "list")
    expect_named(sim, c("T_sampling"
        , "T_sum_events"
        , "trajectory"
        , "obs_events"))
    # expect_gte(min(sim$T_sum_events), 0)
    expect_gte(sim$T_sampling, 0)
    expect_equal(sim$trajectory[1], 0)
    expect_equal(sim$trajectory[length(sim$trajectory)], sim$obs_events)
})

test_that("Simulate sample gives the right results", {
    FROM <- sapply(c("WT", "WT", "A", "A, C", "B", "B, E")
        , str2int)
    TO <- sapply(c("A", "B", "A, C", "A, C, D", "B, E", "B, E, F")
        , str2int)
    GENE_MUTATED <- sapply(c("A", "B", "C", "D", "E", "F"),  
        function(x) log2(str2int(x)) + 1) 
    PROBS <- c(0.5, 1, 1, 1, 1, 1)
    tt <- data.frame(FROM = FROM
        , TO = TO
        , GENE_MUTATED = GENE_MUTATED
        , PROBS = PROBS)

    rownames(tt) <- NULL
    all_genotypes <- generate_sorted_genotypes(6)
    ordering <- order(vapply(all_genotypes, str2int, numeric(1)))
    state_transitions <- runif(length(all_genotypes), 0, 20)
    names(state_transitions) <- all_genotypes
    state_transitions <- state_transitions[ordering]
    
    T_sampling <- 0.15
    set.seed(1)
    sim_1 <- simulate_sample(state_transitions, tt, 6, T_sampling)
    expect_equal(sim_1$T_sampling, T_sampling)

    tmp_sum_events <- rep(-1, 6)
    time_events <- c()
    set.seed(1)
    for (i in 1:(length(sim_1$trajectory) - 1)){
        gene_mutated <- log2(sim_1$trajectory[i + 1] - sim_1$trajectory[i]) + 1
        time_events <- c(time_events, rexp(1, state_transitions[ sim_1$trajectory[i] + 1 ]))
        tmp_sum_events[gene_mutated] <- sum(time_events)
        runif(1)
    }
    expect_equal(sim_1$T_sum_events, tmp_sum_events)
    expect_equal(sim_1$obs_events, str2int("A, C"))
    expect_equal(sim_1$trajectory, as.vector(sapply(c("WT", "A", "A, C"), str2int)))
    # browser()

    T_sampling <- 20
    set.seed(1)
    sim_1 <- simulate_sample(state_transitions, tt, 6, T_sampling)
    tmp_sum_events <- rep(-1, 6)
    time_events <- c()
    set.seed(1)
    for (i in 1:(length(sim_1$trajectory) - 1)){
        gene_mutated <- log2(sim_1$trajectory[i + 1] - sim_1$trajectory[i]) + 1
        time_events <- c(time_events, rexp(1, state_transitions[ sim_1$trajectory[i] + 1 ]))
        tmp_sum_events[gene_mutated] <- sum(time_events)
        runif(1)
    }
    expect_equal(sim_1$T_sum_events, tmp_sum_events)
    expect_equal(sim_1$T_sampling, T_sampling)
    expect_equal(sim_1$obs_events, str2int("A, C, D"))
    expect_equal(sim_1$trajectory, as.vector(sapply(c("WT", "A", "A, C", "A, C, D"), str2int)))

    T_sampling <- 20
    set.seed(42)
    sim_1 <- simulate_sample(state_transitions, tt, 6, T_sampling)

    tmp_sum_events <- rep(-1, 6)
    time_events <- c()
    set.seed(42)
    for (i in 1:(length(sim_1$trajectory) - 1)){
        gene_mutated <- log2(sim_1$trajectory[i + 1] - sim_1$trajectory[i]) + 1
        time_events <- c(time_events, rexp(1, state_transitions[ sim_1$trajectory[i] + 1 ]))
        tmp_sum_events[gene_mutated] <- sum(time_events)
        runif(1)
    }
    expect_equal(sim_1$T_sum_events, tmp_sum_events)
    expect_equal(sim_1$T_sampling, T_sampling)
    expect_equal(sim_1$obs_events, str2int("B, E, F"))
    expect_equal(sim_1$trajectory, as.vector(sapply(c("WT", "B", "B, E", "B, E, F"), str2int)))

    T_sampling <- 0.02
    set.seed(42)
    sim_1 <- simulate_sample(state_transitions, tt, 6, T_sampling)

    tmp_sum_events <- rep(-1, 6)
    time_events <- c()
    set.seed(42)
    for (i in 1:(length(sim_1$trajectory) - 1)){
        gene_mutated <- log2(sim_1$trajectory[i + 1] - sim_1$trajectory[i]) + 1
        time_events <- c(time_events, rexp(1, state_transitions[ sim_1$trajectory[i] + 1 ]))
        tmp_sum_events[gene_mutated] <- sum(time_events)
        runif(1)
    }
    expect_equal(sim_1$T_sum_events, tmp_sum_events)
    expect_equal(sim_1$T_sampling, T_sampling)
    expect_equal(sim_1$obs_events, str2int("B"))
    expect_equal(sim_1$trajectory, as.vector(sapply(c("WT", "B"), str2int)))

 })

test_that("Output behaves_properly", {
    #Expect similar frequencies
    expect_error(simulate_population(trm[-1, ], n_samples = 100), "Transition matrix should be squared")
    sim <- simulate_population(trm, n_samples = 100)
    expect_type(sim, "list")
    expect_named(sim, c("T_sampling"
        , "T_sum_events"
        , "trans_table"
        , "trajectory"
        , "obs_events"
        ))
    expect_gte(min(sim$T_sampling), 0)
    expect_length(sim$T_sampling, 100)
    expect_length(sim$obs_events, 100)
    expect_length(sim$T_sum_events, 100 * 4)

    tt <- sim$trans_table
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
