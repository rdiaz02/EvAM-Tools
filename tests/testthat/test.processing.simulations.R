
pwd0 <- getwd()
setwd("../../code_from_what_genotype_next/")
# source("cbn-process.R")
source("utils.R")
source("simulations.R")
setwd(pwd0)
rm(pwd0)
# library(mccbn)

# set.seed(1)
# true_p1 <- mccbn::random_poset(6)
# true_p1
# lambda_s <- 1
# lambdas <- runif(6, 1/6*lambda_s, 6*lambda_s)
# set.seed(1)
# simGenotypes <- mccbn::sample_genotypes(10, true_p1,
#                                         sampling_param = lambda_s,
#                                         lambdas = lambdas)

sample1 <- c(1, 2, 3, 4, 5)
observed_wt <- c(0, 0, 0, 0, 0)
out1_wt <- unlist(lapply(c("WT"), str2int, sep = ""))
observed1_1 <- c(1, 0, 0, 0, 0)
out1_1 <- unlist(lapply(c("WT", "A"), str2int, sep = ""))
observed1_2 <- c(1, 1, 0, 0, 0)
out1_2 <- unlist(lapply(c("WT", "A", "AB"), str2int, sep = ""))
observed1_3 <- c(1, 1, 1, 0, 0)
out1_3 <- unlist(lapply(c("WT", "A",  "AB", "ABC"), str2int, sep = ""))
observed1_4 <- c(1, 1, 1, 1, 1)
out1_4 <- unlist(lapply(c("WT", "A", "AB", "ABC", "ABCD", "ABCDE"), str2int, sep = ""))

sample2 <- c(5, 4, 3, 2, 1)
observed2_1 <- c(1, 0, 0, 0, 1) ## Warning
out2_1 <- unlist(lapply(c("WT", "A"), str2int, sep = ""))
observed2_2 <- c(0, 0, 0, 0, 1)
out2_2 <- unlist(lapply(c("WT", "E"), str2int, sep = ""))
observed2_3 <- c(0, 0, 1, 1, 1)
out2_3 <- unlist(lapply(c("WT", "E", "ED", "EDC"), str2int, sep = ""))
observed2_4 <- c(1, 1, 1, 1, 1)
out2_4 <- unlist(lapply(c("WT", "E", "DE", "CDE", "BCDE", "ABCDE"), str2int, sep = ""))

sample3 <- c(-5, 4, -3, 2, 1)
observed_error <- c(1, 1, -3, 1, 1)

test_that("A sample is correctly transformed to a trajectory", {
    expect_equal(sample2trajectory(sample1, observed_wt)
        , out1_wt)
    expect_equal(sample2trajectory(sample1, observed1_1)
        , out1_1)
    expect_equal(sample2trajectory(sample1, observed1_2)
        , out1_2)
    expect_equal(sample2trajectory(sample1, observed1_3)
        , out1_3)
    expect_equal(sample2trajectory(sample1, observed1_4)
        , out1_4)
    expect_equal(sample2trajectory(sample2, observed2_2)
        , out2_2)
    expect_equal(sample2trajectory(sample2, observed2_3)
        , out2_3)
    expect_equal(sample2trajectory(sample2, observed2_4)
        , out2_4)
})
    
test_that("Bad observations or times fail", {
    expect_error(sample2trajectory(sample3, observed2_1), "Negative sampling times are not allowed")
    expect_error(sample2trajectory(sample2, observed_error, "Observations should be defined with 0 (not present) or 1 (present)"))
    expect_error(sample2trajectory(c(1, 1), c(1, 1, 1)), "Mismatching sizes of observations and sampling times")
})

simGenotypes <- list()
obs_events <- matrix(c(
      c(0, 0, 0, 0)
    , c(1, 0, 0, 0)
    , c(1, 0, 0, 0)
    , c(1, 0, 1, 0)
    , c(1, 1, 1, 1)
    ), ncol = 4, byrow = TRUE)

T_sum_events <- matrix(c(
      c(1, 2, 3, 4)
    , c(1, 2, 3, 4)
    , c(1, 3, 4, 2)
    , c(2, 3, 1, 4)
    , c(1, 4, 2, 3)
    ), ncol = 4, byrow = TRUE)

sorted_states <- sorted_genotypes(4)
t <- matrix(0L, nrow = 16, ncol = 16)
rownames(t) <- sorted_states
colnames(t) <- sorted_states
t["WT", "A"] <- 3
t["WT", "C"] <- 1
t["C", "A, C"] <- 1
t["A", "A, C"] <- 1
t["A, C", "A, C, D"] <- 1
t["A, C, D", "A, B, C, D"] <- 1

state_count <- data.frame(
    Genotype = c("WT", "A", "C", "A, C", 
        "A, C, D", "A, B, C, D"),
    Counts = c(5, 3, 1, 2, 1, 1)
)
rownames(state_count) <- NULL

freq <- data.frame(
    Genotype = c("WT", "A", "A, C", "A, B, C, D"),
    Counts = c(1, 2, 1, 1)
)
rownames(freq) <- NULL

trajectories <- list(
    as.vector(sapply(c("WT"), str2int))
    , as.vector(sapply(c("WT", "A"), str2int))
    ,as.vector(sapply(c("WT", "A"), str2int))
    ,as.vector(sapply(c("WT", "C", "A, C"), str2int))
    ,as.vector(sapply(c("WT", "A", "A, C", "A, C, D", "A, B, C, D"), str2int))
)

test_that("Trajectories are not generated with bad input",{
    x <- simGenotypes
    expect_error(process_simulations(x), 
        "T_sum_events is missing from your simulations")
    x <- simGenotypes
    x$obs_events <- obs_events
    expect_error(process_simulations(x), 
        "T_sum_events is missing from your simulations")
    x <- simGenotypes
    x$T_sum_events <- T_sum_events
    expect_error(process_simulations(x), 
        "obs_events is missing from your simulations")
})

simGenotypes$obs_events <- obs_events
simGenotypes$T_sum_events <- T_sum_events
test_that("Trajectories are created from simulations data", {
    out_sim <- process_simulations(simGenotypes)
    expect_equal(out_sim$transitions, t)
    expect_equal(sum(sapply(out_sim$trajectories, length) - 1), sum(t))
    state_count_not_zero <- out_sim$state_count[out_sim$state_count$Counts > 0, ]
    rownames(state_count_not_zero) <- NULL
    expect_equal(state_count_not_zero, state_count)
    freq_not_zero <- out_sim$frequencies[out_sim$frequencies$Genotype == freq$Genotype, ]
    rownames(freq_not_zero) <- NULL
    expect_equal(freq_not_zero, freq)
    expect_equal(out_sim$trajectories, trajectories)
})


simGenotypes$obs_events <- apply(simGenotypes$obs_events, 1, binary2int)

test_that("Trajectories are created from simulations data with obs_events as integers", {
    out_sim <- process_simulations(simGenotypes)
    expect_equal(out_sim$transitions, t)
    expect_equal(sum(sapply(out_sim$trajectories, length) - 1), sum(t))
    state_count_not_zero <- out_sim$state_count[out_sim$state_count$Counts > 0, ]
    rownames(state_count_not_zero) <- NULL
    expect_equal(state_count_not_zero, state_count)
    freq_not_zero <- out_sim$frequencies[out_sim$frequencies$Genotype == freq$Genotype, ]
    rownames(freq_not_zero) <- NULL
    expect_equal(freq_not_zero, freq)
    expect_equal(out_sim$trajectories, trajectories)
})