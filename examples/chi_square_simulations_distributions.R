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

K <- 20 # Random scenarios
M <- 320 # Number of replicate
# M <- 2 # Number of replicate
N <- 100000 # Number of simulated samples
# N <- 100 # Number of simulated samples
p_values_1 <- matrix(0, K, M)
p_values_2 <- matrix(0, K, M)

fun_1 <- function(trm){
    sim <- simulate_population(trm, N)
    freqs <- process_simulations(sim, output = c("frequencies"))$frequencies

    ordering <- order(vapply(freqs$Genotype, str2int, numeric(1)))
    chi_sq_test_1 <- chisq.test(freqs$Counts[ordering], p = p_distribution)
    # p_values_1[k, m] <- chi_sq_test$p.value
    chi_sq_test_2 <- chisq.test(freqs$Counts[ordering]/sum(freqs$Counts), p = p_distribution)
    # p_values_2[k, m] <- chi_sq_test$p.value
    return(c(chi_sq_test_1$p.value, chi_sq_test_2$p.value))
}

for (k in 1:K){
    ngenes <- round(runif(1, 2, 12))
    Theta.true <- Random.Theta(n = ngenes, sparsity = 0.50)
    rownames(Theta.true) <- LETTERS[1:ngenes]
    colnames(Theta.true) <- LETTERS[1:ngenes]
    
    p_distribution <- Generate.pTh(Theta.true)
    
    transition_rate_matrix <- theta_to_trans_rate_3(Theta.true)
    print(sprintf("Sample %s", k))
    p_values <- mclapply(list(transition_rate_matrix, transition_rate_matrix)
        , fun_1
        # , MoreArgs = list(trm = transition_rate_matrix)
    )
    p_values <- matrix(unlist(p_values), ncol = 2, byrow = TRUE)

    p_values_1[k, ] <- p_values[, 1]
    p_values_2[k, ] <- p_values[, 2]
    saveRDS(p_values_1, "p_values_1.rds")
    saveRDS(p_values_2, "p_values_2.rds")
}





