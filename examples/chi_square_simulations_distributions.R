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

K <- 100 # Random scenarios
M <- 3200 # Number of replicate
# M <- 8 # Number of replicate
N <- 100000 # Number of simulated samples
# N <- 100 # Number of simulated samples
pval <- matrix(0, K, M)
p_values_2 <- matrix(0, K, M)
n_genes <- rep(0, K)
all_results <- list(rep(0, K * M))
fun_1 <- function(trm){
    sim <- simulate_population(trm, N)
    freqs <- process_simulations(sim, output = c("frequencies"))$frequencies

    ordering <- order(vapply(freqs$Genotype, str2int, numeric(1)))
    chi_sq_test_1 <- chisq.test(freqs$Counts[ordering], p = p_distribution)
    # chi_sq_test_2 <- chisq.test(freqs$Counts[ordering]/sum(freqs$Counts), p = p_distribution)
    # return(c(chi_sq_test_1$p.value, chi_sq_test_2$p.value))
    return(list(
        xi1 = chi_sq_test_1$p.value
        # , xi2 = chi_sq_test_2$p.value
        , sim = sim
        , freqs = freqs
        , ordering = ordering))
}
idx <- 1
p_distributions <- list(rep(0, K))
for (k in 2:5){
    # ngenes <- round(runif(1, 2, 12))
    ngenes <- k
    Theta.true <- Random.Theta(n = ngenes, sparsity = 0.50)
    rownames(Theta.true) <- LETTERS[1:ngenes]
    colnames(Theta.true) <- LETTERS[1:ngenes]
    
    p_distribution <- Generate.pTh(Theta.true)
    p_distributions[[k - 1]] <- p_distribution
    transition_rate_matrix <- theta_to_trans_rate_3(Theta.true)
    print(sprintf("N genes %s", k))
    p_values <- mclapply(sapply(1:M, function(x) list(x = transition_rate_matrix))
        , fun_1, mc.cores = 32
    )
    # browser()
    # p_values <- matrix(unlist(p_values), ncol = 2, byrow = TRUE)

    pval[(k - 1), ] <- as.vector(vapply(p_values, function(x)x$xi1, numeric(1.0)))
    I <- length(p_values)
    for (i in 1:I) {
        print(idx)
        all_results[[idx]] <- p_values[[i]]
        idx <- idx + 1
    }
    # all_results[ (k - 1), ] <- list(xi1 = p_)
    # p_values_2[k, ] <- as.vector(vapply(p_values, function(x)x$xi2), numeric(1.0))
    # n_genes[k] <- ngenes
    saveRDS(p_val, "p_val.rds")
    # saveRDS(p_values_2, "p_values_2.rds")
}

report <- function(idx, ngenes){
    print(sprintf("N genes %s", ngenes))
    idx_2 <- ngenes - 1
    # print(pval[idx_2, ])
    x <- all_results[[(idx_2-1) * 8 + idx]]
    x$freqs["Analytical"] <- floor(p_distributions[[idx_2]] * N)[x$ordering]
    x$freqs["Diff"] <- x$freqs$Counts - x$freqs$Analytical
    x$freqs["Percent"] <- x$freqs$Diff/x$freqs$Analytical
    print(x$freqs[order(abs(x$freqs$Percent)),])
    print(x$xi1)
    print(sum(x$freqs$Diff))
}

report(1, 5)
report(1, 4)
report(1, 3)
report(1, 2)


