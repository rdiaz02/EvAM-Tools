pwd0 <- getwd()
setwd("../R")
source("simulations.R")
# source("simulations_2.R")
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

fun_1 <- function(x, trm, myfun, n, p_distribution){
    sim <- myfun(trm, n)
    freqs <- process_simulations(sim, output = c("frequencies"))$frequencies

    ordering <- order(vapply(freqs$Genotype, str2int, numeric(1)))
    chi_sq_test <- chisq.test(freqs$Counts[ordering], p = p_distribution)
    return(list(
        xi1 = chi_sq_test$p.value
        , freqs = freqs
        , ordering = ordering)
        )
}

chi_sq <- function(myfun
    , K = 20
    , M = 32
    , N_reps = c(10000)
    , N_GENES = NULL
    , label="chi_sq"
    , max_genes = 9){
    
    idx <- 1
    p_distributions <- list(rep(0, K))
    p_val <- list()

    for (x in N_reps) p_val[[as.character(x)]] <- matrix(0, K, M)
    
    n_genes <- rep(0, K)
    all_results <- list(rep(0, K * M))
    dir.create(label)
    
    if(is.null(N_GENES)) N_GENES <- round(runif(K, 2, max_genes))
    else N_GENES <- rep(N_GENES, floor(K / length(N_GENES)))
    saveRDS(n_genes, sprintf("%s/n_genes.rds", label))

    for (k in 1:K){

        ngenes <- N_GENES[k]

        Theta.true <- Random.Theta(n = ngenes, sparsity = 0.50)
        rownames(Theta.true) <- LETTERS[1:ngenes]
        colnames(Theta.true) <- LETTERS[1:ngenes]
        
        p_distribution <- Generate.pTh(Theta.true)
        p_distributions[[k]] <- p_distribution
        saveRDS(p_distributions, sprintf("%s/p_distributions.rds", label))
        transition_rate_matrix <- theta_to_trans_rate_3(Theta.true)
        print(sprintf("N genes %s", ngenes))
      
        for (N in N_reps){
            print(sprintf("%s \t samples", N))
            p_values <- mclapply(
                1:M
                , fun_1
                , n = N
                , trm = transition_rate_matrix
                , myfun = myfun
                , p_distribution = p_distribution
                , mc.cores = 32
            )
            
            tryCatch({
                p_val[[as.character(N)]][k, ] <- as.vector(vapply(p_values
                    , function(x){
                        if(is.null(x)) pvalue <- -1
                        else pvalue <- x$xi1 
                            
                        return(pvalue)
                    }
                    , numeric(1.0)))
            }, error = function(x) browser())

            null_values <- sum(p_val[[as.character(N)]][k, ] == -1)
            if(null_values > 0 ) print(sprintf("Null values: %s", null_values))

            I <- length(p_values)
            for (i in 1:I) {
                all_results[[idx]] <- p_values[[i]]
                idx <- idx + 1
            }
            print("Size of all_results")
            print(object.size(all_results), units = "MB")
            saveRDS(p_val, sprintf("%s/p_val.rds", label))
            saveRDS(all_results, sprintf("%s/all_results.rds", label))

            png(sprintf("./%s/p_values.png", label))
            par(mfrow = c(length(names(p_val)), 1))
            for (i in names(p_val)) {
                tmp_val <- p_val[[i]][1:k, ] 
                hist(tmp_val
                    , breaks = 100
                    , xlab = ""
                    , ylab = ""
                    , main = sprintf("Samples %s", i))
                abline(v = 0.05, col = "blue", lty = 2)
                abline(v = 0.01, col = "red", lty = 2)
            }
            dev.off()
        }
    }
}

K <- 20 # Random scenarios
M <- 32 # Number of replicate
N_reps <- c(10000, 100000, 200000)
N_GENES <- c(3, 5, 7)

# label <- "sim_1"
# chi_sq(simulate_population, K, M*3, N_reps, N_GENES, label)
label <- "sim_2"
chi_sq(simulate_population_2, K, M*3, N_reps, N_GENES, label)


K <- 20 # Random scenarios
M <- 32 # Number of replicate
N_reps <- c(10000, 100000, 200000)
N_GENES <- c(3, 5, 9)

# label_1 <- "tests_1"
# chi_sq(simulate_population,   K, M, N_reps, N_GENES, label_1)
# label_2 <- "tests_2"
# chi_sq(simulate_population, K, M*3, N_reps, N_GENES, label_2)
# label_3 <- "tests_3"
# chi_sq(simulate_population, K, M*7, N_reps, N_GENES, label_3)


# report <- function(idx, ngenes){
#     print(sprintf("N genes %s", ngenes))
#     idx_2 <- ngenes - 1
#     # print(pval[idx_2, ])
#     x <- all_results[[(idx_2-1) * 8 + idx]]
#     x$freqs["Analytical"] <- floor(p_distributions[[idx_2]] * N)[x$ordering]
#     x$freqs["Diff"] <- x$freqs$Counts - x$freqs$Analytical
#     x$freqs["Percent"] <- x$freqs$Diff/x$freqs$Analytical
#     print(x$freqs[order(abs(x$freqs$Percent)),])
#     print(x$xi1)
#     print(sum(x$freqs$Diff))
# }

# report(1, 5)
# report(1, 4)
# report(1, 3)
# report(1, 2)


