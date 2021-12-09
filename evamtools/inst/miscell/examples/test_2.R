RNGkind("L'Ecuyer-CMRG")

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

fun_1 <- function(dummy, n_genes, myfun, n = 2e5){
    print(sprintf("Start %s", dummy))
    Theta.true <- Random.Theta(n = n_genes, sparsity = runif(1, 0.2, 0.6))
    rownames(Theta.true) <- LETTERS[1:n_genes]
    colnames(Theta.true) <- LETTERS[1:n_genes]
    trm <- theta_to_trans_rate_3_SM(Theta.true)
    p_distribution <- Generate.pTh(Theta.true)
    print(sprintf("Theta %s", dummy))
    sim <- myfun(trm, n)
    freqs <- process_simulations(sim, output = c("frequencies"))$frequencies
    print(sprintf("Sim %s", dummy))
    ordering <- order(vapply(freqs$Genotype, str2int, numeric(1)))
    chi_sq_test <- chisq.test(freqs$Counts[ordering], p = p_distribution)
    print(sprintf("Return %s p_value %s", dummy, chi_sq_test$p.value))
    return(chi_sq_test$p.value)
}

chi_sq <- function(myfun
    , K = 10000
    , N_reps = 2e5
    , N_GENES = c(3, 5, 6, 7, 8, 10, 12, 13)
    , label="chi_sq"){
    
    p_val <- list()

    # for (x in N_GENES) p_val[[as.character(x)]] <- rep(0, K)
    
    core_corr <- data.frame(
        n_genes = 2:16
        , n_cores = c(rep(32, 7), rep(20, 5), rep(10, 3)))

    dir.create(label)
    
    saveRDS(N_GENES, sprintf("%s/n_genes.rds", label))
    for(ngenes in N_GENES){
        print(sprintf("N_GENES %s using %s cores", ngenes, core_corr$n_cores[core_corr$n_genes == ngenes]))
        p_values <- mclapply(
                        1:K
                        , fun_1
                        , n = N_reps
                        , n_genes = ngenes
                        , myfun = myfun
                        , mc.cores = core_corr$n_cores[core_corr$n_genes == ngenes]
                    )
        
        tryCatch({
            p_val[[as.character(ngenes)]] <- as.vector(vapply(p_values
                , function(x){
                    if(is.null(x)) return(-1)
                    return(x)
                }
                , numeric(1.0)))
        }, error = function(x) browser())
        
        saveRDS(p_val, sprintf("%s/p_val.rds", label))

        png(sprintf("./%s/p_values.png", label))
        par(mfrow = c(ceiling(length(names(p_val))/2), 2))
        for (i in names(p_val)) {
            tmp_val <- p_val[[i]]
            hist(tmp_val
                , breaks = 100
                , xlab = ""
                , ylab = ""
                , main = sprintf("N_genes %s", i))
            abline(v = 0.05, col = "blue", lty = 2)
            abline(v = 0.01, col = "red", lty = 2)
        }
        dev.off()

    }
}

K <- 15000 # Random scenarios
N_reps <- 200000
# N_reps <- 1000
N_GENES <- rev(c(3, 5, 6, 7, 8, 10, 12, 13))
N_GENES <- rev(c(3, 5, 6, 7, 8))
# N_GENES <- rev(c(3, 5))

label_2 <- "test_2_sim_2"
chi_sq(simulate_population_2, K, N_reps, N_GENES, label_2)
# label_1 <- "test_2_sim_1"
# chi_sq(simulate_population,   K, N_reps, N_GENES, label_1)






