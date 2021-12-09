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

fun_1 <- function(dummy, n_genes, myfun){
    print(sprintf("Start %s", dummy))
    Theta.true <- Random.Theta(n = n_genes, sparsity = runif(1, 0.2, 0.6))
    rownames(Theta.true) <- LETTERS[1:n_genes]
    colnames(Theta.true) <- LETTERS[1:n_genes]
    trm_time <- as.double(system.time({
        trm <- theta_to_trans_rate_3_SM(Theta.true)
    })["elapsed"])

    sim_10000_time <- as.double(system.time({
        myfun(trm, 10000)
    })["elapsed"])

    sim_100000_time <- as.double(system.time({
        myfun(trm, 100000)
    })["elapsed"])

    return(c(trm_time, sim_10000_time, sim_100000_time))
}

chi_sq <- function(myfun
    , K = 10000
    , N_GENES = c(3, 5, 6, 7, 8, 10, 12, 13)
    , label="chi_sq"){
    
    p_val <- list()

    
    core_corr <- data.frame(
        n_genes = 2:16
        , n_cores = c(rep(32, 7), rep(20, 5), rep(10, 3)))

    dir.create(label)
    
    saveRDS(N_GENES, sprintf("%s/n_genes.rds", label))
    for(ngenes in N_GENES){
        print(sprintf("N_GENES %s using %s cores", ngenes, core_corr$n_cores[core_corr$n_genes == ngenes]))
        my_times <- mclapply(
                        1:K
                        , fun_1
                        , n_genes = ngenes
                        , myfun = myfun
                        , mc.cores = core_corr$n_cores[core_corr$n_genes == ngenes]
                    )

        tryCatch({
            p_val[[as.character(ngenes)]] <- my_times
        }, error = function(x) browser())
        
        saveRDS(p_val, sprintf("%s/times.rds", label))

        # png(sprintf("./%s/p_values.png", label))
        # par(mfrow = c(ceiling(length(names(p_val))/2), 2))
        # for (i in names(p_val)) {
        #     tmp_val <- p_val[[i]]
        #     hist(tmp_val
        #         , breaks = 100
        #         , xlab = ""
        #         , ylab = ""
        #         , main = sprintf("N_genes %s", i))
        #     abline(v = 0.05, col = "blue", lty = 2)
        #     abline(v = 0.01, col = "red", lty = 2)
        # }
        # dev.off()

    }
}

K <- 50 # Random scenarios
N_GENES <- c(3, 5, 7, 9, 11, 13)
label_2 <- "b1"
chi_sq(simulate_population_2, K, N_GENES, label_2)

# x <- readRDS(sprintf("%s/times.rds"))
# idx <- 
# for(i in names(x)){
#     tmp_data <- x[[tmp_data]]
#     tmp_mean <- 
# }


# plot_data <- sapply(x, function(x){
#     tmp_data <- matrivx(unlist(x), ncol = 3, byrow = TRUE)
#     return(c(
#         apply(tmp_data, 1, mean)
#         , apply(tmp_data, 1, std)
#     ))
# })



