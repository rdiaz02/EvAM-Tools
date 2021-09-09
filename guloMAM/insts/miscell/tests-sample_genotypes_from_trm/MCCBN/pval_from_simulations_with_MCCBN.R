library(mccbn)
# setwd("../../")
# source("sample_genotypes_from_trm.R")
# setwd("tests-sample_genotypes_from_trm/MCCBN")
mccbn_vs_comp <- function(ngenes, n_samples, B = 10000) {
    true_p1 <- mccbn::random_poset(ngenes)
    rownames(true_p1) <- colnames(true_p1) <- LETTERS[1:ncol(true_p1)]
    lambda_s <- 1
    lambdas <- runif(ngenes, 1/ngenes*lambda_s, ngenes*lambda_s)

    # Simulations with MCCBN
    simGenotypes <- mccbn::sample_genotypes(n_samples, true_p1,
                                        sampling_param = lambda_s,
                                        lambdas = lambdas)
    spop_mccbn <- apply(simGenotypes$obs_events
        , 1
        , function(x) paste(LETTERS[1:ngenes][x == 1], collapse = ", "))
    spop_mccbn_o <- sample_to_pD_order(spop_mccbn, ngenes)

    #Build trm
    ## Adapted from mccbn-process.R
    df1 <- igraph::as_data_frame(graph_from_adjacency_matrix(true_p1))
    colnames(df1) <- c("From", "To")
    no_parent <- setdiff(colnames(true_p1), df1[, 2])
    dfr <- rbind(
        data.frame(From = "Root", To = no_parent,
                   stringsAsFactors = FALSE),
        df1)
    dfr$edge = paste(dfr[, "From"],
                     dfr[, "To"],
                     sep = " -> ")

    names(lambdas) <- colnames(true_p1)
    dfr$lambda <- lambdas[dfr$To]

    trm <- cpm_access_genots_paths_w_simplified(list(edges = dfr))$weighted_fgraph
    
    # Simulations with our code
    spop <- suppressMessages(
        population_sample_from_trm(trm, n_samples = n_samples, cores = 1))
    spop_o <- sample_to_pD_order(spop$obs_events, ngenes)
    
    genot_by_method <- rbind(spop_o, spop_mccbn_o)
    which_cols_both_0 <- colSums(genot_by_method) == 0
    genot_by_method <- genot_by_method[, !which_cols_both_0]

    pv <- chisq.test(genot_by_method, simulate.p.value = TRUE, B = B)$p.value
    ## this allows to catch cases with small values
    ## if(pv < 0.01)
    ##     browser()
    return(pv)
}

for (i in c(5, 6, 7, 8, 9, 10)){
    print(sprintf("GENES %s", i))
    print(date())
    M <- 10000
    Ngenes <- i
    Nsampl <- 50000

    system.time(
        p_values <- unlist(mclapply(1:M,
                            function(x) {
                                print(sprintf("%s %s", x, date()))
                                mccbn_vs_comp(Ngenes, Nsampl)
                            },
                            mc.cores = detectCores() - 1
                            ))
    )
    save(p_values, file = sprintf("p2_values%s_mccbn.RData", i))
}
