## Trying to generate a reproducible segfault
## Will parallelize later; for now, launch a bunch.

library(evamtools)

i <- 0

while (TRUE) {
    for (model in c(-1, 1, 2, 3, 4)) {
        for (penalty in c(0, 0.1, 1)) {
            i <- i + 1
            length <- ifelse(i %% 3, 2, 3)
            mcpm <- ifelse(i %% 2, "CBN", "MHN")
            the_seed <- round(runif(1, 1, 1e9))
            cat("\n\n\n ***************** i = ", i,
                "; the_seed = ", the_seed,
                "; model = ", model,
                "; penalty = ", penalty,
                "; length = ", length,
                "; mcpm = ", mcpm,
                "*******************\n\n")
            the_N <- runif(1, 50, 500)
            data <- sample_evam((random_evam(6, model = mcpm)),
                                N = the_N)
            data <- data[[paste0(mcpm, "_sampled_genotype_counts_as_data")]]
            out_h_1 <- evam(data,
                            methods = c("HyperTraPS"),
                            hyper_traps_opts = list(length = length,
                                                    model = model,
                                                    seed = the_seed,
                                                    penalty = penalty))
        }
    }
}
