#' @title Compute p-values for MHN sample vs our sampling
#' 
#' For testing purposes 
#' 
#' Generates a random transitions rate matrix
#' and run samples using the MHN engine and ours
#' Returns the p-value of a chis-sq test to compare if
#' both genotype distribution come from the same initial
#' distribution
#' 
#' @param ngenes Number of genes to sample
#' @param n_samples How many random samples do we want
#' @param B an integer specifying the number of replicates used in the
#'          Monte Carlo test.
#' 
#' @return list of p-values (one for each sample)
pv_one_comp <- function(ngenes, n_samples, B = 2000) {

    theta <- Random.Theta(n = ngenes, sparsity = runif(1, 0.2, 0.8))
    pTh <- Generate.pTh(theta)

    trmx <- cpm_access_genots_paths_w_simplified(theta)

    spop <- suppressMessages(
        population_sample_from_trm(trmx, n_samples = n_samples, cores = 1))
    spop_o <- sample_to_pD_order(spop$obs_events, ngenes)
    pv <- chisq.test(x = spop_o, p = pTh,
                     simulate.p.value = TRUE, B = B)$p.value
    ## this allows to catch cases with small values
    ## if(pv < 0.01)
    ##     browser()
    return(pv)
}

# print(sum(p_values5 < 0.01)/M) ## 0.0108
# print(sum(p_values5 < 0.05)/M) ## 0.0488
# print(sum(p_values5 < 0.005)/M) ## 0.0053, though questionable this can be estimated well with B = 2000
# print(sum(p_values5 < 0.001)/M) ## 0.0012: again, expect a lot of noise here with M and B used.

# save(file = "p_values5_mccbn.RData", p_values5)
# stop()
## This is very slow, because what is slow is simulating the p value
## in the chi-square test
## system.time(
## for(i in 1:M) {
##     cat("\n Doing i = ", i, "\n")
##     p_values[i] <- pv_one_comp(4, 10000)
## }
## )


## M <- 144 ## 213 seconds; so about 52 secs per set.
# M <- 10000
# Ngenes <- 8
# Nsampl <- 50000
# system.time(
#     p_values8 <- unlist(mclapply(1:M,
#                          function(x) pv_one_comp(Ngenes, Nsampl),
#                          mc.cores = detectCores()
#                          ))
# )


# sum(p_values8 < 0.05)/M ## 0.0488
# sum(p_values8 < 0.01)/M ## 0.0108
# sum(p_values8 < 0.005)/M ## 0.0053, though questionable this can be estimated well with B = 2000
# sum(p_values7 < 0.001)/M ## 0.0012: again, expect a lot of noise here with M and B used.

# save(file = "p_values8_test.RData", p_values8)


# M <- 10000 ## 144 in 230 seconds; so about 57 seconds per set. 
# Ngenes <- 7
# Nsampl <- 50000
# system.time(
#     p_values7 <- unlist(mclapply(1:M,
#                          function(x) pv_one_comp(Ngenes, Nsampl, B = 4000),
#                          mc.cores = detectCores()
#                          ))
# )


# sum(p_values7 < 0.05)/M ## 0.0497
# sum(p_values7 < 0.01)/M ## 0.0103
# sum(p_values7 < 0.005)/M ## 0.0054
# sum(p_values7 < 0.001)/M ## 7e-4, but this will not be well estimated?

# save(file = "p_values7_test.RData", p_values7)


# ## both look OK
# hist(p_values8)
# hist(p_values7)
# ## For ks test, recall we are using permutation test, so min p.value is not 0
# ## but 1/(B + 1). Though this minor thing makes no difference
# ks.test(p_values8, "punif", 1/2001, 1) ## p-value = 0.5
# ks.test(p_values7, "punif", 1/4001, 1) ## p-value = 1

# ## From https://stats.stackexchange.com/a/406717
# plot(ecdf(p_values8))
# curve(punif(x, 1/2001, 1), add = TRUE, col = "blue")

# plot(ecdf(p_values7))
# curve(punif(x, 1/4001, 1), add = TRUE, col = "blue")


# if(FALSE) {
#     ## For the hell of it, if we want, run this later
#     ## This will be very slow!!
#     M <- 30000 
#     Ngenes <- 10
#     Nsampl <- 200000
#     system.time(
#         p_values10 <- unlist(mclapply(1:M,
#                                       function(x) pv_one_comp(Ngenes, Nsampl,
#                                                               B = 10000),
#                                       mc.cores = detectCores()
#                                       ))
#     )


#     sum(p_values10 < 0.05)/M ## 
#     sum(p_values10 < 0.01)/M ## 
#     sum(p_values10 < 0.005)/M ## 
#     sum(p_values10 < 0.001)/M ## 

#     save(file = "p_values10_test.RData", p_values10)

#     hist(p_values10)
#     ks.test(p_values10, "punif", 1/10001, 1) 
#     plot(ecdf(p_values10))
#     curve(punif(x, 1/10001, 1), add = TRUE, col = "blue")
# }





## Can run it as
## nohup R --vanilla -f sample_genotypes_from_trm.R &> sample_genotypes_from_trm.Rout & 
## or, with changes in systemd
## loginctl enable-linger
## systemd-run --scope --user R --vanilla --slave -f sample_genotypes_from_trm.R &> sample_genotypes_from_trm.Rout &

################################################################

###  Examples and older stuff

## colnames(theta4) <- rownames(theta4) <- LETTERS[1:4]
## colnames(theta8) <- rownames(theta8) <- LETTERS[1:8]
## colnames(theta10) <- rownames(theta10) <- LETTERS[1:10]

## trm4 <- theta_to_trans_rate_3_SM(theta4)
## trm8 <- theta_to_trans_rate_3_SM(theta8)
## trm10 <- theta_to_trans_rate_3_SM(theta10)
## trm9 <- theta_to_trans_rate_3_SM(theta9)

## save(file = "three_trm.Data", trm4, trm8, trm10, trm9)


## ## Examples

## population_sample_from_trm(trm4, 12)

## population_sample_from_trm(trm8, T_sampling = rexp(23))

## population_sample_from_trm(trm8, T_sampling = rexp(10000))

## uu <- population_sample_from_trm(trm10, T_sampling = rexp(10000))


## ## 2.6 secs
## system.time(null <- population_sample_from_trm(trm8, T_sampling = rexp(10000)) )

## ## 46 secs; 6 seconds in Draco
## system.time(null <- population_sample_from_trm(trm8, T_sampling = rexp(200000)) )


## ## 3.3 secs
## system.time(null <- population_sample_from_trm(trm10, T_sampling = rexp(10000)) )


## ## 47 secs; 9 seconds in Draco
## system.time(null <- population_sample_from_trm(trm10, T_sampling = rexp(200000)))

## ## 7 seconds in Draco
## system.time(null <- population_sample_from_trm(trm9, T_sampling = rexp(200000)))
