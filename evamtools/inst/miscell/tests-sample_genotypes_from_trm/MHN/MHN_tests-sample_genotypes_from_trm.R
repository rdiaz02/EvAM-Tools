## File for testing the code works.
## the name p_valuesx indicates, with x, number of genes.
## Run in different machines.

## Permutation tests sometimes with B=2000 or B=4000. Not very relevant, but the
## minimum values were observed in the data.

## The permutation tests used code in
## sample_genotypes_from_trm.R
## The files with the output use code in generate_p_values_data_sets.R
## The commented out code shows, for example, how the p_values8 file
## was generated.
## That uses function pv_one_comp

## Using the KS test is not 100% OK, since we have
## permutation tests, and thus p-values are discrete


## Using a chi-square to test for a uniform distribution
my_chi_unif <- function(pv) {
    b <- round((1/min(pv))) - 1
    if(! (b %in% c(2000, 4000, 10000) )) stop("eh??!!")
    seqo <- seq(from = 1, to = b + 1, by = 1)
    seqp <- seqo/(b + 1)
    tt <- table(pv)

    tt_all <- rep(0, length(seqp))
    names(tt_all) <- seqp
    tt_all[names(tt)] <- tt

    chisq.test(tt_all,
               p = rep(1, length(seqp))/length(seqp),
               simulate.p.value = TRUE
               )
}


check_dist <- function(pv) {
    minp <- min(pv)
    cat("prop. < 0.05 : ", sum(pv < 0.05)/length(pv), "\n")
    cat("prop. < 0.01 : ", sum(pv < 0.01)/length(pv), "\n")
    cat("prop. < 0.005: ", sum(pv < 0.005)/length(pv), "\n")
    cat("prop. < 0.001: ", sum(pv < 0.001)/length(pv), "\n")
    cat("\n")
    print(ks.test(pv, "punif", minp, 1))
    cat("\n")
    print(my_chi_unif(pv))
    hist(pv)
    plot(ecdf(pv))
    curve(punif(x, minp, 1), add = TRUE, col = "blue")
}


par(mfrow = c(1, 2))

load("p_values7_test.RData")
check_dist(p_values7)
load("p_values8_test.RData")
check_dist(p_values8)
rm(p_values7, p_values8)

load("p_values7_test_gallotia.RData")
check_dist(p_values7)
load("p_values8_test_gallotia.RData")
check_dist(p_values8)

load("p_values5_test_lacerta.RData")
check_dist(p_values5)
load("p_values6_test_lacerta.RData")
check_dist(p_values6)



