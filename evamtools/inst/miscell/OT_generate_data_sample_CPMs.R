## Showing the behavior of Oncotree::generate.data
## sample_evam, and random_evam to illustrate error model

######################################################################
##
##     Show how D1 and D2 give same results with estimated
##
######################################################################
library(evamtools)
library(Oncotree)
data(ov.cgh)
ov.tree <- oncotree.fit(ov.cgh)

d1 <- generate.data(1e4, ov.tree, with.errors = TRUE, method = "D1",
                    edge.weights = "estimated")
d1t <- evamtools:::data_to_counts(d1, out = "vector")
d1tp <- d1t/nrow(d1)

d2 <- generate.data(1e4, ov.tree, with.errors = TRUE, method = "D2",
                    edge.weights = "observed")
d2t <- evamtools:::data_to_counts(d2, out = "vector")
d2tp <- d2t/nrow(d1)

d3 <- generate.data(3e5, ov.tree, with.errors = TRUE, method = "D2",
                    edge.weights = "estimated")
d3t <- evamtools:::data_to_counts(d3, out = "vector")
d3tp <- d3t/nrow(d1)

chisq.test(cbind(d1t, d2t), simulate.p.value = TRUE) ## different
chisq.test(cbind(d1t, d3t), simulate.p.value = TRUE) ## similar




######################################################################
##
##     random_evam + sample_evam
##     compared to generate.data
##
######################################################################


## Show how a model with random_evam ot_oncobn_eps = x
## and sample_CPM with obs_noise = y
## gives same predictions as
## and ot.fit with those estimated.weights and
##                epos = x + y - (1/2) x y
##                eneg = y

## ot_oncobn_epos : value of that argument in random_evam
## obs_noise: value of that argument in sample_evam
## gn: number of genes
## N: size of sample
compare_sample_generate <- function(ot_oncobn_epos, obs_noise,
                                    gn = 4, N = 1e4) {
    rx <- random_evam(gn, model = "OT",
                               ot_oncobn_epos = ot_oncobn_epos)

    ## Sample using sample_evam
    s_rx <- sample_evam(rx, N = N, genotype_freqs_as_data = FALSE,
                        obs_noise = obs_noise)
    
    ## Sample using generate.data
    otfit0 <- rx$OT_fit
    new_epos <- rx$OT_fit$eps[["epos"]] + obs_noise -
                       2 * (rx$OT_fit$eps[["epos"]] * obs_noise)
    otfit0$eps <- c(epos = new_epos, eneg = obs_noise)

    gd_otfit0 <- generate.data(N, otfit0,
                               with.errors = TRUE, method = "D1",
                               edge.weights = "estimated")
    f_gd_otfit0 <-
        evamtools:::reorder_to_pD(evamtools:::data_to_counts(gd_otfit0,
                                                             out = "vector"))

    
    ## Compare them
    stopifnot(identical(names(f_gd_otfit0), names(s_rx[[1]])))

    m1 <- cbind(rx = s_rx[[1]], otfit = f_gd_otfit0)
    rs0 <- which(rowSums(m1) == 0)
    if(length(rs0))
        m1_nona <- m1[-rs0, ]
    else
        m1_nona <- m1
    if (nrow(m1_nona) <= 2) browser()
    return(list(cbind(m1),
                chisq.test(m1_nona,
                           simulate.p.value = TRUE)))
}

## Remember: under Ho a 5% prob of p-value <= 0.05, etc.
compare_sample_generate(0.12, .07)
compare_sample_generate(0.02, .13)
compare_sample_generate(0.2, .05)
compare_sample_generate(0.1, .1)








######################################################################
##
##     For a given fitted object, show that
##     distribution oncotree on modified eps + my noise adding
##     vs
##     distribution oncotree on modified eps + D2 noise adding
##     vs
##     D1 generate data on modified eps
##     are all the same
##
######################################################################



## distribution oncotree on modified eps + my noise adding
## vs
## distribution oncotree on modified eps + D2 noise adding
## vs
## (from former)
## D1 generate data on modified eps
dot_noise_gd_2 <- function(of, N = 1e4) {

    ## Emulate random_evam + sample_evam
    of2 <- of
    of2$eps <- c(epos = of$eps[["epos"]],
                 eneg = 0)
    dof2 <- distribution.oncotree(of2,
                                  with.probs = TRUE,
                                  with.errors = TRUE,
                                  edge.weights = "estimated")
    r_sample_odf2 <- sample(1:nrow(dof2), size = N, replace =TRUE,
                          prob = dof2$Prob)
    sample_odf2 <- dof2[r_sample_odf2, -c(1, ncol(dof2))]
    noised_sample_odf2 <- evamtools:::add_noise(as.matrix(sample_odf2),
                                                properr = of$eps[["eneg"]])
    freq_sample_odf2 <-
        evamtools:::reorder_to_pD(evamtools:::data_to_counts(sample_odf2,
                                                             out = "vector"))
    freq_noised_sample_odf2 <-
        evamtools:::reorder_to_pD(evamtools:::data_to_counts(noised_sample_odf2,
                                                             out = "vector"))

    ## generate.data via D1 on modified object
    of3 <- of
    of3$eps <- c(epos = of$eps[["epos"]] + of$eps[["eneg"]] -
                     2 * (of$eps[["epos"]] * of$eps[["eneg"]]),
                 eneg = of$eps[["eneg"]] 
                 )
    gd_of3 <- generate.data(N = N, of3,
                            with.errors = TRUE, method = "D1",
                            edge.weights = "estimated")

    freq_gd_of3 <-
        evamtools:::reorder_to_pD(evamtools:::data_to_counts(gd_of3,
                                                             out = "vector"))
    stopifnot(identical(names(freq_noised_sample_odf2),
                        names(freq_gd_of3)))

    ## Compare the internal error adding procedure of D2 and my noise procedure
    ## But we start from the distribution.oncotree as in D1: using with.errors
    ## So this is not D2

    ## This has eneg = 0
    distr <- dof2

    ran.idx <- sample(1:nrow(distr), size = N, prob = distr$Prob, 
                      replace = TRUE)
    ran.data0 <- distr[ran.idx, 2:of3$nmut]
    rownames(ran.data0) <- 1:N

    ee <- of$eps[["eneg"]]
    
    ran.data <- matrix(rbinom(prod(dim(ran.data0)), size = 1, 
                              prob = ifelse(ran.data0 == 0, ee, 1 - ee)),
                       nrow = nrow(ran.data0),
                       ncol = ncol(ran.data0),
                       dimnames = dimnames(ran.data0))
    
    ran.data2 <- evamtools:::add_noise(as.matrix(ran.data0), ee)

    dxt <- evamtools:::reorder_to_pD(evamtools:::data_to_counts(ran.data,
                                                                out = "vector"))


    ## This is almost the same as freq_noised_sample_odf2 ,
    ## we redo sampling and add noise
    ## But even more similar to dxt, since it shares sampling with dxt
    dr2 <- evamtools:::reorder_to_pD(evamtools:::data_to_counts(ran.data2,
                                                                out = "vector"))

    mycbind <- function(x, y) {
        x1 <- cbind(x, y)
        rr <- which(rowSums(x1) == 0)
        if (length(rr)) x1 <- x1[-rr, ]
        return(x1)
    }
    
    ## comparisons
    m1 <- mycbind(freq_noised_sample_odf2, freq_gd_of3) ## OK, from _1
    m2 <- mycbind(freq_noised_sample_odf2, dxt) ## just use m3
    m3 <- mycbind(freq_gd_of3, dxt)
    m4 <- mycbind(dxt, dr2) ## comparing noise adding: same thing
    m5 <- mycbind(freq_noised_sample_odf2, dr2) ## obvious; basically same thing

    return(list(## m1,
        ## m2,
        ## m3,
        ## m4,
        ## generate random evam dist onc + sample CPMs vs D1 on modified
        dist_onc_eneg_0_plus_my_noise_vs_D1_on_modified = 
        chisq.test(m1, simulate.p.value = TRUE),
        ## generate random evam + sample CPMs vs symmetrized D2 on random evam dist.onc.
        dist_onc_eneg_0_plus_my_noise_vs_symmtrized_D2_on_dist_onc_eneg_0 = 
        chisq.test(m2, simulate.p.value = TRUE),
        ## D1 vs symmetrized D2 on random evam dist.onc.
        symmtrized_D2_on_dist_onc_eneg_0_vs_D1_on_modified = 
        chisq.test(m3, simulate.p.value = TRUE),
        ## symmetrized D2 on random evam dist.onc vs. my noise on random evam dist.onc
        ## Simply comparing the noise adding protocols
        noise_adding_procedures = 
        chisq.test(m4, simulate.p.value = TRUE),
        ## Like repeating the same procedure twice
        dist_onc_eneg_0_plus_my_noise_vs_same_thing =
            chisq.test(m5, simulate.p.value = TRUE),
    ))
}



data(examples_csd)
ex_and <- oncotree.fit(examples_csd$csd$AND$data)
ex_linear <- oncotree.fit(examples_csd$csd$Linear$data)
ex_or <- oncotree.fit(examples_csd$csd$OR$data)
ex_xor <- oncotree.fit(examples_csd$csd$XOR$data)
ex_c1 <- oncotree.fit(examples_csd$csd$c1$data)
ex_c3 <- oncotree.fit(examples_csd$csd$c3$data)
ex_c4c2 <- oncotree.fit(examples_csd$csd$c4c2$data)

ovf <- oncotree.fit(ov.cgh)

dot_noise_gd_2(ovf)
dot_noise_gd_2(ex_linear)
dot_noise_gd_2(ex_or)
dot_noise_gd_2(ex_xor)
dot_noise_gd_2(ex_c1)
dot_noise_gd_2(ex_c3)
dot_noise_gd_2(ex_c4c2)




######################################################################
##
##     For a given fitted object, show that
##     generate.data same as
##     distribution oncotree on modified epos and eneg = 0
##     plus noise addition
##     If possible, we use symmetric noise addition
##     (when epos > eneg)
## 
######################################################################

## prop: how much of the original epos is taken to be noise, when we cannot
## use symmetric noise

dot_noise_gd_4 <- function(of, N = 1e4, prop = 0.5) {

    ## Plain generate data on original object
    oogd <- generate.data(N = N, of,
                          with.errors = TRUE, method = "D1",
                          edge.weights = "estimated")
    freq_oogd <-
        evamtools:::reorder_to_pD(evamtools:::data_to_counts(oogd,
                                                             out = "vector"))

    ## Change epos and eneg and sample with asymmetrical noise
    of3 <- of
    if(of$eps[["epos"]] > of$eps[["eneg"]]) {
        message("Using symmetric noise")
        of3$eps <- c(epos = of$eps[["epos"]] ,
                     eneg = 0)
        noise_p <- noise_n <- of$eps[["eneg"]]
    } else {
        of3$eps <- c(epos = prop * of$eps[["epos"]],
                     eneg = 0)
        noise_p <- (1 - prop) * of$eps[["epos"]]
        noise_n <- of$eps[["eneg"]]
    }

    ## First sample from distribution oncotree, then add asymmetrical error.

    distr <- distribution.oncotree(of3, with.probs = TRUE,
                                   with.errors = TRUE,
                                   edge.weights = "estimated")

    ran.idx <- sample(1:nrow(distr), size = N, prob = distr$Prob,
                      replace = TRUE)
    ran.data0 <- distr[ran.idx, 2:of3$nmut]
    rownames(ran.data0) <- 1:N
   

    ran.data <- matrix(
        rbinom(prod(dim(ran.data0)), size = 1, 
               prob = ifelse(ran.data0 == 0, noise_p, 1 - noise_n)),
        nrow = nrow(ran.data0),
        ncol = ncol(ran.data0),
        dimnames = dimnames(ran.data0))

    dxt <- evamtools:::reorder_to_pD(evamtools:::data_to_counts(ran.data,
                                                                out = "vector"))

    mycbind <- function(x, y) {
        x1 <- cbind(x, y)
        rr <- which(rowSums(x1) == 0)
        if (length(rr)) x1 <- x1[-rr, ]
        return(x1)
    }
    m1 <- mycbind(freq_oogd, dxt)
    return(list(m1,
        dist_onc_eneg_0_plus_my_noise_vs_gen_data_oring =
                chisq.test(m1, simulate.p.value = TRUE)))
}


dot_noise_gd_4(ovf)
dot_noise_gd_4(ex_linear, prop = 1/3)
dot_noise_gd_4(ex_or)
dot_noise_gd_4(ex_xor, prop = 0.9)
dot_noise_gd_4(ex_xor, prop = 0.1)
dot_noise_gd_4(ex_c1)
dot_noise_gd_4(ex_c3, prop = 0.22)
dot_noise_gd_4(ex_c4c2)




######################################################################
##
##     For a given fitted object, if epos > eneg show that
##     generate.data same as
##     distribution oncotree on modified epos and eneg = 0
##     plus symmetrized noise addition
##    
##     But this will not work as it requires
##     that original epos and eneg satisfy
##     0 < 2 * (epos -eneg)/eneg < 1
##     so eneg < epos < (3/2) eneg
######################################################################


dot_noise_gd_3 <- function(of, N = 1e4) {
    if (of$eps[["epos"]] <= of$eps[["eneg"]])
        stop("eneg >= epos")

    ## Plain generate data on original object
    oogd <- generate.data(N = N, of,
                          with.errors = TRUE, method = "D1",
                          edge.weights = "estimated")
    freq_oogd <-
        evamtools:::reorder_to_pD(evamtools:::data_to_counts(oogd,
                                                             out = "vector"))
    
    browser()
    ## Change epos and eneg and sample with symmetrical noise

    of2 <- of
    of2$eps <- c(
        epos = (2 * (of$eps[["epos"]] - of$eps[["eneg"]]))/of$eps[["eneg"]],
        eneg = 0)
    dof2 <- distribution.oncotree(of2,
                                  with.probs = TRUE,
                                  with.errors = TRUE,
                                  edge.weights = "estimated")
    r_sample_odf2 <- sample(1:nrow(dof2), size = N, replace = TRUE,
                            prob = dof2$Prob)
    sample_odf2 <- dof2[r_sample_odf2, -c(1, ncol(dof2))]
    noised_sample_odf2 <- evamtools:::add_noise(as.matrix(sample_odf2),
                                                properr = of$eps[["eneg"]])
    freq_noised_sample_odf2 <-
        evamtools:::reorder_to_pD(evamtools:::data_to_counts(noised_sample_odf2,
                                                             out = "vector"))
    
    mycbind <- function(x, y) {
        x1 <- cbind(x, y)
        rr <- which(rowSums(x1) == 0)
        if (length(rr)) x1 <- x1[-rr, ]
        return(x1)
    }
    m1 <- mycbind(freq_noised_sample_odf2, freq_oogd)
    return(list(m1,
        dist_onc_eneg_0_plus_my_noise_vs_gen_data_oring =
                chisq.test(m1, simulate.p.value = TRUE)))
}




dot_noise_gd_3(ovf)
dot_noise_gd_3(ex_linear)
dot_noise_gd_3(ex_or)
dot_noise_gd_3(ex_xor)
dot_noise_gd_3(ex_c1)
dot_noise_gd_3(ex_c3)
dot_noise_gd_3(ex_c4c2)
