test_that("We get requested output, by the specified means", {
    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:6]
    out <- suppressMessages(evam(Dat1,
                                 methods = c("CBN", "OT", "OncoBN",
                                             "MHN", "HESBCN", "MCCBN")))

    ## Sample from the predicted genotype frequencies
    ## for all methods in the output out
    outS1 <- sample_CPMs(out, N = 100)

    expect_true(all(
        c("OT_sampled_genotype_freqs",
          "OncoBN_sampled_genotype_freqs",
          "CBN_sampled_genotype_freqs",
          "MCCBN_sampled_genotype_freqs",
          "MHN_sampled_genotype_freqs",
          "HESBCN_sampled_genotype_freqs") %in% names(outS1)))

    expect_true(sum(is.na(unlist(outS1))) == 0)
        
    ## Only CBN and will simulate sampling from the transition
    ## rate matrix. 

    expect_message(outS2 <- sample_CPMs(out, N = 100, methods = "CBN",
                                         output = "obs_genotype_transitions"),
                    "For the requested output we will need to simulate",
                   fixed = TRUE)
    ## true, but might want to rm the NA components
    ## expect_true(is.na(outS2$CBN_sampled_genotype_freqs))
    ## expect_true(is.na(outS2$CBN_state_counts))
    expect_true(sum(is.na(outS2$CBN_obs_genotype_transitions)) == 0)
    ## expect_identical(names(outS2), c("CBN_sampled_genotype_freqs",
    ##                                  "CBN_obs_genotype_transitions",
    ##                                  "CBN_state_counts"))
    

    ## No output available for OT
    ## For CBN and MHN simulate from the transition rate matrix
    expect_message(outS3 <- sample_CPMs(out, N = 100, methods = c("CBN", "OT", "MHN"), 
                         output = c("obs_genotype_transitions",
                                    "state_counts")),
                   "For the requested output we will need to simulate",
                    fixed = TRUE)

    ## true, but might want to rm the NA components
    ## expect_true(is.na(outS3$CBN_sampled_genotype_freqs))
    ## expect_true(is.na(outS3$CBN_state_counts))
    expect_true(sum(is.na(outS3$CBN_obs_genotype_transitions)) == 0)
    expect_true(sum(is.na(outS3$MHN_obs_genotype_transitions)) == 0)
    expect_true(sum(is.na(outS3$CBN_state_counts)) == 0)
    expect_true(sum(is.na(outS3$MHN_state_counts)) == 0)

    ## expect_identical(names(outS3), c("CBN_sampled_genotype_freqs",
    ##                                  "CBN_obs_genotype_transitions",
    ##                                  "CBN_state_counts"))

    ## OT sampled from the predicted genotype frequencies
    ## No obs_genotype_transitions available for OT
    ## CBN and OT simulate from the transition rate matrix, for consistency

    expect_message(outS4 <- sample_CPMs(out, N = 100, methods = c("CBN", "OT", "MHN"), 
                         output = c("obs_genotype_transitions",
                                    "sampled_genotype_freqs")),
                   "For the requested output we will need to simulate",
                    fixed = TRUE)

    expect_true(sum(is.na(outS4$CBN_obs_genotype_transitions)) == 0)
    expect_true(sum(is.na(outS4$MHN_obs_genotype_transitions)) == 0)
    expect_true(sum(is.na(outS4$CBN_sampled_genotype_freqs)) == 0)
    expect_true(sum(is.na(outS4$MHN_sampled_genotype_freqs)) == 0)
    expect_true(sum(is.na(outS4$OT_sampled_genotype_freqs)) == 0)
})


test_that("Exercise generate_random_evam and sampling", {
    for(i in 1:5) {
    rmhn <- generate_random_evam(model = "MHN", ngenes = 5)
    rcbn <- generate_random_evam(model = "CBN", ngenes = 5,
                             graph_density = 0.5)
    
    sample_mhn <- sample_CPMs(rmhn, N = 1000)
    sample_cbn <- sample_CPMs(rcbn, N = 40)

    sample_mhne <- sample_CPMs(rmhn, N = 1000, obs_noise = 0.2)
    sample_cbne <- sample_CPMs(rcbn, N = 40, obs_noise = 0.1)

    ## d1 <- genotypeCounts_to_data(sample_mhn$MHN_sampled_genotype_freqs, 0)
    ## d2 <- genotypeCounts_to_data(sample_cbn$CBN_sampled_genotype_freqs, 0)

}
})


## Playing around with Oncotree. This is left here, but not part of our
## testing

library(Oncotree)
data(ov.cgh)
ov.tree <- oncotree.fit(ov.cgh)

d1 <- generate.data(3e5, ov.tree, with.errors = TRUE, method = "D1",
                    edge.weights = "estimated")
d1t <- evamtools:::data_to_counts(d1)
d1tp <- d1t/nrow(d1)

d2 <- generate.data(3e5, ov.tree, with.errors = TRUE, method = "D2",
                    edge.weights = "observed")
d2t <- evamtools:::data_to_counts(d2)
d2tp <- d2t/nrow(d1)

d3 <- generate.data(3e5, ov.tree, with.errors = TRUE, method = "D2",
                    edge.weights = "estimated")
d3t <- evamtools:::data_to_counts(d3)
d3tp <- d3t/nrow(d1)

chisq.test(cbind(d1t, d2t)) ## different
chisq.test(cbind(d1t, d3t)) ## similar





ot1 <- oncotree.fit(ov.cgh)
ote <- evam(ov.cgh, methods = "OT")

s_ot1 <- generate.data(1e5, ot1, with.errors = TRUE, method = "D1",
                       edge.weights = "estimated")
f_ot1 <- evamtools:::reorder_to_pD(evamtools:::data_to_counts(s_ot1))

s_ote <- sample_CPMs(ote, N = 1e5, genotype_freqs_as_data = FALSE, obs_noise = 0)

stopifnot(identical(names(f_ot1), names(s_ote[[1]])))

chisq.test(s_ote[[1]], f_ot1)





#############################
r1 <- generate_random_evam(5, model = "OT")

debug(evamtools:::OT_model_2_predict_genots)
## get the otfit, call it otff
## This seems to be it.
en <- 0.03
N <- 1e4
otff2 <- otff
## add .1 to both
otff2$eps <- c(epos = 0.1 + en - 2 * (.1 * en), eneg = en)
gdotff2 <- generate.data(N, otff2, with.errors = TRUE, method = "D1",
                         edge.weights = "estimated")
f_gdotff2 <- evamtools:::reorder_to_pD(evamtools:::data_to_counts(gdotff2))
s_r1 <- sample_CPMs(r1, N = N, genotype_freqs_as_data = FALSE,
                    obs_noise = en)
stopifnot(identical(names(f_gdotff2), names(s_r1[[1]])))
chisq.test(cbind(s_r1[[1]], f_gdotff2))
chisq.test(cbind(s_r1[[1]], f_gdotff2), simulate.p.value = 2000)
cbind(s_r1[[1]], f_gdotff2)





r2 <- generate_random_evam(5, model = "OT", ot_oncobn_epos = 0.02)

debug(evamtools:::OT_model_2_predict_genots)
## get the otfit, call it otff
## This seems to be it.

en <- 0.001
N <- 1e4
ott2 <- ott
## add .1 to both
ott2$eps <- c(epos = unname(ott$eps["epos"] + en - 2 * (ott$eps["epos"] * en)), eneg = en)
gdott2 <- generate.data(N, ott2, with.errors = TRUE, method = "D1",
                         edge.weights = "estimated")
f_gdott2 <- evamtools:::reorder_to_pD(evamtools:::data_to_counts(gdott2))
s_r2 <- sample_CPMs(r2, N = N, genotype_freqs_as_data = FALSE,
                    obs_noise = en)
stopifnot(identical(names(f_gdott2), names(s_r2[[1]])))
chisq.test(cbind(s_r2[[1]], f_gdott2))
chisq.test(cbind(s_r2[[1]], f_gdott2), simulate.p.value = 2000)
cbind(s_r2[[1]], f_gdott2)








cat("\n Done test.sample-genotypes.R \n")






