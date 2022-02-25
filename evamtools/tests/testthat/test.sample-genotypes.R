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

cat("\n Done test.sample-genotypes.R \n")
