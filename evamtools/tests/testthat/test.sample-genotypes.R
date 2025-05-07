t1 <- Sys.time()

test_that("We get requested output, by the specified means", {
    MCCBN_INSTALLED <- requireNamespace("mccbn", quietly = TRUE)
    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:6]
    out <- suppressMessages(
        evam(Dat1,
             methods = c("CBN", "OT", "OncoBN",
                         "MHN", "HESBCN", "MCCBN")[c(rep(TRUE, 5), MCCBN_INSTALLED)]
             ))

    ## Sample from the predicted genotype frequencies
    ## for all methods in the output out
    outS1 <- sample_evam(out, N = 100)

    expect_true(all(
        c("OT_sampled_genotype_counts",
          "OncoBN_sampled_genotype_counts",
          "CBN_sampled_genotype_counts",
          "MCCBN_sampled_genotype_counts",
          "MHN_sampled_genotype_counts",
          "HESBCN_sampled_genotype_counts")[c(rep(TRUE, 3), MCCBN_INSTALLED, rep(TRUE, 2))]
        %in% names(outS1)))

    expect_true(sum(is.na(unlist(outS1))) == 0)
        
    ## Only CBN and will simulate sampling from the transition
    ## rate matrix. 

    expect_message(outS2 <- sample_evam(out, N = 100, methods = "CBN",
                                         output = "obs_genotype_transitions"),
                    "For the requested output we will need to simulate",
                   fixed = TRUE)
    ## true, but might want to rm the NA components
    ## expect_true(is.na(outS2$CBN_sampled_genotype_counts))
    ## expect_true(is.na(outS2$CBN_state_counts))
    expect_true(sum(is.na(outS2$CBN_obs_genotype_transitions)) == 0)
    ## expect_identical(names(outS2), c("CBN_sampled_genotype_counts",
    ##                                  "CBN_obs_genotype_transitions",
    ##                                  "CBN_state_counts"))
    

    ## No output available for OT
    ## For CBN and MHN simulate from the transition rate matrix
    expect_message(outS3 <- sample_evam(out, N = 100, methods = c("CBN", "OT", "MHN"), 
                         output = c("obs_genotype_transitions",
                                    "state_counts")),
                   "For the requested output we will need to simulate",
                    fixed = TRUE)

    ## true, but might want to rm the NA components
    ## expect_true(is.na(outS3$CBN_sampled_genotype_counts))
    ## expect_true(is.na(outS3$CBN_state_counts))
    expect_true(sum(is.na(outS3$CBN_obs_genotype_transitions)) == 0)
    expect_true(sum(is.na(outS3$MHN_obs_genotype_transitions)) == 0)
    expect_true(sum(is.na(outS3$CBN_state_counts)) == 0)
    expect_true(sum(is.na(outS3$MHN_state_counts)) == 0)

    ## expect_identical(names(outS3), c("CBN_sampled_genotype_counts",
    ##                                  "CBN_obs_genotype_transitions",
    ##                                  "CBN_state_counts"))

    ## OT sampled from the predicted genotype frequencies
    ## No obs_genotype_transitions available for OT
    ## CBN and OT simulate from the transition rate matrix, for consistency

    expect_message(outS4 <- sample_evam(out, N = 100, methods = c("CBN", "OT", "MHN"), 
                         output = c("obs_genotype_transitions",
                                    "sampled_genotype_counts")),
                   "For the requested output we will need to simulate",
                    fixed = TRUE)

    expect_true(sum(is.na(outS4$CBN_obs_genotype_transitions)) == 0)
    expect_true(sum(is.na(outS4$MHN_obs_genotype_transitions)) == 0)
    expect_true(sum(is.na(outS4$CBN_sampled_genotype_counts)) == 0)
    expect_true(sum(is.na(outS4$MHN_sampled_genotype_counts)) == 0)
    expect_true(sum(is.na(outS4$OT_sampled_genotype_counts)) == 0)
})


test_that("Exercise random_evam and sampling and check N is respected", {
    for (i in 1:5) {
        rmhn <- random_evam(model = "MHN", ngenes = 5)
        rcbn <- random_evam(model = "CBN", ngenes = 5,
                            graph_density = 0.5)
        rot <- random_evam(model = "OT", ngenes = 5,
                           graph_density = 0.5)
        robn <- random_evam(model = "OncoBN", ngenes = 5,
                            graph_density = 0.5)
        rhe <- random_evam(model = "HESBCN", ngenes = 5,
                           graph_density = 0.325)
        
        sample_mhn <- sample_evam(rmhn, N = 333)
        sample_cbn <- sample_evam(rcbn, N = 40)
        sample_ot <- sample_evam(rot, N = 55)
        sample_obn <- sample_evam(robn, N = 71)
        sample_he <- sample_evam(rhe, N = 102)

        expect_equal(nrow(sample_mhn$MHN_sampled_genotype_counts_as_data), 333)
        expect_equal(nrow(sample_cbn$CBN_sampled_genotype_counts_as_data), 40)
        expect_equal(nrow(sample_ot$OT_sampled_genotype_counts_as_data), 55)
        expect_equal(nrow(sample_obn$OncoBN_sampled_genotype_counts_as_data), 71)
        expect_equal(nrow(sample_he$HESBCN_sampled_genotype_counts_as_data), 102)

        sample_mhne <- sample_evam(rmhn, N = 333, obs_noise = 0.2)
        sample_cbne <- sample_evam(rcbn, N = 40, obs_noise = 0.1)
    }
})


test_that("standard_rank/order_genots", {
    g1 <- c("F, B", "B, M", "U, A", "C, F", "H, D, T", "E, A, B", "WT")
    expect_identical(standard_rank_genots_1(g1),
                     as.integer(c(3, 4, 2, 5, 7, 6, 1)))
    expect_identical(standard_order_genots_1(g1),
                     as.integer(c(7, 3, 1, 2, 4, 6, 5)))
    expect_identical(canonicalize_genotype_names(g1)[standard_order_genots_1(g1)],
                     c("WT", "A, U", "B, F", "B, M", "C, F", "A, B, E", "D, H, T"))


    g5 <- c("TUV", "PLK, AM1", "RB1, ADP2B2", "ORT, BMN", "CK9, SDN",
            "F, B", "B, M", "I, C", "J, H", "H, I")
    expect_identical(standard_rank_genots_1(g5),
                     as.integer(c(1, 3, 2, 6, 8, 4, 5, 7, 10, 9)))
    expect_identical(standard_order_genots_1(g5),
                     as.integer(c(1, 3, 2, 6, 7, 4, 8, 5, 10, 9)))

    ## Also tested in shiny utils, the code for reorder_genotypes
})


test_that("Miscell corner cases", {
    expect_identical(canonicalize_genotype_names(character(0)),
                     character(0))
    expect_identical(standard_rank_genots_2(character(0), character(0)),
                     integer(0))
    expect_identical(standard_rank_genots_1(character(0)),
                     integer(0))
})

cat("\n Done test.sample-genotypes.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
