test_that("Exercise random_evam with different options", {

    for (mm in c("MHN", "CBN", "HESBCN", "OT", "OncoBN")) {
        ou <- paste0(mm, "_predicted_genotype_freqs")
        o1 <- random_evam(ngenes = 5, model = mm)
        expect_true(!any(is.na(o1[[ou]])))
        expect_true(!is.null(o1[[ou]]))
        }

    for (mm in c("MHN", "CBN", "HESBCN", "OT", "OncoBN")) {
        ng <- sample(2:6, 1)
        gnam <- sample(LETTERS, size = ng)
        ou <- paste0(mm, "_predicted_genotype_freqs")
        o1 <- random_evam(gene_names = gnam, model = mm,
                                   graph_density = 0.7,
                                   cbn_hesbcn_lambda_min = 1/4,
                                   cbn_hesbcn_lambda_max = 6,
                                   hesbcn_probs = c(AND = 0.1,
                                                    OR = 0.6,
                                                    XOR = 0.3),
                                   ot_oncobn_weight_min = 0.2,
                                   ot_oncobn_weight_max = 0.7,
                                   ot_oncobn_epos = 0.1,
                                   oncobn_model = "CBN"
                                   )
        expect_true(!any(is.na(o1[[ou]])))
        expect_true(!is.null(o1[[ou]]))
    }

    
    ng <- sample(2:6, 1)
    gnam <- sample(LETTERS, size = ng)
    ou <- paste0("HESBCN", "_predicted_genotype_freqs")
    expect_error(o1 <- random_evam(gene_names = gnam, model = "HESBCN",
                                            hesbcn_probs = c(AND = 0.1,
                                                UR = 0.6,
                                                XOR = 0.3)
                               ),
                 'identical(sort(names(hesbcn_probs)), c("AND", "OR", "XOR")) is not TRUE',
                 fixed = TRUE)
})


test_that("Test OncoBN thetas in right order", {
    ## Both a test that I am doing it right and that OncoBN is doing it right.
    ## The m2b are about the procedures in evam, but they cannot work
    ## unless m2 do, which is testing OncoBN itself.
    ## And that is how I found a bug in OncoBN
    m1 <- data.frame(From = c("Root", "Root", "A", "B"),
                     To = c("A", "B", "C", "C"),
                     Thetas = c(0.3, 0.4, 0.6, 0.6),
                     Relation = c("Single", "Single", "OR", "OR"))
    
    m1o <- OncoBN_model_2_output(m1, 0)
    expect_true(m1o$OncoBN_predicted_genotype_freqs["C"] == 0)
    expect_equal(unname(m1o$OncoBN_predicted_genotype_freqs["A, C"]),
                 .3 * .6 * (1 - .4))

    m1b <- m1[c(4, 2, 3, 1), ]
    m1bo <- OncoBN_model_2_output(m1b, 0)
    expect_true(m1bo$OncoBN_predicted_genotype_freqs["C"] == 0)
    expect_equal(unname(m1bo$OncoBN_predicted_genotype_freqs["A, C"]),
                 .3 * .6 * (1 - .4))

    
    m2 <- data.frame(From = c("Root", "Root", "A", "B"),
                     To = c("A", "B", "C", "C"),
                     Thetas = c(0.3, 0.4, 0.6, 0.6),
                     Relation = c("Single", "Single", "AND", "AND"))
    
    m2o <- OncoBN_model_2_output(m2, 0)
    expect_true(m2o$OncoBN_predicted_genotype_freqs["C"] == 0)
    expect_equal(unname(m2o$OncoBN_predicted_genotype_freqs["A, C"]),
                 0)
    expect_equal(unname(m2o$OncoBN_predicted_genotype_freqs["A, B"]),
                 .3 * .4 * (1 - 0.6))
    expect_equal(unname(m2o$OncoBN_predicted_genotype_freqs["A, B, C"]),
                 .3 * .6 * .4)
    expect_equal(unname(m2o$OncoBN_predicted_genotype_freqs["A"]),
                 .3 * (1 - .4))
    expect_equal(unname(m2o$OncoBN_predicted_genotype_freqs["B"]),
                 .4 * (1 - .3))
    expect_equal(unname(m2o$OncoBN_predicted_genotype_freqs["C"]),
                 0)
    expect_equal(unname(m2o$OncoBN_predicted_genotype_freqs["B, C"]),
                 0)
    expect_equal(unname(m2o$OncoBN_predicted_genotype_freqs["WT"]),
                 (1 - .3) * (1 - .4))

    
    ## Reorder the model data frame, to check we reorder it internally
    ## as it should
    m2b <- m2[c(4, 2, 3, 1), ]
    m2bo <- OncoBN_model_2_output(m2b, 0)
    expect_true(m2bo$OncoBN_predicted_genotype_freqs["C"] == 0)
    expect_equal(unname(m2bo$OncoBN_predicted_genotype_freqs["A, C"]),
                 0)
    expect_equal(unname(m2bo$OncoBN_predicted_genotype_freqs["A, B"]),
                 .3 * .4 * (1 - 0.6))
    expect_equal(unname(m2bo$OncoBN_predicted_genotype_freqs["A, B, C"]),
                 .3 * .6 * .4)
    expect_equal(unname(m2bo$OncoBN_predicted_genotype_freqs["A"]),
                 .3 * (1 - .4))
    expect_equal(unname(m2bo$OncoBN_predicted_genotype_freqs["B"]),
                 .4 * (1 - .3))
    expect_equal(unname(m2bo$OncoBN_predicted_genotype_freqs["C"]),
                 0)
    expect_equal(unname(m2bo$OncoBN_predicted_genotype_freqs["B, C"]),
                 0)
    expect_equal(unname(m2bo$OncoBN_predicted_genotype_freqs["WT"]),
                 (1 - .3) * (1 - .4))
})
