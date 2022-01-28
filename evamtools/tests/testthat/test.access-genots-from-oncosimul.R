## Testing functionality in access_genots_from_oncosimul.R



test_that("Genotype not accessible if no increase in fitness wrt to ancestor", {
    x1 <- c("WT" = 1, "A" = 2, "B" = 1, "A, B" = 2.5)
    ox1 <- evamtools:::genots_2_fgraph_and_trans_mat(x1)
    expect_equal(ox1$accessible_genotypes, c("A" = 2.0, "A, B" = 2.5))

    x2 <- c("WT" = 1, "A" = 2, "B" = 1, "A, B" = 2)
    ox2 <- evamtools:::genots_2_fgraph_and_trans_mat(x2)
    expect_equal(ox2$accessible_genotypes, c("A" = 2.0))

    x3 <- c("WT" = 1, "A" = 2, "B" = 1, "A, B" = 2)
    ox3 <- evamtools:::genots_2_fgraph_and_trans_mat(x3)
    expect_equal(ox3$accessible_genotypes, c("A" = 2.0))

    ## Also test the equality
    x4 <- c("WT" = 1, "A" = 2, "B" = 0.1, "C" = 1.1,
            "A, B" = 0.2, "B, C" = 1.1, "A, C" = 4,
            "A, B, C" = 4)
    ox4 <- evamtools:::genots_2_fgraph_and_trans_mat(x4)
    expect_equal(ox4$accessible_genotypes,
                 c("A" = 2.0, "C" = 1.1, "A, C" = 4.0))
})


test_that("Genotype not accessible if no path to it from WT",  {
    x3 <- c("WT" = 1, "A" = 2, "B" = 0.1, "A, B" = 0.2)
    ox3 <- evamtools:::genots_2_fgraph_and_trans_mat(x3)
    expect_equal(ox3$accessible_genotypes, c("A" = 2.0))

    ## Also test the equality
    x4 <- c("WT" = 1, "A" = 2, "B" = 0.1, "C" = 0.9,
            "A, B" = 0.2, "B, C" = 3, "A, C" = 4,
            "A, B, C" = 4)
    ox4 <- evamtools:::genots_2_fgraph_and_trans_mat(x4)
    expect_equal(ox4$accessible_genotypes, c("A" = 2.0, "A, C" = 4.0))

    x5 <- c("WT" = 1, "A" = 2, "B" = 0.1, "C" = 0.9, "D" = 1.2,
            "A, B" = 2.2, "A, B, C" = 2.4, "A, B, C, D" = 9,
            "B, C" = 3,
            "C, D" = 1.1,
            "B, D" = 1.1,
            "A, D, C" = 5,
            "B, C, D" = 6)
    ox5 <- evamtools:::genots_2_fgraph_and_trans_mat(x5)
    expect_equal(ox5$accessible_genotypes, c("A" = 2.0,
                                             "D" = 1.2,
                                             "A, B" = 2.2,
                                             "A, B, C" = 2.4,
                                             "A, B, C, D" = 9))
})




test_that("Minimal tests for evamtools:::genots_2_fgraph_and_trans_mat
under general fitness landscapes", {
    ## A minimal set of tests for evamtools:::genots_2_fgraph_and_trans_mat
    ## under general fitness landscapes

    x1 <- c(WT = 1, A = 2.5, B = 1.5,
            "A, B" = 2, D = 4,
            "A, B, C" = 3, "A, B, C, D" = 3.5)

    x1o <- evamtools:::genots_2_fgraph_and_trans_mat(x1)

    expect_equal(unname(x1o$transition_matrix)[1, ],
                 c(0, 1.5 / (1.5 + .5 + 3), .5 / (1.5 + .5 + 3),
                   3 / (1.5 + .5 + 3)
                  , 0, 0, 0)
                 )
    expect_equal(unname(x1o$transition_matrix["D", ]), rep(0, 7))
    expect_equal(unname(x1o$transition_matrix["A", ]), rep(0, 7))

    expect_equal(names(x1o$accessible_genotypes), c("A", "B", "D", "A, B",
                                                    "A, B, C", "A, B, C, D"))


    x2 <- c(WT = 2, A = 2.5, B = 1.5,
            "A, B" = 2, D = 4,
            "A, B, C" = 3, "A, B, C, D" = 3.5)

    x2o <- suppressMessages(evamtools:::genots_2_fgraph_and_trans_mat(x2))

    expect_equal(unname(x2o$transition_matrix)[1, ],
                 c(0, 0.5 / (0.5 + 2), 2 / (0.5 + 2))
                 )

    expect_equal(unname(x2o$transition_matrix["D", ]), rep(0, 3))
    expect_equal(unname(x2o$transition_matrix["A", ]), rep(0, 3))

    expect_equal(names(x2o$accessible_genotypes), c("A", "D"))
})

cat("\n Done test.access-genots-from-oncosimul.R \n")
