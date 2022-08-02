t1 <- Sys.time()


test_that("data_to_counts correct output, including WT", {

    m00 <- matrix(nrow = 12, ncol = 0)
    m0_2c <- matrix(0, nrow = 4, ncol = 2)
    m0_2c_named <- m0_2c
    colnames(m0_2c_named) <- c("cucu", "coco")
    m0_1c_named <- m0_2c_named[, 1, drop = FALSE]
    m0_1c <- m0_2c[, 1, drop = FALSE]
    m1_unnamed <- m0_2c
    m1_unnamed[3, 1] <- 1
    m1_nammed <- m1_unnamed
    colnames(m1_nammed) <- c("aei", "coco")

    out_wt <- 4
    names(out_wt) <- "WT"

    out_both <- c(3, 1)
    names(out_both) <- c("WT", "aei")

    expect_equal(data_to_counts(m00, "vector"), out_wt)
    expect_equal(data_to_counts(m0_2c, "vector"), out_wt)
    expect_equal(data_to_counts(m0_1c, "vector"), out_wt)
    expect_equal(data_to_counts(m0_2c_named, "vector"), out_wt)
    expect_equal(data_to_counts(m0_1c_named, "vector"), out_wt)    

    expect_equal(data_to_counts(m1_nammed, "vector"), out_both)
    expect_error(data_to_counts(m1_unnamed, "vector"),
                 "!is.null(colnames(data)) is not TRUE", fixed = TRUE)
})


cat("\n Done test.data_to_counts.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
