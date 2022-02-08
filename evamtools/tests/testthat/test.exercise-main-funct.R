## This test fails if you try to use MCCBN
test_that("Minimal test: we can run", {
    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:6]
    out <- suppressMessages(evam(Dat1,
                                 methods = c("CBN", "OT",
                                             "MHN", "HESBCN",
                                             "DBN")))
    expect_true(exists("OT_model", where = out))
})




test_that("We can deal with duplicated columns and columns without events and constant columns", {

    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:6]
    Dat1[, "rep_5"] <- Dat1[, 5]
    Dat1[, "no_event"] <- rep(0, nrow(Dat1))
    Dat1[, "constant"] <- rep(1, nrow(Dat1))
    
    out1 <- suppressMessages(evam(Dat1[, 1:6],
                                  methods = c("CBN", "OT",
                                              "MHN", "HESBCN",
                                              "DBN")))
    expect_true(exists("OT_model", where = out1))

    out2 <- suppressMessages(evam(Dat1[, c(1:5, 7)],
                                  methods = c("CBN", "OT",
                                              "MHN", "HESBCN",
                                              "DBN")))
    expect_true(exists("OT_model", where = out2))

    out3 <- suppressMessages(evam(Dat1[, c(1:5, 8)],
                                  methods = c("CBN", "OT",
                                              "MHN", "HESBCN",
                                              "DBN")))
    expect_true(exists("OT_model", where = out3))


    out4 <- suppressMessages(evam(Dat1,
                                  methods = c("CBN", "OT",
                                              "MHN", "HESBCN",
                                              "DBN")))
    expect_true(exists("OT_model", where = out4))

})


test_that("Examples from initial-simple-examples", {
    ## Very weird behaviour with CBN
    ## very sensitive to the n00
    N <- 100
    na <- N + round( 10 * runif(1))
    nc <- N  + round( 10 * runif(1))
    nac <- .5 * N  + round( 10 * runif(1))
    nabc <- .2 * N + round( 10 * runif(1))
    n00 <- 19 ## 1.95 * N  ## 1000 + round( 10 * runif(1))
    dB <- matrix(
        c(
            rep(c(1, 0, 0), na) 
          , rep(c(0, 0, 1), nc)
          , rep(c(1, 0, 1), nac)
          , rep(c(1, 1, 1), nabc)
          , rep(c(0, 0, 0), n00)
        ), ncol = 3, byrow = TRUE
    )
    colnames(dB) <- LETTERS[1:3]
    db3 <- evamtools:::remove_WT(dB, 1)
    db4 <- evamtools:::add_WT(db3, 10 * nrow(db3))
    db5 <- evamtools:::add_WT(db3, 5 * nrow(db3))

    out3 <- suppressMessages(evam(db3,
                                  methods = c("CBN", "OT",
                                              "MHN", "HESBCN",
                                              "DBN")))
    expect_true(exists("OT_model", where = out3))


    out4 <- suppressMessages(evam(db4,
                                  methods = c("CBN", "OT",
                                              "MHN", "HESBCN",
                                              "DBN")))
    expect_true(exists("OT_model", where = out4))


    out5 <- suppressMessages(evam(db5,
                                  methods = c("CBN", "OT",
                                              "MHN", "HESBCN",
                                              "DBN")))
    expect_true(exists("OT_model", where = out5))
})



cat("\n Done test.exercise-main-funct.R \n")
