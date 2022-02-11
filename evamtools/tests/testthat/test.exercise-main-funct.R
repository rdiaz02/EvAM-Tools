## This test fails if you try to use MCCBN
test_that("Minimal test: we can run", {
    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:6]
    out <- suppressMessages(evam(Dat1,
                                 methods = c("CBN", "OT", "OncoBN",
                                             "MHN", "HESBCN")))
    expect_true(exists("OT_model", where = out))
})


exercise_sample_all_CPMs <- function(out) {
    samp <- evamtools:::sample_all_CPMs(out, 1000)
    se <- paste0(c("CBN", "OT", "OncoBN", "MHN", "HESBCN"),
                 "_genotype_freqs")
    expect_true(all(vapply(se, function(x) exists(x, samp), TRUE)))  
}


test_that("We can deal with duplicated columns and columns without events and constant columns", {

    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:6]
    Dat1[, "rep_1"] <- Dat1[, 1]
    Dat1[, "no_event"] <- rep(0, nrow(Dat1))
    Dat1[, "constant"] <- rep(1, nrow(Dat1))
    
    out1 <- suppressMessages(evam(Dat1[, 1:6],
                                  methods = c("CBN", "OT", "OncoBN",
                                              "MHN", "HESBCN")))
    expect_true(exists("OT_model", where = out1))

    out2 <- suppressMessages(evam(Dat1[, c(1:5, 7)],
                                  methods = c("CBN", "OT", "OncoBN",
                                              "MHN", "HESBCN")))
    expect_true(exists("OT_model", where = out2))

    out3 <- suppressMessages(evam(Dat1[, c(1:5, 8)],
                                  methods = c("CBN", "OT", "OncoBN",
                                              "MHN", "HESBCN")))
    expect_true(exists("OT_model", where = out3))


    out4 <- suppressMessages(evam(Dat1,
                                  methods = c("CBN", "OT", "OncoBN",
                                              "MHN", "HESBCN")))
    expect_true(exists("OT_model", where = out4))

    out5 <- suppressMessages(evam(Dat1,
                                  methods = c("CBN", "OT", "OncoBN",
                                              "MHN", "HESBCN"),
                                  max_cols = 3))
    expect_true(exists("OT_model", where = out5))


    exercise_sample_all_CPMs(out1)
    exercise_sample_all_CPMs(out2)
    exercise_sample_all_CPMs(out3)
    exercise_sample_all_CPMs(out4)
    exercise_sample_all_CPMs(out5)

    ## Some examples with random, weird names
    ## Could probably bring to fewer
    iters <- 5
    dd <- Dat1
    dd <- dd[1:30, ]
    for(i in 1:iters) {
        ngenes <- ncol(dd)
        gn <- vector(mode = "character", length = ngenes)

        ## create weird gene names
        for(g in seq_len(ngenes)) {
            l1 <- sample(c(LETTERS, letters), 1)
            rest <-  c(sample(
                c(
                    sample(c("_",  "=", "?", "#", "@", "%", "&", "!"), 4, replace = TRUE)
                    , sample(0:9, 4, replace = TRUE)
                    , sample(c(letters, LETTERS), 4, replace = TRUE)
                )))
            
            rest <- paste(rest, sep = "", collapse = "")
            gn[g] <- paste0(l1, rest)
        }
        colnames(dd) <- gn
        outdd <- suppressMessages(evam(dd,
                                  methods = c("CBN", "OT", "OncoBN",
                                              "MHN", "HESBCN"),
                                  max_cols = 4))
        expect_true(exists("OT_model", where = outdd))
        exercise_sample_all_CPMs(outdd)
    }
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
                                  methods = c("CBN", "OT", "OncoBN",
                                              "MHN", "HESBCN")))
    expect_true(exists("OT_model", where = out3))


    out4 <- suppressMessages(evam(db4,
                                  methods = c("CBN", "OT", "OncoBN",
                                              "MHN", "HESBCN")))
    expect_true(exists("OT_model", where = out4))


    out5 <- suppressMessages(evam(db5,
                                  methods = c("CBN", "OT", "OncoBN",
                                              "MHN", "HESBCN")))
    expect_true(exists("OT_model", where = out5))


    exercise_sample_all_CPMs(out3)
    exercise_sample_all_CPMs(out4)
    exercise_sample_all_CPMs(out5)    
   
})


test_that("We can run evam with non-default arguments", {
    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:8]
    out <- suppressMessages(evam(Dat1,
                                  methods = c("CBN", "OT", "OncoBN",
                                              "MHN", "HESBCN"),
                                 max_cols = 4,
                                 cbn_opts = list(cores = 2),
                                 hesbcn_opts = list(steps = 20000, seed = 2),
                                 mhn_opts = list(lambda = 1/10),
                                 oncobn_opts = list(model = "CBN", epsilon = 0.01),
                                 ot_opts = list(with_errors_dist_ot = FALSE),
                                 mccbn_opts = list(model = "H-CBN2",
                                                   max.iter = 10L,
                                                   L = 6,
                                                   max.iter.asa = 20L)
                                 ))
    expect_true(exists("OT_model", where = out))

    exercise_sample_all_CPMs(out)
})

cat("\n Done test.exercise-main-funct.R \n")
