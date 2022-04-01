t1 <- Sys.time()

test_that("Minimal test: we can run", {
    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:6]
    out <- suppressMessages(evam(Dat1,
                                 methods = c("CBN", "OT", "OncoBN",
                                             "MHN", "HESBCN", "MCCBN")))
    expect_true(exists("OT_model", where = out))

    out2 <- evam(Dat1[, 1:3],
                 methods = c("OT", "OncoBN",
                             "MHN", "CBN"),
                 paths_max = TRUE)
    expect_true(is.na(out2$HESBCN_paths_max))
    expect_true(all(!is.na(out2$CBN_paths_max)))
})


exercise_sample_CPMs <- function(out) {
    samp <- sample_CPMs(out, 1000,
                                    output = c("sampled_genotype_counts",
                                               "obs_genotype_transitions"))
    samp2 <- sample_CPMs(out, 1000, output = "sampled_genotype_counts")
    se <- paste0(c("CBN", "OT", "OncoBN", "MHN", "HESBCN"),
                 "_sampled_genotype_counts")
    expect_true(all(vapply(se, function(x) exists(x, samp), TRUE)))
    expect_true(exists("CBN_obs_genotype_transitions", samp))
    expect_true(all(vapply(se, function(x) exists(x, samp2), TRUE)))
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


    exercise_sample_CPMs(out1)
    exercise_sample_CPMs(out2)
    exercise_sample_CPMs(out3)
    exercise_sample_CPMs(out4)
    exercise_sample_CPMs(out5)

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
        exercise_sample_CPMs(outdd)
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
    db3 <- remove_WT(dB, 1)
    db4 <- add_WT(db3, 10 * nrow(db3))
    db5 <- add_WT(db3, 5 * nrow(db3))

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


    exercise_sample_CPMs(out3)
    exercise_sample_CPMs(out4)
    exercise_sample_CPMs(out5)    
   
})


test_that("We can run evam with non-default arguments", {
    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:8]
    out <- suppressMessages(evam(Dat1,
                                  methods = c("CBN", "OT", "OncoBN",
                                              "MHN", "HESBCN", "MCCBN"),
                                 max_cols = 4,
                                 cores = 1,
                                 hesbcn_opts = list(steps = 20000,
                                                    seed = 2,
                                                    reg = "aic",
                                                    silent = FALSE),
                                 mhn_opts = list(lambda = 1/10),
                                 oncobn_opts = list(model = "CBN", epsilon = 0.01),
                                 ot_opts = list(with_errors_dist_ot = FALSE),
                                 mccbn_opts = list(model = "H-CBN2",
                                                   max.iter = 10L,
                                                   L = 6,
                                                   max.iter.asa = 20L)
                                 ))
    expect_true(exists("OT_model", where = out))
    exercise_sample_CPMs(out)

    out2 <- suppressMessages(evam(Dat1,
                                  methods = c("CBN", "OT", "OncoBN",
                                              "MHN", "HESBCN", "MCCBN"),
                                 max_cols = 4,
                                 cbn_opts = list(omp_threads = 2),
                                 hesbcn_opts = list(steps = 20000, seed = 2),
                                 mhn_opts = list(lambda = 1/10),
                                 oncobn_opts = list(model = "CBN", epsilon = 0.01),
                                 ot_opts = list(with_errors_dist_ot = FALSE),
                                 mccbn_opts = list(model = "OT-CBN",
                                                   max.iter = 10L,
                                                   L = 6,
                                                   max.iter.asa = 20L)
                                 ))
    expect_true(exists("OT_model", where = out2))
    exercise_sample_CPMs(out2)


    expect_warning(out3 <- suppressMessages(evam(Dat1,
                                  methods = c("CBN", "OT", "OncoBN",
                                              "MHN", "HESBCN", "MCCBN"),
                                  max_cols = 4,
                                  cores = 1, ## if >, will hang because of MCCBN thrds > 1
                                 cbn_opts = list(omp_threads = 2, cucu = 9),
                                 hesbcn_opts = list(steps = 20000, seed = 2),
                                 mhn_opts = list(lambda = 1/10),
                                 oncobn_opts = list(model = "CBN", epsilon = 0.01),
                                 ot_opts = list(with_errors_dist_ot = FALSE),
                                 mccbn_opts = list(model = "H-CBN2",
                                                      thrds = 8, max.iter = 20,
                                                   max.iter.asa = 25))),
                   "Option(s) cucu",
                   fixed = TRUE)

    
    expect_true(exists("OT_model", where = out3))
    exercise_sample_CPMs(out3)
   
})


test_that("Handling invalid methods and no single valid method", {
    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:25, 2:8]
    expect_warning(out3 <- suppressMessages(evam(Dat1,
                                  methods = c("OT", "OncoBN", "cucu", "OT"),
                                  max_cols = 4,
                                  cores = 1, ## if >, will hang because of MCCBN thrds > 1
                                 hesbcn_opts = list(steps = 20000, seed = 2),
                                 mhn_opts = list(lambda = 1/10),
                                 oncobn_opts = list(model = "CBN", epsilon = 0.01),
                                 ot_opts = list(with_errors_dist_ot = FALSE),
                                 mccbn_opts = list(model = "H-CBN2",
                                                      thrds = 8, max.iter = 20,
                                                   max.iter.asa = 25))),
                   "Method(s) cucu not among",
                   fixed = TRUE)
    expect_true(exists("OT_model", where = out3))

    sample_CPMs(out3, 1000)
    
    expect_error(out4 <- suppressWarnings(evam(Dat1,
                                  methods = c("coco", "cucu"),
                                  max_cols = 4,
                                  cores = 1, ## if >, will hang because of MCCBN thrds > 1
                                 hesbcn_opts = list(steps = 20000, seed = 2),
                                 mhn_opts = list(lambda = 1/10),
                                 oncobn_opts = list(model = "CBN", epsilon = 0.01),
                                 ot_opts = list(with_errors_dist_ot = FALSE),
                                 mccbn_opts = list(model = "H-CBN2",
                                                      thrds = 8, max.iter = 20,
                                                   max.iter.asa = 25))),
                   "No valid methods given",
                   fixed = TRUE)
    

})


test_that("Options are what they should when we change them", {
    dB_c1 <- matrix(
 c(
     rep(c(1, 0, 0, 0, 0), 30) #A
   , rep(c(0, 0, 1, 0, 0), 30) #C
   , rep(c(1, 1, 0, 0, 0), 20) #AB
   , rep(c(0, 0, 1, 1, 0), 20) #CD
   , rep(c(1, 1, 1, 0, 0), 10) #ABC
   , rep(c(1, 0, 1, 1, 0), 10) #ACD
   , rep(c(1, 1, 0, 0, 1), 10) #ABE
   , rep(c(0, 0, 1, 1, 1), 10) #CDE
   , rep(c(1, 1, 1, 0, 1), 10) #ABCE
   , rep(c(1, 0, 1, 1, 1), 10) #ACDE
   , rep(c(1, 1, 1, 1, 0), 5) # ABCD
   , rep(c(0, 0, 0, 0, 0), 1) # WT
 ), ncol = 5, byrow = TRUE
)
    colnames(dB_c1) <- LETTERS[1:5]

out3 <- suppressMessages(evam(dB_c1,
                                  methods = c("OT", "OncoBN",
                                              "MHN"),
                                  max_cols = 4,
                                  cores = 1, ## if >, will hang because of MCCBN thrds > 1
                                 cbn_opts = list(omp_threads = 2),
                                 hesbcn_opts = list(steps = 20000, seed = 2),
                                 mhn_opts = list(lambda = 100/nrow(dB_c1)),
                              oncobn_opts = list(model = "CBN",
                                                 epsilon = max(colMeans(dB_c1)),
                                                 silent = FALSE),
                                 ot_opts = list(with_errors_dist_ot = FALSE),
                                 mccbn_opts = list(model = "H-CBN2",
                                                   thrds = 8, max.iter = 20,
                                                   max.iter.asa = 25,
                                                   silent = FALSE,
                                                   tmp_dir = "tmpdmccbn")))
    
    expect_true(out3$all_options$mhn_opts$lambda == 100/nrow(dB_c1))
    expect_true(out3$all_options$oncobn_opts$epsilon == max(colMeans(dB_c1)))
})


test_that("Exercise some internal code, inaccessible o.w." ,{
    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:6]
    o1 <- ot_proc(Dat1, nboot = 3, distribution.oncotree = FALSE)
    expect_true(!is.na(o1$ot.boot.original))
    expect_true(is.na(o1$predicted_genotype_freqs))

    o2 <- cbn_proc(Dat1[, 1:3],
                   addname = "tmpi",
                   parall = TRUE,
                   mc.cores = 2,
                   omp_threads = NULL,
                   nboot = 2)
    expect_true(all(!is.na(o2$edges)))

    o3 <- cbn_proc(Dat1[, 1:3],
                   addname = NULL,
                   init.poset = "linear",
                   parall = FALSE,
                   mc.cores = 2,
                   verbose = TRUE,
                   nboot = 2)
    expect_true(all(!is.na(o3$edges)))
    o4 <- cbn_proc(Dat1[, 1:3],
                               addname = NULL,
                               init.poset = "linear",
                               parall = FALSE,
                               verbose = TRUE,
                               silent = FALSE,
                               mc.cores = 2,
                               nboot = 0)
    expect_true(all(is.na(o4$edges$CBN_edgeBootFreq)))

    expect_error(df_2_mat_integer(matrix(c(1, 0, 1e-9, 1), ncol = 2)),
                 "Not in 0L, 1L", fixed = TRUE)
    expect_error(df_2_mat_integer(matrix(c(1, 0, "1", 1), ncol = 2)))
    expect_error(df_2_mat_integer(matrix(c(1, 0, "a", 1), ncol = 2)))
    expect_error(df_2_mat_integer(data.frame(a = c(1, 0), b = c("1", 0))))
    expect_error(df_2_mat_integer(data.frame(a = c(1, 0), b = c(1e-10, 0))))

    expect_error(pre_process(matrix(1:4, ncol = 2)),
                 "Values in x not in 0L, 1L", fixed = TRUE)
    expect_silent(null1 <- pre_process(cbind(c1 = c(1L, 1L), c2 = c(0L, 1L)),
                                       remove.constant = TRUE))
    expect_identical(null1,
                 matrix(c(0L, 1L), dimnames = list(NULL, c("c2")), ncol = 1))
    expect_silent(null2 <- pre_process(cbind(c1 = rep(c(0L, 1L), 50),
                                             c2 = c(rep(0L, 97), rep(1L, 3))),
                                       remove.constant = TRUE,
                                       min.freq = 0.04))
    expect_identical(null2,
                     matrix(rep(c(0L, 1L), 50),
                            dimnames = list(NULL, c("c1")), ncol = 1))
    expect_error(null2 <- pre_process(cbind(c1 = rep(c(0L, 1L), 50),
                                             c2 = c(rep(0L, 97), rep(1L, 3))),
                                       remove.constant = TRUE,
                                      min.freq = -1e-9),
                 "min.freq has to be positive or 0", fixed = TRUE)

    m1 <- matrix(1:9, ncol = 3, dimnames = list(c("23", "Root", "17"),
                                                c("23", "Root", "17")))
    expect_silent(m1s <- sortAdjMat(m1))
    m1se <- matrix(c(5L, 6L, 4L,
                     8L, 9L, 7L,
                     2L, 3L, 1L),
                   ncol = 3,
                   dimnames = list(c("Root", "17", "23"),
                                   c("Root", "17", "23")))
    expect_identical(m1s, m1se)
    
    
    ## Fail graciously
    expect_true(cpm2tm(NULL)$fgraph == "ERROR_CPM_ANALYSIS")
    expect_true(cpm2tm(NA)$fgraph == "ERROR_CPM_ANALYSIS")
    expect_true(
        cpm2tm(try(log("a"), silent = TRUE))$fgraph == "ERROR_CPM_ANALYSIS")

    dfe <- list(edges = data.frame(From = c("Root", "A"),
                                   To = c("A", "B"),
                                   lambda = 3,
                                   OT_edgeWeight = 2))

    expect_error(cpm2tm(dfe),
                 "more than one column with weights",
                 fixed = TRUE)

    dfe2 <-  list(edges = data.frame(From = c("A"),
                                     To = c("B"),
                                     lambda = 3,
                                     OT_edgeWeight = 2))
    expect_error(cpm2tm(dfe2),
                 "Some error here",
                 fixed = TRUE)

    ## The error is  "using 'as.environment(NULL)'
    ## because we pass a NULL object
    ## The actual error is irrelevant. What matters is we fail
    expect_error(cpm2tm(dfe$edges))
    
    dfe3 <- list(edges = data.frame(From = c("Root", "A"),
                                   To = c("A", "B"),
                                   lambda = 3),
                 parent_set = c(A = "Single", B = "Single"))

    expect_error(cpm2tm(dfe3),
                 "x has parent_set but no Relation",
                 fixed = TRUE)
})


test_that("Miscell error conditions", {
   data(every_which_way_data)
   Dat1 <- every_which_way_data[[16]][1:40, 2:6]
   d1e <- Dat1
   colnames(d1e)[1] <- "cucu,tras"
   expect_error(evam(d1e),
                "At least one of your gene names has a comma",
                fixed = TRUE)
   colnames(d1e)[1] <- "cucu\tras"
   expect_error(evam(d1e),
                "At least one of your gene names has a space",
                fixed = TRUE)
   colnames(d1e)[1] <- "cucu\\tras"
   expect_error(evam(d1e),
                "At least one of your gene names has a backslash",
                fixed = TRUE)
   colnames(d1e)[1] <- "cucu tras"
   expect_error(evam(d1e),
                "At least one of your gene names has a space",
                fixed = TRUE)
   colnames(d1e)[1] <- "WT"
   expect_error(evam(d1e),
                "One of your genes is called WT",
                fixed = TRUE)
   expect_error(evam(Dat1[, 1, drop = FALSE]),
                "Fewer than 2 columns",
                fixed = TRUE)
   expect_error(evam(Dat1[, 1]),
                "colnames must exist",
                fixed = TRUE)
   
 expect_error(evam(Dat1[, 1:3], methods = "CBN",
                      cbn_opts = list(init_poset = "custom")),
                 "CBN's init_poset must be one of OT or linear",
                 fixed = TRUE)
})


cat("\n Done test.exercise-main-funct.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
