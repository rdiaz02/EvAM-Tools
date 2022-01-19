## We will use this hazard matrix and transition rates several times below.
create_MHN_test_data <- function() {
    ## from rounding the first 4 x 4 entries of Schill's
    ## cancer example

    me1 <- structure(c(
        0.5, 0, 0, 0,
        0, -1.2, 0, 0,
        0, 0, -1.1, 0,
        -0.1, 0, 0.3, 0.1
    ),
    .Dim = c(4L, 4L),
    .Dimnames =
        list(
            c("A", "B", "C", "D"),
            c("A", "B", "C", "D")
        )
    )

    eme1 <- exp(me1)

    ## off-diagonal terms first
    trm0 <- matrix(0, ncol = 16, nrow = 16)

    geneNames <- colnames(eme1)
    k <- 4
    genots <- c(list(rep(0, k)), allGenotypes_former(k))
    genotNames <- unlist(
        lapply(
            genots,
            function(x) {
                  paste(geneNames[which(x == 1L)], sep = "", collapse = ", ")
              }
        )
    )
    genotNames[genotNames == ""] <- "WT"

    colnames(trm0) <- rownames(trm0) <- genotNames

    trm0["WT", "A"] <- eme1["A", "A"]
    trm0["WT", "B"] <- eme1["B", "B"]
    trm0["WT", "C"] <- eme1["C", "C"]
    trm0["WT", "D"] <- eme1["D", "D"]

    trm0["A", "A, B"] <- eme1["B", "B"] * eme1["B", "A"]
    trm0["A", "A, C"] <- eme1["C", "C"] * eme1["C", "A"]
    trm0["A", "A, D"] <- eme1["D", "D"] * eme1["D", "A"]

    trm0["B", "A, B"] <- eme1["A", "A"] * eme1["A", "B"]
    trm0["B", "B, C"] <- eme1["C", "C"] * eme1["C", "B"]
    trm0["B", "B, D"] <- eme1["D", "D"] * eme1["D", "B"]

    trm0["C", "A, C"] <- eme1["A", "A"] * eme1["A", "C"]
    trm0["C", "B, C"] <- eme1["B", "B"] * eme1["B", "C"]
    trm0["C", "C, D"] <- eme1["D", "D"] * eme1["D", "C"]

    trm0["D", "A, D"] <- eme1["A", "A"] * eme1["A", "D"]
    trm0["D", "B, D"] <- eme1["B", "B"] * eme1["B", "D"]
    trm0["D", "C, D"] <- eme1["C", "C"] * eme1["C", "D"]

    trm0["A, B", "A, B, C"] <- eme1["C", "C"] * eme1["C", "A"] * eme1["C", "B"]
    trm0["A, B", "A, B, D"] <- eme1["D", "D"] * eme1["D", "A"] * eme1["D", "B"]

    trm0["A, C", "A, B, C"] <- eme1["B", "B"] * eme1["B", "A"] * eme1["B", "C"]
    trm0["A, C", "A, C, D"] <- eme1["D", "D"] * eme1["D", "A"] * eme1["D", "C"]

    trm0["A, D", "A, B, D"] <- eme1["B", "B"] * eme1["B", "A"] * eme1["B", "D"]
    trm0["A, D", "A, C, D"] <- eme1["C", "C"] * eme1["C", "A"] * eme1["C", "D"]

    trm0["B, C", "A, B, C"] <- eme1["A", "A"] * eme1["A", "B"] * eme1["A", "C"]
    trm0["B, C", "B, C, D"] <- eme1["D", "D"] * eme1["D", "B"] * eme1["D", "C"]

    trm0["B, D", "A, B, D"] <- eme1["A", "A"] * eme1["A", "B"] * eme1["A", "D"]
    trm0["B, D", "B, C, D"] <- eme1["C", "C"] * eme1["C", "B"] * eme1["C", "D"]

    trm0["C, D", "A, C, D"] <- eme1["A", "A"] * eme1["A", "C"] * eme1["A", "D"]
    trm0["C, D", "B, C, D"] <- eme1["B", "B"] * eme1["B", "C"] * eme1["B", "D"]

    trm0["A, B, C", "A, B, C, D"] <-
        eme1["D", "D"] * eme1["D", "A"] * eme1["D", "B"] * eme1["D", "C"]

    trm0["A, B, D", "A, B, C, D"] <-
        eme1["C", "C"] * eme1["C", "A"] * eme1["C", "B"] * eme1["C", "D"]

    trm0["A, C, D", "A, B, C, D"] <-
        eme1["B", "B"] * eme1["B", "A"] * eme1["B", "C"] * eme1["B", "D"]

    trm0["B, C, D", "A, B, C, D"] <-
        eme1["A", "A"] * eme1["A", "B"] * eme1["A", "C"] * eme1["A", "D"]

    trm1 <- trm0
    diag(trm1) <- -1 * rowSums(trm0)

    return(list(
        me1 = me1,
        eme1 = eme1,
        trm1 = trm1,
        trm0 = trm0
    ))
}


test_that("Diagonals of MHN trm", {
    ## paranoid check
    ## Adding the diagonal
    ## Q_ii = - \Sum Q_ji
    suppressWarnings(try(rm(tmp, trm1, trm0, me1, eme1), silent = TRUE))
    tmp <- create_MHN_test_data()
    expect_equal(rep(0, 16), rowSums(tmp$trm1),
        check.attributes = FALSE
        )
    rm(tmp)
})


test_that("transition rates between different genotypes are the
  same from time discretized and using the competing exponentials
  approach", {
      suppressWarnings(try(rm(tmp, trm1, trm0, me1, eme1), silent = TRUE))
      tmp <- create_MHN_test_data()
      trm1 <- tmp$trm1
      trm0 <- tmp$trm0
      rm(tmp)

      ## Uniformization
      gamma <- max(abs(diag(trm1)))
      time_discr_trm <- diag(16) + trm1 / gamma
      
      trans_mat_time_discr <- time_discr_trm
      diag(trans_mat_time_discr) <- 0
      s_tmtd <- rowSums(trans_mat_time_discr)
      trans_mat_time_discr <- sweep(trans_mat_time_discr, 1, s_tmtd, "/")
      ## note: last row is NaN
      trans_mat_time_discr[16, ] <- 0
      ## paranoid check
      expect_equal(c(rep(1, 15), 0),
                   rowSums(trans_mat_time_discr),
                   check.attributes = FALSE
                   )

      ## Competing exponentials
      trans_mat_comp_exp <- trm0
      s_tmce <- rowSums(trans_mat_comp_exp)
      trans_mat_comp_exp <- sweep(trans_mat_comp_exp, 1, s_tmce, "/")
      trans_mat_comp_exp[16, ] <- 0 ## last row is 0
      ## paranoid check
      expect_equal(c(rep(1, 15), 0),
                   rowSums(trans_mat_comp_exp),
                   check.attributes = FALSE
                   )

      ## Verify identity, as it should be
      expect_equal(
          trans_mat_comp_exp,
          trans_mat_time_discr
      )

      ## For us, now, it is faster and simpler if we just use the competing
      ## exponentials one.

      tm_u <- trans_rate_to_trans_mat(trm0, method = "uniformization")
      tm_ce <- trans_rate_to_trans_mat(trm0, method = "competingExponentials")

      expect_equal(tm_u[nrow(tm_u), ncol(tm_u)], 1)
      expect_equal(tm_u, time_discr_trm)
      expect_equal(tm_ce, trans_mat_comp_exp)
})



test_that("Verify code to obtain transition rate matrix, against
a reference and using different algorithms", {
    suppressWarnings(try(rm(tmp, trm1, trm0, me1, eme1), silent = TRUE))
    tmp <- create_MHN_test_data()
    me1 <- tmp$me1
    eme1 <- tmp$eme1
    trm0 <- tmp$trm0
    rm(tmp)

    ttr_1 <- theta_to_trans_rate_1(me1)
    ttr_2 <- theta_to_trans_rate(me1, inner_transition = inner_transitionRate_1)
    ttr_22 <- theta_to_trans_rate(me1,
                                  inner_transition = inner_transitionRate_2)
    ttr_31 <- theta_to_trans_rate_3(me1,
                                    inner_transition = inner_transitionRate_3_1)
    ttr_32 <- theta_to_trans_rate_3(me1,
                                    inner_transition = inner_transitionRate_3_2)

    expect_identical(ttr_1, ttr_2)
    expect_identical(ttr_2, ttr_22)
    expect_identical(ttr_2, ttr_31)
    expect_identical(ttr_2, ttr_32)

    expect_identical(ttr_1, trm0)
    expect_identical(ttr_2, trm0)
    expect_identical(ttr_31, trm0)
    expect_identical(ttr_32, trm0)

    ## In particular, note the asymmetrical changes
    ## Figure 2, right
    ## From D mutated to both C and D mutated
    ## theta_33 * theta_34
    expect_equal(eme1[3, 3] * eme1[3, 4], trm0["D", "C, D"])
    ## This is different, of course
    expect_false(eme1[3, 3] * eme1[4, 3] == trm0["D", "C, D"])
})



test_that("Compute TRM gives the same results with sparse matrices", {
    theta0 <- structure(c(
        0.26, -0.57, 0, -0.97, -0.43, -1.18, 0, -0.29, 0,
        0, -1.84, 0, -0.23, -0.5, -0.04, -0.24
    ), .Dim = c(4L, 4L), .Dimnames = list(
        c("A", "B", "C", "D"), c("A", "B", "C", "D")
    ))

    trm0 <- theta_to_trans_rate_3(theta0,
        inner_transition = inner_transitionRate_3_1
        )

    trmSM <-
        theta_to_trans_rate_3_SM(theta0,
                                 inner_transition = inner_transitionRate_3_1
                                 )

    trmSMm <- as.matrix(trmSM)
    expect_identical(trm0, trmSMm)
})
    

test_that("we are using the indices of theta correctly 1", {

    ## Create some data
    N <- 200
    na <- N
    nc <- N + 21
    nab <- N + 5
    nbc <- N + 7 + round(10 * runif(1))
    n00 <- N / 10 + round(10 * runif(1))
    dB <- matrix(
        c(
            rep(c(1, 0, 0, 0), na)
          , rep(c(0, 0, 1, 0), nc)
          , rep(c(1, 1, 0, 0), nab)
          , rep(c(0, 1, 1, 0), nbc)
          , rep(c(0, 0, 0, 0), n00)
        ), ncol = 4, byrow = TRUE
    )
    colnames(dB) <- LETTERS[1:4]
    sampledGenotypes(dB) ## from OncoSimulR
    ## Note :the next crashes. Of course: don't do silly things like passing
    ## a data frame with one or more colSums == 0
    try(do_MHN(dB), silent = TRUE)
    dB2 <- dB[, -4]
    mm1 <- do_MHN(dB2)
    ## From A to A, B:
    expect_equal(mm1$transitionRateMatrix["A", "A, B"],
                 exp(mm1$theta["B", "B"]) * exp(mm1$theta["B", "A"]))

    ## The next are left as interesting examples of why checking
    ## for inequality does not always work.
    ## If symmetry and thetas are identical A, B and B, A, if na = nab. That
    ## is why we make sure different N above
    ## stopifnot(mm1$transitionRateMatrix["A", "A, B"] != 
    ##           exp(mm1$theta["B", "B"]) * exp(mm1$theta["A", "B"]))
    
    ## stopifnot(mm1$transitionRateMatrix["A", "A, B"] != 
    ##           exp(mm1$theta["A", "A"]) * exp(mm1$theta["A", "B"]))
    
    ## stopifnot(mm1$transitionRateMatrix["A", "A, B"] != 
    ##           exp(mm1$theta["A", "A"]) * exp(mm1$theta["B", "A"]))
    
    ## if(isTRUE(all.equal(
    ##   mm1$transitionRateMatrix["A", "A, B"], 
    ##   exp(mm1$theta["B", "B"]) * exp(mm1$theta["A", "B"])))) {
    ##   cat("\n \n Here they are identical to 0_1")
    ##   print(exp(mm1$theta["A", "B"]))
    ##   print(exp(mm1$theta["B", "A"]))
    ##   print(exp(mm1$theta["B", "B"]))
    ##   print(sampledGenotypes(dB))
    ## }
    
    ## if(isTRUE(all.equal(
    ##   mm1$transitionRateMatrix["A", "A, B"], 
    ##   exp(mm1$theta["A", "A"]) * exp(mm1$theta["A", "B"])))) {
    ##   cat("\n \n Here they are identical to 0_2")
    ##   print(exp(mm1$theta["A", "B"]))
    ##   print(exp(mm1$theta["B", "A"]))
    ##   print(exp(mm1$theta["A", "A"]))
    ##   print(sampledGenotypes(dB))
    ## }
    
    
    ## if(isTRUE(all.equal(
    ##   mm1$transitionRateMatrix["A", "A, B"], 
    ##   exp(mm1$theta["A", "A"]) * exp(mm1$theta["B", "A"])))) {
    ##   cat("\n \n Here they are identical to 0_3")
    ##   print(exp(mm1$theta["A", "B"]))
    ##   print(exp(mm1$theta["B", "A"]))
    ##   print(exp(mm1$theta["A", "A"]))
    ##   print(sampledGenotypes(dB))
    ## }
  

})


test_that("we are using the indices of theta correctly 2", {

    N <- 200
    na <- N
    nc <- N + 3 + round(10 * runif(1))
    nab <- N + 5 + round(10 * runif(1))
    ncd <- N + 7 + round(10 * runif(1))
    n00 <- N / 10 + round(10 * runif(1))
    dB <- matrix(
        c(
            rep(c(1, 0, 0, 0), na) 
          , rep(c(0, 0, 1, 0), nc)
          , rep(c(1, 1, 0, 0), nab)
          , rep(c(0, 0, 1, 1), ncd)
          , rep(c(0, 0, 0, 0), n00)
        ), ncol = 4, byrow = TRUE
    )
    colnames(dB) <- LETTERS[1:4]
    sampledGenotypes(dB)
    mm4 <- do_MHN(dB)

    round(mm4$transitionMatrixCompExp, 3)

    ## From A to A, B:
    expect_equal(mm4$transitionRateMatrix["A", "A, B"],
                 exp(mm4$theta["B", "B"]) * exp(mm4$theta["B", "A"]))
    
    expect_equal( mm4$transitionRateMatrix["B", "A, B"],
                 exp(mm4$theta["A", "A"]) * exp(mm4$theta["A", "B"]))
    
    
    expect_equal( mm4$transitionRateMatrix["C", "C, D"],
                 exp(mm4$theta["D", "D"]) * exp(mm4$theta["D", "C"]))
    
    
    expect_equal( mm4$transitionRateMatrix["D", "C, D"],
                 exp(mm4$theta["C", "C"]) * exp(mm4$theta["C", "D"]))
    
    
    expect_equal( mm4$transitionRateMatrix["D", "B, D"],
                 exp(mm4$theta["B", "B"]) * exp(mm4$theta["B", "D"]))
    
    
    expect_equal( mm4$transitionRateMatrix["C, D", "A, C, D"],
                 exp(mm4$theta["A", "A"]) *
                 exp(mm4$theta["A", "D"]) *
                 exp(mm4$theta["A", "C"]))
    
    
    expect_equal( mm4$transitionRateMatrix["B, C", "A, B, C"],
                 exp(mm4$theta["A", "A"]) *
                 exp(mm4$theta["A", "B"]) *
                 exp(mm4$theta["A", "C"]))
    
    expect_equal( mm4$transitionRateMatrix["B, C", "B, C, D"],
                 exp(mm4$theta["D", "D"]) *
                 exp(mm4$theta["D", "B"]) *
                 exp(mm4$theta["D", "C"]))
    
    
    expect_equal( mm4$transitionRateMatrix["B, C, D", "A, B, C, D"],
                 exp(mm4$theta["A", "A"]) *
                 exp(mm4$theta["A", "B"]) *
                 exp(mm4$theta["A", "D"]) *
                 exp(mm4$theta["A", "C"]))
    
    
    ## A qualitative assessment: output matches the observed
    ## From A we move mostly to A, B
    ## From B, we can move to either A, B or B, C
    ## From C, we move to C, D
    
})


test_that("we are using the indices of theta correctly 3", {
    N <- 200
    na <- N
    nc <- 5 * N + 13 + round( 10 * runif(1))
    nd <- round(1.5 * N) + round( 10 * runif(1))
    nab <- 2 * N + round( 10 * runif(1))
    nabc <- 3 * N + round(10 * runif(1))
    nabd <- 5 + round(10 * runif(1))    
    n00 <- N/10 + round( 10 * runif(1))
    dB <- matrix(
        c(
            rep(c(1, 0, 0, 0), na) 
          , rep(c(0, 0, 1, 0), nc)
          , rep(c(0, 0, 0, 1), nd)            
          , rep(c(1, 1, 0, 0), nab)
          , rep(c(1, 1, 1, 0), nabc)
          , rep(c(1, 1, 0, 1), nabd)            
          , rep(c(0, 0, 0, 0), n00)
        ), ncol = 4, byrow = TRUE
    )
    colnames(dB) <- LETTERS[1:4]
    sampledGenotypes(dB)
    mm77 <- do_MHN(dB)
    options(width = 300)
    
    round(mm77$transitionMatrixCompExp, 3)
    round(mm77$transitionMatrixTimeDiscretized, 3)
    
    ## From A to A, B:
    mm77$transitionMatrixCompExp["A", "A, B"]
    mm77$transitionRateMatrix["A", "A, B"]
    exp(mm77$theta["B", "B"]) * exp(mm77$theta["B", "A"])
    

    expect_equal( mm77$transitionRateMatrix["A", "A, B"],
                 exp(mm77$theta["B", "B"]) * exp(mm77$theta["B", "A"]))

    expect_equal( mm77$transitionRateMatrix["A", "A, C"],
                 exp(mm77$theta["C", "C"]) * exp(mm77$theta["C", "A"]))
    
    expect_equal( mm77$transitionRateMatrix["A, B", "A, B, C"],
                 exp(mm77$theta["C", "C"]) *
                 exp(mm77$theta["C", "A"]) *
                 exp(mm77$theta["C", "B"]) 
                 )
    
    expect_equal( mm77$transitionRateMatrix["A, B", "A, B, D"],
                 exp(mm77$theta["D", "D"]) *
                 exp(mm77$theta["D", "A"]) *
                 exp(mm77$theta["D", "B"]) 
                 )
    
    ## A very extreme one, from a non-existing genotype
    expect_equal( mm77$transitionRateMatrix["A, C", "A, B, C"],
                 exp(mm77$theta["B", "B"]) *
                 exp(mm77$theta["B", "A"]) *
                 exp(mm77$theta["B", "C"]) 
                 )
    
    expect_equal( mm77$transitionRateMatrix["B, C, D", "A, B, C, D"],
                 exp(mm77$theta["A", "A"]) *
                 exp(mm77$theta["A", "B"]) *
                 exp(mm77$theta["A", "D"]) *
                 exp(mm77$theta["A", "C"]))
    
    expect_equal( mm77$transitionRateMatrix["A, B, C", "A, B, C, D"],
                 exp(mm77$theta["D", "D"]) *
                 exp(mm77$theta["D", "B"]) *
                 exp(mm77$theta["D", "A"]) *
                 exp(mm77$theta["D", "C"]))
    
    expect_equal( mm77$transitionRateMatrix["A, B, D", "A, B, C, D"],
                 exp(mm77$theta["C", "C"]) *
                 exp(mm77$theta["C", "B"]) *
                 exp(mm77$theta["C", "D"]) *
                 exp(mm77$theta["C", "A"]))
    
    expect_equal( mm77$transitionRateMatrix["A, C, D", "A, B, C, D"],
                 exp(mm77$theta["B", "B"]) *
                 exp(mm77$theta["B", "A"]) *
                 exp(mm77$theta["B", "D"]) *
                 exp(mm77$theta["B", "C"]))
  
    
    round(mm77$transitionRateMatrix, 3)
    
    ## And output matches the observed
    ## In particular: A, B to A,B,C vs A, B, D

    ## As above, these could, by chance, be equal
    ## if(isTRUE(all.equal(
    ##   mm77$transitionRateMatrix["A", "A, B"], 
    ##   exp(mm77$theta["B", "B"]) * exp(mm77$theta["A", "B"])))) {
    ##   cat("\n \n Here they are identical too")
    ##   print(exp(mm77$theta["A", "B"]))
    ##   print(exp(mm77$theta["B", "A"]))
    ##   print(exp(mm77$theta["B", "B"]))
    ##   print(sampledGenotypes(dB))
    ## }
    
    ## stopifnot(!isTRUE(all.equal(
    ##     mm77$transitionRateMatrix["A", "A, B"], 
    ##     exp(mm77$theta["A", "A"]) * exp(mm77$theta["A", "B"]))))
    
    ## stopifnot(!isTRUE(all.equal(
    ##     mm77$transitionRateMatrix["A", "A, B"], 
    ##     exp(mm77$theta["A", "A"]) * exp(mm77$theta["B", "A"]))))
})



test_that("we are using the indices of theta correctly 4", {
    ## And a six gene example. Also testing the sparse matrix implementation
    N <- 200
    na <- N
    nc <- 3 * N + 2 + round( 10 * runif(1))
    nd <- round(1.5 * N) + round( 10 * runif(1))
    nab <- 2 * N + round( 10 * runif(1))
    nabc <- 3 * N + round(10 * runif(1))
    nabd <- 5 * N + round(10 * runif(1))
    nabce <- 5 * N + round(10 * runif(1))
    nabcf <- 3 * N + round(10 * runif(1))
    nabcde <- 3 * N + round(10 * runif(1))
    nabcef <- 1 * N + round(10 * runif(1))
    nabcdef <- 2 * N + round(10 * runif(1))
    n00 <- N/10 + round( 10 * runif(1))
    dB <- matrix(
        c(
            rep(c(1, 0, 0, 0, 0, 0), na) 
          , rep(c(0, 0, 1, 0, 0, 0), nc)
          , rep(c(0, 0, 0, 1, 0, 0), nd)            
          , rep(c(1, 1, 0, 0, 0, 0), nab)
          , rep(c(1, 1, 1, 0, 0, 0), nabc)
          , rep(c(1, 1, 0, 1, 0, 0), nabd)            
          , rep(c(0, 0, 0, 0, 0, 0), n00)
          , rep(c(1, 1, 1, 0, 1, 0), nabce)
          , rep(c(1, 1, 1, 0, 0, 1), nabcf)
          , rep(c(1, 1, 1, 1, 1, 0), nabcde)
          , rep(c(1, 1, 1, 0, 1, 1), nabcef)            
          , rep(c(1, 1, 1, 1, 1, 1), nabcdef)            
        ), ncol = 6, byrow = TRUE
    )
    colnames(dB) <- LETTERS[1:6]
    sampledGenotypes(dB)
    mm2 <- do_MHN(dB)
    mm22 <- do_MHN2(dB)
    options(width = 350)
    
    round(mm2$transitionMatrixCompExp, 2)
    round(mm2$transitionMatrixTimeDiscretized, 2)
    
    
    ## Yes, sparse ones
    expect_s4_class(mm22$transitionRateMatrix, "dgCMatrix")
    expect_s4_class(mm22$transitionMatrixTimeDiscretized, "dgCMatrix")
    expect_s4_class(mm22$transitionMatrixCompExp, "dgCMatrix")
    
    ## To compare same values, need to convert into matrices
    expect_identical(mm2$transitionRateMatrix,
                     as.matrix(mm22$transitionRateMatrix))
    ## Identical would fail
    expect_equal(mm2$transitionMatrixCompExp,
                 as.matrix(mm22$transitionMatrixCompExp))
    
    expect_equal(mm2$transitionMatrixTimeDiscretized,
                 as.matrix(mm22$transitionMatrixTimeDiscretized))
    
    
    ## This would fail, unless we turned them into matrices
    expect_false(identical(mm2$transitionRateMatrix,
                           mm22$transitionRateMatrix))
    expect_false(identical(mm2$transitionMatrixTimeDiscretized,
                           mm22$transitionMatrixTimeDiscretized))
    expect_false(identical(mm2$transitionMatrixCompExp,
                           mm22$transitionMatrixCompExp))
        
    ## Incidentally, inference does not seem great for how we go
    ## ABCE to ABCDE and ABCEF. But probably a hard model anyway,
    
    ## From A to A, B:
    mm2$transitionMatrixCompExp["A", "A, B"]
    mm2$transitionRateMatrix["A", "A, B"]
    exp(mm2$theta["B", "B"]) * exp(mm2$theta["B", "A"])
    

    expect_equal( mm2$transitionRateMatrix["A", "A, B"],
                 exp(mm2$theta["B", "B"]) * exp(mm2$theta["B", "A"]))

    expect_equal( mm2$transitionRateMatrix["A", "A, C"],
                        exp(mm2$theta["C", "C"]) * exp(mm2$theta["C", "A"]))
   
    expect_equal( mm2$transitionRateMatrix["A, B", "A, B, C"],
                 exp(mm2$theta["C", "C"]) *
                 exp(mm2$theta["C", "A"]) *
                 exp(mm2$theta["C", "B"]) 
                 )
    
    expect_equal( mm2$transitionRateMatrix["A, B", "A, B, D"],
                 exp(mm2$theta["D", "D"]) *
                 exp(mm2$theta["D", "A"]) *
                 exp(mm2$theta["D", "B"]) 
                 )
    
    ## A very extreme one, from a non-existing genotype
    expect_equal( mm2$transitionRateMatrix["A, C", "A, B, C"],
                 exp(mm2$theta["B", "B"]) *
                 exp(mm2$theta["B", "A"]) *
                 exp(mm2$theta["B", "C"]) 
                 )
    
    
    expect_equal( mm2$transitionRateMatrix["B, C, D", "A, B, C, D"],
                 exp(mm2$theta["A", "A"]) *
                 exp(mm2$theta["A", "B"]) *
                 exp(mm2$theta["A", "D"]) *
                 exp(mm2$theta["A", "C"]))
    
    expect_equal( mm2$transitionRateMatrix["A, B, C", "A, B, C, D"],
                 exp(mm2$theta["D", "D"]) *
                 exp(mm2$theta["D", "B"]) *
                 exp(mm2$theta["D", "A"]) *
                 exp(mm2$theta["D", "C"]))
    
    expect_equal( mm2$transitionRateMatrix["A, B, D", "A, B, C, D"],
                 exp(mm2$theta["C", "C"]) *
                 exp(mm2$theta["C", "B"]) *
                 exp(mm2$theta["C", "D"]) *
                 exp(mm2$theta["C", "A"]))
    
    expect_equal( mm2$transitionRateMatrix["A, C, D", "A, B, C, D"],
                 exp(mm2$theta["B", "B"]) *
                 exp(mm2$theta["B", "A"]) *
                 exp(mm2$theta["B", "D"]) *
                 exp(mm2$theta["B", "C"]))
    
    
    ## 5 genes
    expect_equal( mm2$transitionRateMatrix["A, B, C, E", "A, B, C, D, E"],
                 exp(mm2$theta["D", "D"]) *
                 exp(mm2$theta["D", "A"]) *
                 exp(mm2$theta["D", "B"]) *
                 exp(mm2$theta["D", "C"]) *                        
                 exp(mm2$theta["D", "E"]))
    
    expect_equal( mm2$transitionRateMatrix["A, B, C, E", "A, B, C, E, F"],
                 exp(mm2$theta["F", "F"]) *
                 exp(mm2$theta["F", "A"]) *
                 exp(mm2$theta["F", "B"]) *
                 exp(mm2$theta["F", "C"]) *                        
                 exp(mm2$theta["F", "E"]))
    
    expect_equal( mm2$transitionRateMatrix["A, B, C, F", "A, B, C, E, F"],
                 exp(mm2$theta["E", "E"]) *
                 exp(mm2$theta["E", "A"]) *
                 exp(mm2$theta["E", "B"]) *
                 exp(mm2$theta["E", "C"]) *                        
                 exp(mm2$theta["E", "F"]))
    
    expect_equal( mm2$transitionRateMatrix["A, B, C, F", "A, B, C, D, F"],
                 exp(mm2$theta["D", "D"]) *
                 exp(mm2$theta["D", "A"]) *
                 exp(mm2$theta["D", "B"]) *
                 exp(mm2$theta["D", "C"]) *                        
                 exp(mm2$theta["D", "F"]))
    
    ## six genes
    expect_equal( mm2$transitionRateMatrix["A, B, C, D, E", "A, B, C, D, E, F"],
                 exp(mm2$theta["F", "F"]) *
                 exp(mm2$theta["F", "A"]) *
                 exp(mm2$theta["F", "B"]) *
                 exp(mm2$theta["F", "C"]) *
                 exp(mm2$theta["F", "D"]) *                        
                 exp(mm2$theta["F", "E"]))
    
    expect_equal( mm2$transitionRateMatrix["A, B, C, D, F", "A, B, C, D, E, F"],
                 exp(mm2$theta["E", "E"]) *
                 exp(mm2$theta["E", "A"]) *
                 exp(mm2$theta["E", "B"]) *
                 exp(mm2$theta["E", "C"]) *
                 exp(mm2$theta["E", "D"]) *                        
                 exp(mm2$theta["E", "F"]))
    
    
    expect_equal( mm2$transitionRateMatrix["A, B, C, E, F", "A, B, C, D, E, F"],
                 exp(mm2$theta["D", "D"]) *
                 exp(mm2$theta["D", "A"]) *
                 exp(mm2$theta["D", "B"]) *
                 exp(mm2$theta["D", "C"]) *
                 exp(mm2$theta["D", "E"]) *                        
                 exp(mm2$theta["D", "F"]))
    
    ## six genes, with sparse. Totally unnecessary since they were identical to non-sparse
    expect_equal( mm22$transitionRateMatrix["A, B, C, D, E", "A, B, C, D, E, F"],
                 exp(mm22$theta["F", "F"]) *
                 exp(mm22$theta["F", "A"]) *
                 exp(mm22$theta["F", "B"]) *
                 exp(mm22$theta["F", "C"]) *
                 exp(mm22$theta["F", "D"]) *                        
                 exp(mm22$theta["F", "E"]))
    
    expect_equal( mm22$transitionRateMatrix["A, B, C, D, F", "A, B, C, D, E, F"],
                 exp(mm22$theta["E", "E"]) *
                 exp(mm22$theta["E", "A"]) *
                 exp(mm22$theta["E", "B"]) *
                 exp(mm22$theta["E", "C"]) *
                 exp(mm22$theta["E", "D"]) *                        
                 exp(mm22$theta["E", "F"]))
    
    
    expect_equal( mm22$transitionRateMatrix["A, B, C, E, F", "A, B, C, D, E, F"],
                 exp(mm22$theta["D", "D"]) *
                 exp(mm22$theta["D", "A"]) *
                 exp(mm22$theta["D", "B"]) *
                 exp(mm22$theta["D", "C"]) *
                 exp(mm22$theta["D", "E"]) *                        
                 exp(mm22$theta["D", "F"]))
    
    ## tmp <- round(mm2$transitionRateMatrix, 3)

   ## Same issues as in examples above
 
    ## if(isTRUE(all.equal(
    ##   mm2$transitionRateMatrix["A", "A, B"], 
    ##   exp(mm2$theta["B", "B"]) * exp(mm2$theta["A", "B"])))) {
    ##   cat("\n \n Here they are identical too 3_1")
    ##   print(exp(mm2$theta["A", "B"]))
    ##   print(exp(mm2$theta["B", "A"]))
    ##   print(exp(mm2$theta["B", "B"]))
    ##   print(sampledGenotypes(dB))
    ## }
    ## if(isTRUE(all.equal(
    ##   mm2$transitionRateMatrix["A", "A, B"], 
    ##   exp(mm2$theta["A", "A"]) * exp(mm2$theta["A", "B"])))) {
    ##   cat("\n \n Here they are identical too 3_2")
    ##   print(exp(mm2$theta["A", "B"]))
    ##   print(exp(mm2$theta["B", "A"]))
    ##   print(exp(mm2$theta["A", "A"]))
    ##   print(sampledGenotypes(dB))
    ## }
    ## if(isTRUE(all.equal(
    ##   mm2$transitionRateMatrix["A", "A, B"], 
    ##   exp(mm2$theta["A", "A"]) * exp(mm2$theta["B", "A"])))) {
    ##   cat("\n \n Here they are identical too 3_3")
    ##   print(exp(mm2$theta["A", "B"]))
    ##   print(exp(mm2$theta["B", "A"]))
    ##   print(exp(mm2$theta["A", "A"]))
    ##   print(sampledGenotypes(dB))
    ## }
  
    ## if(isTRUE(all.equal(
    ##   mm2$transitionRateMatrix["A", "A, C"], 
    ##   exp(mm2$theta["C", "C"]) * exp(mm2$theta["A", "C"])))) {
    ##   cat("\n \n Here they are identical too 3_4")
    ##   print(exp(mm2$theta["A", "C"]))
    ##   print(exp(mm2$theta["C", "A"]))
    ##   print(exp(mm2$theta["C", "C"]))
    ##   print(sampledGenotypes(dB))
    ## }
    
    ## if(isTRUE(all.equal(
    ##   mm2$transitionRateMatrix["A", "A, C"], 
    ##   exp(mm2$theta["A", "A"]) * exp(mm2$theta["A", "C"])))) {
    ##   cat("\n \n Here they are identical too 3_5")
    ##   print(exp(mm2$theta["A", "C"]))
    ##   print(exp(mm2$theta["C", "A"]))
    ##   print(exp(mm2$theta["A", "A"]))
    ##   print(sampledGenotypes(dB))
    ## }
    
    ## if(isTRUE(all.equal(
    ##   mm2$transitionRateMatrix["A", "A, C"], 
    ##   exp(mm2$theta["A", "A"]) * exp(mm2$theta["C", "A"])))) {
    ##   cat("\n \n Here they are identical too 3_6")
    ##   print(exp(mm2$theta["A", "C"]))
    ##   print(exp(mm2$theta["C", "A"]))
    ##   print(exp(mm2$theta["A", "A"]))
    ##   print(sampledGenotypes(dB))
    ## }
})




## Tests with Schill's data
## this is from evamtools/R
## path_to_Schill_data <- "../inst/miscell/MHN_data/"
path_to_Schill_data <- "../../inst/miscell/MHN_data/"

test_that("Working with breast cancer: identical results from different algos", {
    Dat <- readRDS(file = paste0(path_to_Schill_data, "BreastCancer.rds"))

    pD <- Data.to.pD(Dat)
    Theta.BC <- Learn.MHN(pD, lambda = 0.01)

    colnames(Theta.BC) <- colnames(Dat)
    rownames(Theta.BC) <- colnames(Theta.BC)

    ## smaller:
    e1 <- Theta.BC[1:4, 1:4]

    system.time(te1_1 <- theta_to_trans_rate_1(e1))
    system.time(te1_2 <- theta_to_trans_rate(e1))
    system.time(te1_3 <- theta_to_trans_rate_3(e1))

    expect_identical(te1_1, te1_2)
    expect_identical(te1_1, te1_3)

    tt31 <- theta_to_trans_rate_3(Theta.BC,
                                  inner_transition = inner_transitionRate_3_1)

    ## slow
    system.time(tt1 <- theta_to_trans_rate_1(Theta.BC))
    ## these two are very slow
    ## system.time(tt2 <- theta_to_trans_rate(Theta.BC))
    ## system.time(tt22 <- theta_to_trans_rate(Theta.BC,
    ##                                    inner_transition = inner_transitionRate_2))
    ## but see this! Culprit was accessing a data frame?
    system.time(
        tt31 <-
            theta_to_trans_rate_3(Theta.BC,
                                  inner_transition = inner_transitionRate_3_1))
    ## but this is slower by x2. Why? setdiff, I guess
    system.time(
        tt32 <-
            theta_to_trans_rate_3(Theta.BC,
                                  inner_transition = inner_transitionRate_3_2))
    expect_identical(tt1, tt31)
    expect_identical(tt31, tt32)
})


## From Schill's ExampleApplications.R
## additional tests and example calls
  

test_that("do_MHN and do_MHN2 identical in Schill's data", {
    DatBS <- readRDS(file = paste0(path_to_Schill_data,
                                   "BreastCancer.rds")) [1:50, 1:4]
    DatCC <- readRDS(file = paste0(path_to_Schill_data,
                                   "ColorectalCancer.rds"))[1:300, 1:9]
    DatRCC <- readRDS(file = paste0(path_to_Schill_data,
                                    "RenalCellCarcinoma.rds"))[, 1:9]

    datasets <- list(DatBS, DatCC, DatRCC)
    datasets <- sapply(datasets, function(data) {
        colnames(data) <- LETTERS[1:ncol(data)]
        data <- as.matrix(data)
        return(data)
    })


    for (data in datasets) {
        mhn0 <- do_MHN(data)
        mhnSM <- do_MHN2(data)

        expect_identical(mhn0$theta, mhnSM$theta)
        expect_identical(mhn0$transitionRateMatrix,
                         as.matrix(mhnSM$transitionRateMatrix))
        expect_equal(mhn0$transitionMatrixTimeDiscretized,
                     as.matrix(mhnSM$transitionMatrixTimeDiscretized))
        expect_equal(mhn0$transitionMatrixCompExp,
                     as.matrix(mhnSM$transitionMatrixCompExp))
    }


    Dat11 <- readRDS(file = paste0(path_to_Schill_data, "Glioblastoma.rds"))

    for (i in 1:10) {
        Dat1 <- Dat11[, sample(1:ncol(Dat11), 6)]
        Dat1 <- as.matrix(Dat1)
        mhn0 <- do_MHN(Dat1)
        mhnSM <- do_MHN2(Dat1)
        expect_identical(mhn0$theta, mhnSM$theta)
        expect_identical(mhn0$transitionRateMatrix,
                         as.matrix(mhnSM$transitionRateMatrix))
        expect_equal(mhn0$transitionMatrixTimeDiscretized,
                     as.matrix(mhnSM$transitionMatrixTimeDiscretized))
        expect_equal(mhn0$transitionMatrixCompExp,
                     as.matrix(mhnSM$transitionMatrixCompExp))
    }
})
