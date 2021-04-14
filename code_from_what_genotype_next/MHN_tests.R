source("schill-trans-mat.R")


local({
  
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
  stopifnot(inherits(mm22$transitionRateMatrix, "dgCMatrix"))
  stopifnot(inherits(mm22$transitionMatrixTimeDiscretized, "dgCMatrix"))
  stopifnot(inherits(mm22$transitionMatrixCompExp, "dgCMatrix"))
  
  ## To compare same values, need to convert into matrices
  stopifnot(identical(mm2$transitionRateMatrix,
                      as.matrix(mm22$transitionRateMatrix)))
  ## Identical would fail
  stopifnot(all.equal(mm2$transitionMatrixCompExp,
                      as.matrix(mm22$transitionMatrixCompExp)))
  
  stopifnot(all.equal(mm2$transitionMatrixTimeDiscretized,
                      as.matrix(mm22$transitionMatrixTimeDiscretized)))
  
  
  ## This would fail, of course unless we turned them into matrices
  stopifnot(!identical(mm2$transitionRateMatrix,
                       mm22$transitionRateMatrix))
  stopifnot(!identical(mm2$transitionMatrixTimeDiscretized,
                       mm22$transitionMatrixTimeDiscretized))
  stopifnot(!identical(mm2$transitionMatrixCompExp,
                       mm22$transitionMatrixCompExp))
  
  
  ## Does not seem great for how we go
  ## ABCE to ABCDE and ABCEF. But probably a hard model anyway,
  
  ## From A to A, B:
  mm2$transitionMatrixCompExp["A", "A, B"]
  mm2$transitionRateMatrix["A", "A, B"]
  exp(mm2$theta["B", "B"]) * exp(mm2$theta["B", "A"])
  
  ## same, of course
  stopifnot(all.equal( mm2$transitionRateMatrix["A", "A, B"],
                       exp(mm2$theta["B", "B"]) * exp(mm2$theta["B", "A"])))
  
  ## not the same, but not always
  ## stopifnot(!isTRUE(all.equal(
  ##     mm2$transitionRateMatrix["A", "A, B"], 
  ##     exp(mm2$theta["B", "B"]) * exp(mm2$theta["A", "B"]))))
  
  ## stopifnot(!isTRUE(all.equal(
  ##     mm2$transitionRateMatrix["A", "A, B"], 
  ##     exp(mm2$theta["A", "A"]) * exp(mm2$theta["A", "B"]))))
  
  ## stopifnot(!isTRUE(all.equal(
  ##     mm2$transitionRateMatrix["A", "A, B"], 
  ##     exp(mm2$theta["A", "A"]) * exp(mm2$theta["B", "A"]))))
  
  if(isTRUE(all.equal(
    mm2$transitionRateMatrix["A", "A, B"], 
    exp(mm2$theta["B", "B"]) * exp(mm2$theta["A", "B"])))) {
    cat("\n \n Here they are identical too 3_1")
    print(exp(mm2$theta["A", "B"]))
    print(exp(mm2$theta["B", "A"]))
    print(exp(mm2$theta["B", "B"]))
    print(sampledGenotypes(dB))
  }
  if(isTRUE(all.equal(
    mm2$transitionRateMatrix["A", "A, B"], 
    exp(mm2$theta["A", "A"]) * exp(mm2$theta["A", "B"])))) {
    cat("\n \n Here they are identical too 3_2")
    print(exp(mm2$theta["A", "B"]))
    print(exp(mm2$theta["B", "A"]))
    print(exp(mm2$theta["A", "A"]))
    print(sampledGenotypes(dB))
  }
  if(isTRUE(all.equal(
    mm2$transitionRateMatrix["A", "A, B"], 
    exp(mm2$theta["A", "A"]) * exp(mm2$theta["B", "A"])))) {
    cat("\n \n Here they are identical too 3_3")
    print(exp(mm2$theta["A", "B"]))
    print(exp(mm2$theta["B", "A"]))
    print(exp(mm2$theta["A", "A"]))
    print(sampledGenotypes(dB))
  }
  
  
  
  
  
  ## ## same, of course
  ## stopifnot(all.equal( mm2$transitionRateMatrix["A", "A, C"],
  ##                     exp(mm2$theta["C", "C"]) * exp(mm2$theta["C", "A"])))
  ## not the same, but not always
  if(isTRUE(all.equal(
    mm2$transitionRateMatrix["A", "A, C"], 
    exp(mm2$theta["C", "C"]) * exp(mm2$theta["A", "C"])))) {
    cat("\n \n Here they are identical too 3_4")
    print(exp(mm2$theta["A", "C"]))
    print(exp(mm2$theta["C", "A"]))
    print(exp(mm2$theta["C", "C"]))
    print(sampledGenotypes(dB))
  }
  
  if(isTRUE(all.equal(
    mm2$transitionRateMatrix["A", "A, C"], 
    exp(mm2$theta["A", "A"]) * exp(mm2$theta["A", "C"])))) {
    cat("\n \n Here they are identical too 3_5")
    print(exp(mm2$theta["A", "C"]))
    print(exp(mm2$theta["C", "A"]))
    print(exp(mm2$theta["A", "A"]))
    print(sampledGenotypes(dB))
  }
  
  if(isTRUE(all.equal(
    mm2$transitionRateMatrix["A", "A, C"], 
    exp(mm2$theta["A", "A"]) * exp(mm2$theta["C", "A"])))) {
    cat("\n \n Here they are identical too 3_6")
    print(exp(mm2$theta["A", "C"]))
    print(exp(mm2$theta["C", "A"]))
    print(exp(mm2$theta["A", "A"]))
    print(sampledGenotypes(dB))
  }
  
  
  ## stopifnot(!isTRUE(all.equal(
  ##     mm2$transitionRateMatrix["A", "A, C"], 
  ##     exp(mm2$theta["C", "C"]) * exp(mm2$theta["A", "C"]))))
  
  ## Same story as above.
  ## stopifnot(!isTRUE(all.equal(
  ##     mm2$transitionRateMatrix["A", "A, C"], 
  ##     exp(mm2$theta["A", "A"]) * exp(mm2$theta["A", "C"]))))
  
  ## stopifnot(!isTRUE(all.equal(
  ##     mm2$transitionRateMatrix["A", "A, C"], 
  ##     exp(mm2$theta["A", "A"]) * exp(mm2$theta["C", "A"]))))
  
  
  ## same, of course
  stopifnot(all.equal( mm2$transitionRateMatrix["A, B", "A, B, C"],
                       exp(mm2$theta["C", "C"]) *
                         exp(mm2$theta["C", "A"]) *
                         exp(mm2$theta["C", "B"]) 
  ))
  
  stopifnot(all.equal( mm2$transitionRateMatrix["A, B", "A, B, D"],
                       exp(mm2$theta["D", "D"]) *
                         exp(mm2$theta["D", "A"]) *
                         exp(mm2$theta["D", "B"]) 
  ))
  
  ## A very extreme one, from a non-existing genotype
  stopifnot(all.equal( mm2$transitionRateMatrix["A, C", "A, B, C"],
                       exp(mm2$theta["B", "B"]) *
                         exp(mm2$theta["B", "A"]) *
                         exp(mm2$theta["B", "C"]) 
  ))
  
  
  stopifnot(all.equal( mm2$transitionRateMatrix["B, C, D", "A, B, C, D"],
                       exp(mm2$theta["A", "A"]) *
                         exp(mm2$theta["A", "B"]) *
                         exp(mm2$theta["A", "D"]) *
                         exp(mm2$theta["A", "C"])))
  
  stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C", "A, B, C, D"],
                       exp(mm2$theta["D", "D"]) *
                         exp(mm2$theta["D", "B"]) *
                         exp(mm2$theta["D", "A"]) *
                         exp(mm2$theta["D", "C"])))
  
  stopifnot(all.equal( mm2$transitionRateMatrix["A, B, D", "A, B, C, D"],
                       exp(mm2$theta["C", "C"]) *
                         exp(mm2$theta["C", "B"]) *
                         exp(mm2$theta["C", "D"]) *
                         exp(mm2$theta["C", "A"])))
  
  stopifnot(all.equal( mm2$transitionRateMatrix["A, C, D", "A, B, C, D"],
                       exp(mm2$theta["B", "B"]) *
                         exp(mm2$theta["B", "A"]) *
                         exp(mm2$theta["B", "D"]) *
                         exp(mm2$theta["B", "C"])))
  
  
  ## 5 genes
  stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C, E", "A, B, C, D, E"],
                       exp(mm2$theta["D", "D"]) *
                         exp(mm2$theta["D", "A"]) *
                         exp(mm2$theta["D", "B"]) *
                         exp(mm2$theta["D", "C"]) *                        
                         exp(mm2$theta["D", "E"])))
  
  stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C, E", "A, B, C, E, F"],
                       exp(mm2$theta["F", "F"]) *
                         exp(mm2$theta["F", "A"]) *
                         exp(mm2$theta["F", "B"]) *
                         exp(mm2$theta["F", "C"]) *                        
                         exp(mm2$theta["F", "E"])))
  
  stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C, F", "A, B, C, E, F"],
                       exp(mm2$theta["E", "E"]) *
                         exp(mm2$theta["E", "A"]) *
                         exp(mm2$theta["E", "B"]) *
                         exp(mm2$theta["E", "C"]) *                        
                         exp(mm2$theta["E", "F"])))
  
  stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C, F", "A, B, C, D, F"],
                       exp(mm2$theta["D", "D"]) *
                         exp(mm2$theta["D", "A"]) *
                         exp(mm2$theta["D", "B"]) *
                         exp(mm2$theta["D", "C"]) *                        
                         exp(mm2$theta["D", "F"])))
  
  ## six genes
  stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C, D, E", "A, B, C, D, E, F"],
                       exp(mm2$theta["F", "F"]) *
                         exp(mm2$theta["F", "A"]) *
                         exp(mm2$theta["F", "B"]) *
                         exp(mm2$theta["F", "C"]) *
                         exp(mm2$theta["F", "D"]) *                        
                         exp(mm2$theta["F", "E"])))
  
  stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C, D, F", "A, B, C, D, E, F"],
                       exp(mm2$theta["E", "E"]) *
                         exp(mm2$theta["E", "A"]) *
                         exp(mm2$theta["E", "B"]) *
                         exp(mm2$theta["E", "C"]) *
                         exp(mm2$theta["E", "D"]) *                        
                         exp(mm2$theta["E", "F"])))
  
  
  stopifnot(all.equal( mm2$transitionRateMatrix["A, B, C, E, F", "A, B, C, D, E, F"],
                       exp(mm2$theta["D", "D"]) *
                         exp(mm2$theta["D", "A"]) *
                         exp(mm2$theta["D", "B"]) *
                         exp(mm2$theta["D", "C"]) *
                         exp(mm2$theta["D", "E"]) *                        
                         exp(mm2$theta["D", "F"])))
  
  ## six genes, with sparse. Totally unnecessary since they were identical to non-sparse
  stopifnot(all.equal( mm22$transitionRateMatrix["A, B, C, D, E", "A, B, C, D, E, F"],
                       exp(mm22$theta["F", "F"]) *
                         exp(mm22$theta["F", "A"]) *
                         exp(mm22$theta["F", "B"]) *
                         exp(mm22$theta["F", "C"]) *
                         exp(mm22$theta["F", "D"]) *                        
                         exp(mm22$theta["F", "E"])))
  
  stopifnot(all.equal( mm22$transitionRateMatrix["A, B, C, D, F", "A, B, C, D, E, F"],
                       exp(mm22$theta["E", "E"]) *
                         exp(mm22$theta["E", "A"]) *
                         exp(mm22$theta["E", "B"]) *
                         exp(mm22$theta["E", "C"]) *
                         exp(mm22$theta["E", "D"]) *                        
                         exp(mm22$theta["E", "F"])))
  
  
  stopifnot(all.equal( mm22$transitionRateMatrix["A, B, C, E, F", "A, B, C, D, E, F"],
                       exp(mm22$theta["D", "D"]) *
                         exp(mm22$theta["D", "A"]) *
                         exp(mm22$theta["D", "B"]) *
                         exp(mm22$theta["D", "C"]) *
                         exp(mm22$theta["D", "E"]) *                        
                         exp(mm22$theta["D", "F"])))
  
  
  
  tmp <- round(mm2$transitionRateMatrix, 3)
  
  ## And output matches the observed
  ## In particular: A, B to A,B,C vs A, B, D
})

