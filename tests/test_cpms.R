source("../code_from_what_genotype_next/code-all-methods-minimal.R")

####################### Tests #################

## Not using a test_that block here. But doing it locally
local({
  
  
  ex0 <- list(edges = data.frame(
    From = c("Root", "B", "C", "D"),
    To =   c("B", "C", "D", "A"),
    rerun_lambda = c(1, 2, 3, 4)
  ))
  
  
  ex1 <- list(edges = data.frame(
    From = c("Root", "Root", "A", "B", "C"),
    To =   c("A", "B", "C", "C", "D"),
    rerun_lambda = c(1, 2, 3, 3, 4)
  ))
  
  ex2 <- list(edges = data.frame(
    From = c("Root", "Root", "A", "B"),
    To =   c("A", "B", "C", "D"),
    rerun_lambda = c(10, 11, 12, 14)
  ))
  
  ex3 <- list(edges = data.frame(
    From = c("Root", "Root", "A", "C"),
    To =   c("A", "B", "C", "D"),
    rerun_lambda = c(2, 3, 4, 5)
  ))
  
  ex4 <- list(edges = data.frame(
    From = c("Root", "Root", "A", "A"),
    To =   c("A", "B", "C", "D"),
    rerun_lambda = c(1, 2, 3, 4)
  ))
  
  ex5 <- list(edges = data.frame(
    From = c("Root", "Root", "Root", "A", "B", "C"),
    To =   c("A", "B", "C", "D", "D", "D"),
    rerun_lambda = c(1, 2, 3, 4, 4, 4)
  ))
  ## yes. this is wrong. Tested below.
  ex5e <- list(edges = data.frame(
    From = c("Root", "Root", "Root", "A", "B", "C"),
    To =   c("A", "B", "C", "D", "D", "D"),
    rerun_lambda = c(1, 2, 3, 4, 4, 5)
  ))
  
  ex6 <- list(edges = data.frame(
    From = c("Root", "Root", "Root", "A", "B"),
    To =   c("A", "B", "D", "C", "C"),
    rerun_lambda = c(1, 2, 3, 5, 5)
  ))
  
  ex7 <- list(edges = data.frame(
    From = c("Root", "E", "E", "E", "A", "B"),
    To =   c("E",    "A", "B", "D", "C", "C"),
    rerun_lambda = c(1, 2, 3, 4, 5, 5)
  ))
  
  ex8 <- list(edges = data.frame(
    From = c("Root", "Root", "A", "B", "C", "C"),
    To =   c("A",     "B",   "C", "C", "D",  "E"),
    rerun_lambda = c(1, 2, 3, 3, 4, 5)
  ))
  
  ex9 <- list(edges = data.frame(
    From = c("Root", "A", "B", "B", "B"),
    To =   c("A",    "B", "C", "D", "E"),
    rerun_lambda = c(6, 7, 8, 9, 10)
  ))
  ## wrong, triggers error
  ex10e <- list(edges = data.frame(
    From = c("Root", "A", "A", "B", "C", "D"),
    To =   c("A",    "B", "C", "D", "D", "E"),
    rerun_lambda = c(1, 2, 3, 4, 5, 5)
  ))
  
  ex10 <- list(edges = data.frame(
    From = c("Root", "A", "A", "B", "C", "D"),
    To =   c("A",    "B", "C", "D", "D", "E"),
    rerun_lambda = c(1, 2, 3, 4, 4, 5)
  ))
  
  
  ex11 <- list(edges = data.frame(
    From = c(rep("Root", 4), "A", "B", "C", "D"),
    To =   c("A",  "B", "C", "D", rep("E", 4)),
    rerun_lambda = c(3, 4, 5, 6, rep(7, 4))
  ))
  
  
  
  oex0 <- cpm_access_genots_paths_w(ex0)
  oex1 <- cpm_access_genots_paths_w(ex1)
  oex2 <- cpm_access_genots_paths_w(ex2)
  oex3 <- cpm_access_genots_paths_w(ex3)
  oex4 <- cpm_access_genots_paths_w(ex4)
  oex5 <- cpm_access_genots_paths_w(ex5)
  expect_error(cpm_access_genots_paths_w(ex5e), "Different lambda", fixed = TRUE)
  
  oex6 <- cpm_access_genots_paths_w(ex6)
  oex7 <- cpm_access_genots_paths_w(ex7)
  oex8 <- cpm_access_genots_paths_w(ex8)
  oex9 <- cpm_access_genots_paths_w(ex9)
  expect_error(cpm_access_genots_paths_w(ex10e), "Different lambda", fixed = TRUE)
  oex10 <- cpm_access_genots_paths_w(ex10)
  oex11 <- cpm_access_genots_paths_w(ex11)
  
  
  expect_equivalent(oex0$weighted_paths[, 2], 1)
  
  expect_equivalent(oex1$weighted_paths[, 2],
                    c(1/3 * 1 * 1 * 1, 2/3 * 1 * 1 * 1))
  
  expect_equivalent(oex2$weighted_paths[, 2],
                    c(10/21 * 11/23 * 12/26 * 1,
                      10/21 * 11/23 * 14/26 * 1,
                      10/21 * 12/23 * 1     * 1,
                      11/21 * 10/24 * 12/26 * 1,
                      11/21 * 10/24 * 14/26 * 1,
                      11/21 * 14/24 * 1     * 1))
  
  expect_equivalent(oex3$weighted_paths[, 2],
                    c(2/5 * 3/7 * 1 * 1,
                      2/5 * 4/7 * 3/8 * 1,
                      2/5 * 4/7 * 5/8 * 1,
                      3/5 * 1 *   1   * 1
                    ))
  
  expect_equivalent(oex4$weighted_paths[, 2],
                    c(1/3 * 2/9 * 3/7 * 1,
                      1/3 * 2/9 * 4/7 * 1,
                      1/3 * 3/9 * 2/6 * 1,
                      1/3 * 3/9 * 4/6 * 1,
                      1/3 * 4/9 * 2/5 * 1,
                      1/3 * 4/9 * 3/5 * 1,
                      2/3 * 1 *   3/7 * 1,
                      2/3 * 1 *   4/7 * 1
                    ))
  
  expect_equivalent(oex5$weighted_paths[, 2],
                    c(1/6 * 2/5 * 1 * 1,
                      1/6 * 3/5 * 1 * 1,
                      2/6 * 1/4 * 1 * 1,
                      2/6 * 3/4 * 1 * 1,
                      3/6 * 1/3 * 1 * 1,
                      3/6 * 2/3 * 1 * 1
                    ))
  
  expect_equivalent(oex6$weighted_paths[, 2],
                    c(1/6 * 2/5 * 3/8 * 1,
                      1/6 * 2/5 * 5/8 * 1,
                      1/6 * 3/5 * 1 * 1,
                      2/6 * 1/4 * 3/8 * 1,
                      2/6 * 1/4 * 5/8 * 1,
                      2/6 * 3/4 * 1 * 1,
                      3/6 * 1/3 * 1 * 1,
                      3/6 * 2/3 * 1 * 1
                    ))
  
  
  expect_equivalent(oex7$weighted_paths[, 2],
                    c(1 * 2/9 * 3/7 * 4/9,
                      1 * 2/9 * 3/7 * 5/9,
                      1 * 2/9 * 4/7 * 1,
                      1 * 3/9 * 2/6 * 4/9,
                      1 * 3/9 * 2/6 * 5/9,
                      1 * 3/9 * 4/6 * 1,
                      1 * 4/9 * 2/5 * 1,
                      1 * 4/9 * 3/5 * 1
                    ))
  
  expect_equivalent(oex8$weighted_paths[, 2],
                    c(1/3 * 1 * 1 * 4/9,
                      1/3 * 1 * 1 * 5/9,
                      2/3 * 1 * 1 * 4/9,
                      2/3 * 1 * 1 * 5/9
                    ))
  
  
  expect_equivalent(oex9$weighted_paths[, 2],
                    c(6/6 * 7/7 * 8/27 * 9/19 * 10/10,
                      6/6 * 7/7 * 8/27 * 10/19 * 9/9,                    
                      6/6 * 7/7 * 9/27 * 8/18 * 10/10,
                      6/6 * 7/7 * 9/27 * 10/18 * 8/8,                   
                      6/6 * 7/7 * 10/27 * 8/17 * 9/9,
                      6/6 * 7/7 * 10/27 * 9/17 * 8/8
                    ))
  
  expect_equivalent(oex10$weighted_paths[, 2],
                    c(1/1 * 2/5 * 3/3 * 4/4 * 5/5,
                      1/1 * 3/5 * 2/2 * 4/4 * 5/5                                        
                    ))
  
  expect_equivalent(oex11$weighted_paths[, 2],
                    c(3/18 * 4/15 * 5/11 * 7/7,
                      3/18 * 4/15 * 6/11 * 7/7,
                      3/18 * 5/15 * 4/10 * 7/7,
                      3/18 * 5/15 * 6/10 * 7/7,
                      3/18 * 6/15 * 4/9 * 7/7,
                      3/18 * 6/15 * 5/9 * 7/7,
                      4/18 * 3/14 * 5/11 * 7/7,                    
                      4/18 * 3/14 * 6/11 * 7/7,
                      4/18 * 5/14 * 3/9 * 7/7,
                      4/18 * 5/14 * 6/9 * 7/7,
                      4/18 * 6/14 * 3/8 * 7/7,
                      4/18 * 6/14 * 5/8 * 7/7,                    
                      5/18 * 3/13 * 4/10 * 7/7,
                      5/18 * 3/13 * 6/10 * 7/7,
                      5/18 * 4/13 * 3/9 * 7/7,
                      5/18 * 4/13 * 6/9 * 7/7,
                      5/18 * 6/13 * 3/7 * 7/7,
                      5/18 * 6/13 * 4/7 * 7/7,                    
                      6/18 * 3/12 * 4/9 * 7/7,
                      6/18 * 3/12 * 5/9 * 7/7,
                      6/18 * 4/12 * 3/8 * 7/7,
                      6/18 * 4/12 * 5/8 * 7/7,
                      6/18 * 5/12 * 3/7 * 7/7,
                      6/18 * 5/12 * 4/7 * 7/7                    
                    ))
  
  
  
  
  ## Test the simplified code
  oex0_simplified <- cpm_access_genots_paths_w_simplified(ex0)
  oex1_simplified <- cpm_access_genots_paths_w_simplified(ex1)
  oex2_simplified <- cpm_access_genots_paths_w_simplified(ex2)
  oex3_simplified <- cpm_access_genots_paths_w_simplified(ex3)
  oex4_simplified <- cpm_access_genots_paths_w_simplified(ex4)
  oex5_simplified <- cpm_access_genots_paths_w_simplified(ex5)
  expect_error(cpm_access_genots_paths_w_simplified(ex5e), "Different lambda", fixed = TRUE)
  
  ## We are now using sparseMatrices: fuller testing component by component below
  
  oex6_simplified <- cpm_access_genots_paths_w_simplified(ex6)
  oex7_simplified <- cpm_access_genots_paths_w_simplified(ex7)
  oex8_simplified <- cpm_access_genots_paths_w_simplified(ex8)
  oex9_simplified <- cpm_access_genots_paths_w_simplified(ex9)
  expect_error(cpm_access_genots_paths_w_simplified(ex10e), "Different lambda", fixed = TRUE)
  oex10_simplified <- cpm_access_genots_paths_w_simplified(ex10)
  oex11_simplified <- cpm_access_genots_paths_w_simplified(ex11)
  
  expect_equal(oex0[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
               lapply(oex0_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))
  
  expect_equal(oex1[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
               lapply(oex1_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))
  
  expect_equal(oex2[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
               lapply(oex2_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))
  
  expect_equal(oex3[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
               lapply(oex3_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))
  
  expect_equal(oex4[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
               lapply(oex4_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))
  
  expect_equal(oex5[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
               lapply(oex5_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))
  
  expect_equal(oex6[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
               lapply(oex6_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))
  
  expect_equal(oex7[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
               lapply(oex7_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))
  
  expect_equal(oex8[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
               lapply(oex8_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))
  
  expect_equal(oex9[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
               lapply(oex9_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))
  
  expect_equal(oex10[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
               lapply(oex10_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))
  
  expect_equal(oex11[c("fgraph", "weighted_fgraph", "trans_mat_genots")],
               lapply(oex11_simplified[c("fgraph", "weighted_fgraph", "trans_mat_genots")], as.matrix))
  
  gacc1 <- list("A", "B", "D", c("A", "B"), c("A", "D"), c("B", "D"), c("A", 
                                                                        "C"), c("A", "B", "D"), c("A", "B", "C"), c("A", "C", "D"), c("A", 
                                                                                                                                      "B", "C", "D"))
  
  am1 <- unrestricted_fitness_graph(gacc1)
  am1M <- unrestricted_fitness_graph_sparseM(gacc1)
  am1Mm <- as.matrix(am1M)
  storage.mode(am1Mm) <- "integer"
  expect_identical(am1, am1Mm)
  
  wg <- structure(list(To = c("A", "B", "C", "D"), OT_edgeWeight = c(0.525915054637741, 
                                                                     0.101508072999909, 0.108065026764796, 0.302542959038882)), row.names = c("A", 
                                                                                                                                              "B", "C", "D"), class = "data.frame")
  
  tm1 <- transition_fg(am1, wg)
  tm1M <- transition_fg_sparseM(am1M, wg)
  tm1M <- as.matrix(tm1M)
  expect_identical(tm1, tm1M)
  
  ## To test with a full run. Checking sparse matrix implementation
  cpm_out_others1 <- list(OT = list(edges = structure(list(From = c("Root", "Root", 
"A", "Root"), To = c("A", "B", "C", "D"), edge = c("Root -> A", 
"Root -> B", "A -> C", "Root -> D"), OT_edgeBootFreq = c(NA, 
NA, NA, NA), OT_edgeWeight = c(0.525915054637741, 0.101508072999909, 
0.108065026764796, 0.302542959038882), OT_obsMarginal = c(0.56, 
0.18, 0.14, 0.36), OT_predMarginal = c(0.56, 0.18, 0.139999436433903, 
0.36)), class = "data.frame", row.names = c("A", "B", "C", "D"
)), consensus = NA, OT_error.fun = "std", ot.boot.original = NA, 
    genots_predicted = NA, genots_observed = NA, two_way_predicted = NA, 
    two_way_observed = NA), CAPRESE = list(edges = structure(list(
    From = c("Root", "Root", "Root", "A"), To = c("A", "B", "D", 
    "C"), edge = c("Root -> A", "Root -> B", "Root -> D", "A -> C"
    ), CAPRESE_edgeBootFreq = c(NA, NA, NA, NA), CAPRESE_pr = c(NA, 
    NA, NA, 0.223188719001567), CAPRESE_tp = c(NA, NA, NA, 1), 
    CAPRESE_hg = c(NA, NA, NA, 0.0948328267477203)), class = "data.frame", row.names = c(NA, 
-4L)), CAPRESE_eloss = NA, CAPRESE_prederr = NA, CAPRESE_posterr = NA, 
    CAPRESE_useless_extra = list(selective_advantage = list(caprese = structure(list(
        SELECTS = "variant A", SELECTED = "variant C", OBS.SELECTS = 28, 
        OBS.SELECTED = 7, TEMPORAL.PRIORITY = 1, PROBABILITY.RAISING = 0.223188719001567, 
        HYPERGEOMETRIC = 0.0948328267477203), row.names = "1", class = "data.frame")), 
        bootstrap_scores = NA), TRONCO_model_boot = NULL, TRONCO_model = list(
        genotypes = structure(c(0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 
        1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 
        0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 
        0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 
        0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 
        0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 
        0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 
        0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 
        0L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 
        1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L), .Dim = c(50L, 
        4L), .Dimnames = list(c("1", "2", "3", "4", "5", "6", 
        "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", 
        "17", "18", "19", "20", "21", "22", "23", "24", "25", 
        "26", "27", "28", "29", "30", "31", "32", "33", "34", 
        "35", "36", "37", "38", "39", "40", "41", "42", "43", 
        "44", "45", "46", "47", "48", "49", "50"), c("G1", "G2", 
        "G3", "G4"))), annotations = structure(c("variant", "variant", 
        "variant", "variant", "A", "B", "C", "D"), .Dim = c(4L, 
        2L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("type", 
        "event"))), types = structure("Darkgreen", .Dim = c(1L, 
        1L), .Dimnames = list("variant", "color")), hypotheses = NA, 
        confidence = structure(list(structure(c(0, 0.321428571428571, 
        0.25, 0.642857142857143, 1, 0, 0.777777777777778, 1, 
        1, 1, 0, 1, 1, 0.5, 0.388888888888889, 0), .Dim = c(4L, 
        4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("G1", 
        "G2", "G3", "G4"))), structure(c(-1, -0.125943176525384, 
        0.132346627230226, -0.00509762193462174, -0.171557211613104, 
        -1, -0.123178689550371, -0.281478505423865, 0.223188719001567, 
        -0.125943176525384, -1, -0.14193770830939, -0.00649653637701445, 
        -0.255488752279943, -0.123178689550371, -1), .Dim = c(4L, 
        4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("G1", 
        "G2", "G3", "G4"))), structure(c(0, 0.657728939304999, 
        0.0948328267477203, 0.403066884258277, 0.657728939304999, 
        0, 0.369777142376588, 0.707693654795176, 0.0948328267477203, 
        0.369777142376588, 0, 0.494537285101578, 0.403066884258277, 
        0.707693654795176, 0.494537285101578, 0), .Dim = c(4L, 
        4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("G1", 
        "G2", "G3", "G4")))), .Dim = c(3L, 1L), .Dimnames = list(
            c("temporal priority", "probability raising", "hypergeometric test"
            ), "confidence")), model = list(caprese = list(probabilities = list(
            probabilities.observed = list(marginal.probs = structure(c(0.56, 
            0.18, 0.14, 0.36), .Dim = c(4L, 1L), .Dimnames = list(
                c("G1", "G2", "G3", "G4"), "marginal probability")), 
                joint.probs = structure(c(0.56, 0.08, 0.1, 0.2, 
                0.08, 0.18, 0.02, 0.04, 0.1, 0.02, 0.14, 0.04, 
                0.2, 0.04, 0.04, 0.36), .Dim = c(4L, 4L), .Dimnames = list(
                  c("G1", "G2", "G3", "G4"), c("G1", "G2", "G3", 
                  "G4"))), conditional.probs = structure(c(1, 
                1, 0.178571428571429, 1), .Dim = c(4L, 1L), .Dimnames = list(
                  c("G1", "G2", "G3", "G4"), "conditional probability"))), 
            probabilities.fit = list(estimated.marginal.probs = NA, 
                estimated.joint.probs = NA, estimated.conditional.probs = NA)), 
            parents.pos = structure(c(-1, -1, 1, -1), .Dim = c(4L, 
            1L), .Dimnames = list(c("G1", "G2", "G3", "G4"), 
                "parents")), error.rates = list(error.fp = NA, 
                error.fn = NA), adj.matrix = list(adj.matrix.fit = structure(c(0, 
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), .Dim = c(4L, 
            4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), 
                c("G1", "G2", "G3", "G4")))), logLik = -110.377237990856)), 
        parameters = list(algorithm = "CAPRESE", lambda = 0.5, 
            silent = TRUE, error.rates = list(epos = 0, eneg = 0)), 
        execution.time = structure(c(user.self = 0.00199999999998113, 
        sys.self = 0, elapsed = 0.00100000000020373, user.child = 0, 
        sys.child = 0), class = "proc_time"))), CAPRI_BIC = list(
    edges = structure(list(From = c("Root", "Root", "Root", "Root"
    ), To = c("A", "B", "C", "D"), edge = c("Root -> A", "Root -> B", 
    "Root -> C", "Root -> D"), CAPRI_edgeBootFreq = c(NA, NA, 
    NA, NA), CAPRI_pr = c(NA_real_, NA_real_, NA_real_, NA_real_
    ), CAPRI_tp = c(NA_real_, NA_real_, NA_real_, NA_real_), 
        CAPRI_hg = c(NA_real_, NA_real_, NA_real_, NA_real_)), class = "data.frame", row.names = c(NA, 
    -4L)), CAPRI_eloss = NA, CAPRI_prederr = NA, CAPRI_posterr = NA, 
    CAPRI_useless_extra = list(selective_advantage = list(capri_bic = NULL), 
        bootstrap_scores = NA), TRONCO_model_boot = NULL, TRONCO_model = list(
        genotypes = structure(c(0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 
        1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 
        0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 
        0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 
        0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 
        0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 
        0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 
        0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 
        0L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 
        1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L), .Dim = c(50L, 
        4L), .Dimnames = list(c("1", "2", "3", "4", "5", "6", 
        "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", 
        "17", "18", "19", "20", "21", "22", "23", "24", "25", 
        "26", "27", "28", "29", "30", "31", "32", "33", "34", 
        "35", "36", "37", "38", "39", "40", "41", "42", "43", 
        "44", "45", "46", "47", "48", "49", "50"), c("G1", "G2", 
        "G3", "G4"))), annotations = structure(c("variant", "variant", 
        "variant", "variant", "A", "B", "C", "D"), .Dim = c(4L, 
        2L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("type", 
        "event"))), types = structure("Darkgreen", .Dim = c(1L, 
        1L), .Dimnames = list("variant", "color")), hypotheses = NA, 
        adj.matrix.prima.facie = structure(c(0, 0, 0, 0, 0, 0, 
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0), .Dim = c(4L, 4L), .Dimnames = list(
            c("G1", "G2", "G3", "G4"), c("G1", "G2", "G3", "G4"
            ))), confidence = structure(list(structure(c(1, 1, 
        1, 1, 9.39901768425988e-35, 1, 0.999999998490081, 1.54941877188073e-32, 
        9.94334396623681e-35, 1.53273973984097e-09, 1, 2.65146581110054e-34, 
        4.00518485843201e-31, 1, 1, 1), .Dim = c(4L, 4L), .Dimnames = list(
            c("G1", "G2", "G3", "G4"), c("G1", "G2", "G3", "G4"
            ))), structure(c(1, 0.999999999923137, 3.71813560969053e-15, 
        0.846484514790407, 0.999999999999868, 1, 0.999886735833476, 
        1, 8.73588365187455e-17, 0.999958011279254, 1, 0.999997831094051, 
        0.781442879272848, 1, 0.999999534542517, 1), .Dim = c(4L, 
        4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("G1", 
        "G2", "G3", "G4"))), structure(c(1, 0.872936649650717, 
        0.0948328267477203, 0.635507694679078, 0.872936649650717, 
        1, 0.369777142376588, 0.707693654795176, 0.0948328267477203, 
        0.369777142376588, 1, 0.494537285101578, 0.635507694679078, 
        0.707693654795176, 0.494537285101578, 1), .Dim = c(4L, 
        4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("G1", 
        "G2", "G3", "G4")))), .Dim = c(3L, 1L), .Dimnames = list(
            c("temporal priority", "probability raising", "hypergeometric test"
            ), "confidence")), model = list(capri_bic = list(
            probabilities = list(probabilities.observed = list(
                marginal.probs = structure(c(0.5534, 0.183, 0.1418, 
                0.3628), .Dim = c(4L, 1L), .Dimnames = list(c("G1", 
                "G2", "G3", "G4"), "marginal probability")), 
                joint.probs = structure(c(0.5534, 0.0776, 0.1008, 
                0.1978, 0.0776, 0.183, 0.0204, 0.0422, 0.1008, 
                0.0204, 0.1418, 0.041, 0.1978, 0.0422, 0.041, 
                0.3628), .Dim = c(4L, 4L), .Dimnames = list(c("G1", 
                "G2", "G3", "G4"), c("G1", "G2", "G3", "G4"))), 
                conditional.probs = structure(list(1, 1, 1, 1), .Dim = c(4L, 
                1L), .Dimnames = list(c("G1", "G2", "G3", "G4"
                ), "conditional probability"))), probabilities.fit = list(
                estimated.marginal.probs = NA, estimated.joint.probs = NA, 
                estimated.conditional.probs = NA)), parents.pos = structure(list(
                -1, -1, -1, -1), .Dim = c(4L, 1L), .Dimnames = list(
                c("G1", "G2", "G3", "G4"), "parents")), error.rates = list(
                error.fp = NA, error.fn = NA), adj.matrix = list(
                adj.matrix.pf = structure(c(0, 0, 0, 0, 0, 0, 
                0, 0, 1, 0, 0, 0, 0, 0, 0, 0), .Dim = c(4L, 4L
                ), .Dimnames = list(c("G1", "G2", "G3", "G4"), 
                  c("G1", "G2", "G3", "G4"))), adj.matrix.fit = structure(c(0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), .Dim = c(4L, 
                4L), .Dimnames = list(c("G1", "G2", "G3", "G4"
                ), c("G1", "G2", "G3", "G4")))), score = -118.609294356862, 
            logLik = -110.785248346005)), parameters = list(algorithm = "CAPRI", 
            command = "hc", regularization = "bic", do.boot = TRUE, 
            nboot = 100, pvalue = 0.05, min.boot = 3, min.stat = TRUE, 
            boot.seed = NULL, silent = TRUE, error.rates = list(
                epos = 0, eneg = 0), restart = 100), execution.time = structure(c(user.self = 0.0769999999999982, 
        sys.self = 0, elapsed = 0.077000000001135, user.child = 0, 
        sys.child = 0), class = "proc_time"))), CAPRI_AIC = list(
    edges = structure(list(From = c("Root", "Root", "Root", "Root"
    ), To = c("A", "B", "C", "D"), edge = c("Root -> A", "Root -> B", 
    "Root -> C", "Root -> D"), CAPRI_edgeBootFreq = c(NA, NA, 
    NA, NA), CAPRI_pr = c(NA_real_, NA_real_, NA_real_, NA_real_
    ), CAPRI_tp = c(NA_real_, NA_real_, NA_real_, NA_real_), 
        CAPRI_hg = c(NA_real_, NA_real_, NA_real_, NA_real_)), class = "data.frame", row.names = c(NA, 
    -4L)), CAPRI_eloss = NA, CAPRI_prederr = NA, CAPRI_posterr = NA, 
    CAPRI_useless_extra = list(selective_advantage = list(capri_aic = NULL), 
        bootstrap_scores = NA), TRONCO_model_boot = NULL, TRONCO_model = list(
        genotypes = structure(c(0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 
        1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 
        0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 
        0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 
        0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 
        0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 
        0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 
        0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 
        0L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 
        1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L), .Dim = c(50L, 
        4L), .Dimnames = list(c("1", "2", "3", "4", "5", "6", 
        "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", 
        "17", "18", "19", "20", "21", "22", "23", "24", "25", 
        "26", "27", "28", "29", "30", "31", "32", "33", "34", 
        "35", "36", "37", "38", "39", "40", "41", "42", "43", 
        "44", "45", "46", "47", "48", "49", "50"), c("G1", "G2", 
        "G3", "G4"))), annotations = structure(c("variant", "variant", 
        "variant", "variant", "A", "B", "C", "D"), .Dim = c(4L, 
        2L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("type", 
        "event"))), types = structure("Darkgreen", .Dim = c(1L, 
        1L), .Dimnames = list("variant", "color")), hypotheses = NA, 
        adj.matrix.prima.facie = structure(c(0, 0, 0, 0, 0, 0, 
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0), .Dim = c(4L, 4L), .Dimnames = list(
            c("G1", "G2", "G3", "G4"), c("G1", "G2", "G3", "G4"
            ))), confidence = structure(list(structure(c(1, 1, 
        1, 1, 9.81725428393167e-35, 1, 0.999990867323265, 8.96819397729642e-32, 
        9.74566638078606e-35, 9.23440871142156e-06, 1, 3.31407732958239e-34, 
        1.34669152568532e-32, 1, 1, 1), .Dim = c(4L, 4L), .Dimnames = list(
            c("G1", "G2", "G3", "G4"), c("G1", "G2", "G3", "G4"
            ))), structure(c(1, 0.999999998167757, 2.73780178422517e-19, 
        0.463492715389974, 0.999999986477575, 1, 0.998529402441165, 
        1, 5.49487631358828e-17, 0.99916699252184, 1, 0.999998976978537, 
        0.525815417121827, 1, 0.999999778739998, 1), .Dim = c(4L, 
        4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("G1", 
        "G2", "G3", "G4"))), structure(c(1, 0.657728939304999, 
        0.0948328267477203, 0.506373911111967, 0.657728939304999, 
        1, 0.369777142376588, 0.890399114532546, 0.0948328267477203, 
        0.369777142376588, 1, 0.445575084798027, 0.506373911111967, 
        0.890399114532546, 0.445575084798027, 1), .Dim = c(4L, 
        4L), .Dimnames = list(c("G1", "G2", "G3", "G4"), c("G1", 
        "G2", "G3", "G4")))), .Dim = c(3L, 1L), .Dimnames = list(
            c("temporal priority", "probability raising", "hypergeometric test"
            ), "confidence")), model = list(capri_aic = list(
            probabilities = list(probabilities.observed = list(
                marginal.probs = structure(c(0.5678, 0.1802, 
                0.1472, 0.3494), .Dim = c(4L, 1L), .Dimnames = list(
                  c("G1", "G2", "G3", "G4"), "marginal probability")), 
                joint.probs = structure(c(0.5678, 0.0864, 0.1066, 
                0.1978, 0.0864, 0.1802, 0.0222, 0.0374, 0.1066, 
                0.0222, 0.1472, 0.0406, 0.1978, 0.0374, 0.0406, 
                0.3494), .Dim = c(4L, 4L), .Dimnames = list(c("G1", 
                "G2", "G3", "G4"), c("G1", "G2", "G3", "G4"))), 
                conditional.probs = structure(list(1, 1, 1, 1), .Dim = c(4L, 
                1L), .Dimnames = list(c("G1", "G2", "G3", "G4"
                ), "conditional probability"))), probabilities.fit = list(
                estimated.marginal.probs = NA, estimated.joint.probs = NA, 
                estimated.conditional.probs = NA)), parents.pos = structure(list(
                -1, -1, -1, -1), .Dim = c(4L, 1L), .Dimnames = list(
                c("G1", "G2", "G3", "G4"), "parents")), error.rates = list(
                error.fp = NA, error.fn = NA), adj.matrix = list(
                adj.matrix.pf = structure(c(0, 0, 0, 0, 0, 0, 
                0, 0, 1, 0, 0, 0, 0, 0, 0, 0), .Dim = c(4L, 4L
                ), .Dimnames = list(c("G1", "G2", "G3", "G4"), 
                  c("G1", "G2", "G3", "G4"))), adj.matrix.fit = structure(c(0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), .Dim = c(4L, 
                4L), .Dimnames = list(c("G1", "G2", "G3", "G4"
                ), c("G1", "G2", "G3", "G4")))), score = -114.785248346005, 
            logLik = -110.785248346005)), parameters = list(algorithm = "CAPRI", 
            command = "hc", regularization = "aic", do.boot = TRUE, 
            nboot = 100, pvalue = 0.05, min.boot = 3, min.stat = TRUE, 
            boot.seed = NULL, silent = TRUE, error.rates = list(
                epos = 0, eneg = 0), restart = 100), execution.time = structure(c(user.self = 0.0769999999999982, 
        sys.self = 0, elapsed = 0.0769999999993161, user.child = 0, 
        sys.child = 0), class = "proc_time"))), CBN_ot = list(
    edges = structure(list(From = c("Root", "C", "D", "A", "Root"
    ), To = c("A", "B", "B", "C", "D"), edge = c("Root -> A", 
    "C -> B", "D -> B", "A -> C", "Root -> D"), init_lambda = c(1.414067, 
    0.009031, 0.009031, 0.012027, 0.381915), final_lambda = c(1.416377, 
    0.109282, 0.109282, 0.014903, 0.380714), rerun_lambda = c(1.416373, 
    0.110483, 0.110483, 0.01489, 0.380711), CBN_edgeBootFreq = c(NA, 
    NA, NA, NA, NA)), class = "data.frame", row.names = c(NA, 
    -5L)), nboot = 0, init.poset = "OT"), MCCBN = list(edges = structure(list(
    From = c("Root", "Root", "Root", "A"), To = c("A", "B", "D", 
    "C"), edge = c("Root -> A", "Root -> B", "Root -> D", "A -> C"
    ), lambda = c(1.17159022808865, 0.217522208241225, 0.460587464855744, 
    0.25995329754024)), row.names = c(NA, -4L), class = "data.frame")), 
    time_ot = c(elapsed = 0.00399999999899592), time_caprese = c(elapsed = 0.011000000000422), 
    time_capri_bic = c(elapsed = 0.0860000000011496), time_capri_aic = c(elapsed = 0.0869999999995343), 
    time_cbn_ot = c(elapsed = 0.199999999998909), time_mccbn = c(elapsed = 0.514999999999418), 
    input_data = structure(c(0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 0L, 
    0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 
    0L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 
    1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 
    0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 
    1L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 
    0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L), .Dim = c(50L, 
    4L), .Dimnames = list(c("1", "2", "3", "4", "5", "6", "7", 
    "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", 
    "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", 
    "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", 
    "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", 
    "48", "49", "50"), c("A", "B", "C", "D"))), input_data_pseudosamples = structure(c(0L, 
    0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 
    0L, 0L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 
    1L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 
    0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 
    0L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 
    0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 
    1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 
    0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 
    0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 
    0L, 1L, 0L, 0L), .Dim = c(50L, 4L), .Dimnames = list(c("1", 
    "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", 
    "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", 
    "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", 
    "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", 
    "43", "44", "45", "46", "47", "48", "49", "50"), c("A", "B", 
    "C", "D"))))
                                                                                                                                        
  
  mm <- c("OT", "CAPRESE", "CAPRI_BIC", "CAPRI_AIC",
          "CBN_ot", "MCCBN")
  
  mm0  <- lapply(cpm_out_others1[mm], cpm_access_genots_paths_w)
  mmSM <- lapply(cpm_out_others1[mm], cpm_access_genots_paths_w_simplified)
  
  ## Identical unweighted transition matrices (fitness graphs)
  uw0 <- lapply(mm0, function(x) rowScaleMatrix(x$fgraph))
  uwSM <- lapply(mmSM, function(x) rowScaleMatrix(x$fgraph))
  uwSMm <- lapply(uwSM, as.matrix)
  expect_identical(uw0, uwSMm)
  
  ## Identical Weighted transition matrices
  wg0 <- lapply(mm0[c("OT", "MCCBN", "CBN_ot")],
                function(x) x$trans_mat_genots)
  wgSM <- lapply(mmSM[c("OT", "MCCBN", "CBN_ot")],
                 function(x) x$trans_mat_genots)
  wgSMm <- lapply(wgSM, as.matrix)
  expect_identical(wg0, wgSMm)
  
  ## Diagonal
  td0 <- lapply(mm0[c("MCCBN", "CBN_ot")],
                function(x)
                  trans_rate_to_trans_mat(x$weighted_fgraph,
                                          method = "uniformization")) 
  tdSM <- lapply(mmSM[c("MCCBN", "CBN_ot")],
                 function(x)
                   trans_rate_to_trans_mat(x$weighted_fgraph,
                                           method = "uniformization")) 
  tdSMm <- lapply(tdSM, as.matrix)
  expect_identical(td0, tdSMm)
  
  ## recheck competing exponentials done differently
  wg2 <- lapply(mm0[c("OT", "MCCBN", "CBN_ot")],
                function(x) trans_rate_to_trans_mat(x$weighted_fgraph,
                                                    method = "competingExponentials"))
  wg2SM <- lapply(mmSM[c("OT", "MCCBN", "CBN_ot")],
                  function(x) as.matrix(trans_rate_to_trans_mat(x$weighted_fgraph,
                                                                method = "competingExponentials")))
  expect_identical(wg2, wg2SM)
  ## tripwire to check tests run
  ## expect_identical(1, 3)
  
  rm(mm, mm0, mmSM)
  rm(cpm_out_others1)
  
  test_others <- function(data) {
    data <- as.matrix(data)
    data <- df_2_mat_integer(data)
    cpm_out_others2 <- all_methods(data, do_MCCBN = MCCBN_INSTALLED)
    
    mm <- c("OT",
            "CAPRESE", "CAPRI_BIC", "CAPRI_AIC",
            "CBN_ot"
            , "MCCBN")[c(TRUE, rep(FALSE, 3), TRUE, MCCBN_INSTALLED)]
    
    mm0  <- lapply(cpm_out_others2[mm], cpm_access_genots_paths_w)
    mmSM <- lapply(cpm_out_others2[mm], cpm_access_genots_paths_w_simplified)
    
    ## Identical unweighted transition matrices (fitness graphs)
    uw0 <- lapply(mm0, function(x) rowScaleMatrix(x$fgraph))
    uwSM <- lapply(mmSM, function(x) rowScaleMatrix(x$fgraph))
    uwSMm <- lapply(uwSM, as.matrix)
    expect_identical(uw0, uwSMm)
    
    ## Identical Weighted transition matrices
    wg0 <- lapply(mm0[c("OT" ,
                        "MCCBN"
                        , "CBN_ot")[c(TRUE, MCCBN_INSTALLED, TRUE)]],
                  function(x) x$trans_mat_genots)
    wgSM <- lapply(mmSM[c("OT" ,
                          "MCCBN"
                          , "CBN_ot")[c(TRUE, MCCBN_INSTALLED, TRUE)]],
                   function(x) x$trans_mat_genots)
    wgSMm <- lapply(wgSM, as.matrix)
    expect_equal(wg0, wgSMm)
    
    ## Diagonal
    td0 <- lapply(mm0[c("MCCBN",
                        "CBN_ot")[c(MCCBN_INSTALLED, TRUE)]],
                  function(x)
                    trans_rate_to_trans_mat(x$weighted_fgraph,
                                            method = "uniformization")) 
    tdSM <- lapply(mmSM[c("MCCBN",
                          "CBN_ot")[c(MCCBN_INSTALLED, TRUE)]],
                   function(x)
                     trans_rate_to_trans_mat(x$weighted_fgraph,
                                             method = "uniformization")) 
    tdSMm <- lapply(tdSM, as.matrix)
    expect_equal(td0, tdSMm)
    
    ## recheck competing exponentials done differently
    wg2 <- lapply(mm0[c("OT", "MCCBN"
                        , "CBN_ot")[c(TRUE, MCCBN_INSTALLED, TRUE)]],
                  function(x) trans_rate_to_trans_mat(x$weighted_fgraph,
                                                      method = "competingExponentials"))
    wg2SM <- lapply(mmSM[c("OT", "MCCBN"
                           , "CBN_ot")[c(TRUE, MCCBN_INSTALLED, TRUE)]],
                    function(x) as.matrix(trans_rate_to_trans_mat(x$weighted_fgraph,
                                                                  method = "competingExponentials")))
    expect_equal(wg2, wg2SM)
    
  }
  
  ## Additional testing
  
  
  Dat1 <- readRDS(file="./MHN/data/BreastCancer.rds") [1:50, 1:4]
  test_others(Dat1)
  
  
  Dat1 <- readRDS(file="./MHN/data/ColorectalCancer.rds")[1:40, 1:6]
  test_others(Dat1)
  
  Dat1 <- readRDS(file="./MHN/data/RenalCellCarcinoma.rds")[1:30,2:6]
  test_others(Dat1)
  
  Dat1 <- readRDS(file="./MHN/data/Glioblastoma.rds")[1:20, 1:5]
  test_others(Dat1)
  
  rm(Dat1)
  
}) ## end of testing