local({
    ## From ex1.R
ex_cbn_out1 <- structure(list(From = c("Root", "Root"),
                              To = c("A", "B"), edge = c("Root -> A", "Root -> B"),
                              init_lambda = c(88.234297, 268.921382),
                              final_lambda = c(88.234297, 268.921382),
                              rerun_lambda = c(88.234297, 268.921382),
                              CBN_edgeBootFreq = c(NA, NA)),
                         class = "data.frame", row.names = c("A", "B")) 

ex_cbn_out2 <- structure(list(From = c("Root", "A", "Root", "C"),
                              To = c("A", "B", "C", "D"),
                              edge = c("Root -> A", "A -> B", "Root -> C", "C -> D"),
                              init_lambda = c(0.689845, 2.734304, 0.6988, 2.463583),
                              final_lambda = c(0.689845, 2.734304, 0.6988, 2.463583),
                              rerun_lambda = c(0.689845, 2.734304, 0.6988, 2.463583),
                              CBN_edgeBootFreq = c(NA, NA, NA, NA)),
                         class = "data.frame", row.names = c("A", "B", "C", "D"))


## From simulating-posets.R

ex_cbn_out3 <- structure(list(From = c("Root", "A", "B", "C", "Root", "E"), 
    To = c("A", "B", "C", "D", "E", "F"), edge = c("Root -> A", 
    "A -> B", "B -> C", "C -> D", "Root -> E", "E -> F"), init_lambda = c(0.874492, 
    0.920866, 0.748741, 0.678371, 0.900958, 0.842805), final_lambda = c(0.874492, 
    0.920866, 0.748741, 0.678371, 0.900958, 0.842805), rerun_lambda = c(0.874492, 
    0.920866, 0.748741, 0.678371, 0.900958, 0.842805), CBN_edgeBootFreq = c(NA, 
    NA, NA, NA, NA, NA)), class = "data.frame", row.names = c("A", 
                                                              "B", "C", "D", "E", "F"))

ex_cbn_out4 <- structure(list(From = c("Root", "A", "A", "Root", "Root", "D"
), To = c("A", "B", "C", "D", "E", "F"), edge = c("Root -> A", 
"A -> B", "A -> C", "Root -> D", "Root -> E", "D -> F"), init_lambda = c(1.574403, 
1.477307, 1.129679, 1.402693, 1.282986, 1.153868), final_lambda = c(0.763003, 
0.499678, 0.748408, 0.712498, 0.97903, 0.830571), rerun_lambda = c(0.762984, 
0.49971, 0.748445, 0.712476, 0.979017, 0.830669), CBN_edgeBootFreq = c(NA, 
NA, NA, NA, NA, NA)), class = "data.frame", row.names = c("A", 
                                                          "B", "C", "D", "E", "F"))

ex_cbn_out5 <- structure(list(From = c("Root", "Root", "Root", "Root", "A", 
"B", "D", "D", "B"), To = c("A", "B", "C", "D", "E", "E", "E", 
"F", "G"), edge = c("Root -> A", "Root -> B", "Root -> C", "Root -> D", 
"A -> E", "B -> E", "D -> E", "D -> F", "B -> G"), init_lambda = c(7.567088, 
8.029489, 6.199193, 3.9428, 0.276702, 0.276702, 0.276702, 3.772443, 
2.231341), final_lambda = c(6.352607, 5.752454, 3.437221, 2.72262, 
0.248693, 0.248693, 0.248693, 3.045175, 2.448479), rerun_lambda = c(6.34855, 
5.749031, 3.436123, 2.721351, 0.249042, 0.249042, 0.249042, 3.046078, 
2.448528), CBN_edgeBootFreq = c(NA, NA, NA, NA, NA, NA, NA, NA, 
NA)), class = "data.frame", row.names = c(NA, -9L))

ex_cbn_out6 <- structure(list(From = c("Root", "Root", "Root", "Root", "D", 
"A", "E", "A"), To = c("A", "B", "C", "D", "E", "F", "F", "G"
), edge = c("Root -> A", "Root -> B", "Root -> C", "Root -> D", 
"D -> E", "A -> F", "E -> F", "A -> G"), init_lambda = c(9.394563, 
15.91677, 7.428318, 2.40452, 8.45151, 2.603337, 2.603337, 4.274681
), final_lambda = c(6.59639, 4.978562, 2.82143, 1.273479, 4.952848, 
2.36438, 2.36438, 4.095935), rerun_lambda = c(6.598142, 4.978507, 
2.821383, 1.27322, 4.961959, 2.363096, 2.363096, 4.09447), CBN_edgeBootFreq = c(NA, 
NA, NA, NA, NA, NA, NA, NA)), class = "data.frame", row.names = c(NA, 
-8L))



ex_cbn_out7 <- structure(list(From = c("Root", "Root", "B", "B", "B", "Root"
), To = c("A", "B", "C", "D", "E", "F"), edge = c("Root -> A", 
"Root -> B", "B -> C", "B -> D", "B -> E", "Root -> F"), init_lambda = c(5.551964, 
5.782397, 20.123668, 2.31525, 12.844557, 0.844419), final_lambda = c(4.106536, 
2.932381, 5.310502, 1.364737, 5.431533, 0.367361), rerun_lambda = c(4.107311, 
2.931829, 5.312215, 1.365132, 5.434725, 0.367228), CBN_edgeBootFreq = c(NA, 
NA, NA, NA, NA, NA)), class = "data.frame", row.names = c("A", 
"B", "C", "D", "E", "F"))


ex_ot_out1 <-structure(list(From = c("Root", "Root", "B", "B", "F", "D", "B"
), To = c("A", "B", "C", "D", "E", "F", "G"), edge = c("Root -> A", 
"Root -> B", "B -> C", "B -> D", "F -> E", "D -> F", "B -> G"
), OT_edgeBootFreq = c(NA, NA, NA, NA, NA, NA, NA), OT_edgeWeight = c(0.905524143115485, 
0.885838835656453, 0.888232094027126, 0.849455249540266, 0.235215622554582, 
0.777051610224958, 0.740395374420974), OT_obsMarginal = c(0.874, 
0.855, 0.779, 0.74, 0.132, 0.555, 0.611), OT_predMarginal = c(0.874, 
0.855, 0.759438440393193, 0.726284238356928, 0.132746367988167, 
0.564360336896258, 0.633038045129933)), class = "data.frame", row.names = c("A", 
"B", "C", "D", "E", "F", "G"))

    
ex_ot_out2 <- structure(list(From = c("Root", "A", "Root", "C"),
                             To = c("A", "B", "C", "D"),
                             edge = c("Root -> A", "A -> B", "Root -> C", 
                                      "C -> D"),
                             OT_edgeBootFreq = c(NA, NA, NA, NA),
                             OT_edgeWeight = c(0.490601503759398, 0.616858237547893,
                                               0.494360902255639, 0.593155893536122),
                             OT_obsMarginal = c(0.490601503759398, 0.302631578947368,
                                                0.494360902255639, 0.293233082706767),
                             OT_predMarginal = c(0.490601503759398, 0.302631578947368,
                                                 0.494360902255639, 0.293233082706767)),
                        class = "data.frame",
                        row.names = c("A", "B", "C", "D"))

ex_ot_out3 <- structure(list(From = c("Root", "A", "B", "G", "D", "E", "A"), 
    To = c("A", "B", "C", "D", "E", "F", "G"), edge = c("Root -> A", 
    "A -> B", "B -> C", "G -> D", "D -> E", "E -> F", "A -> G"
    ), OT_edgeBootFreq = c(NA, NA, NA, NA, NA, NA, NA), OT_edgeWeight = c(0.911885569159806, 
    0.944051084043843, 0.873782587428252, 0.737946596615086, 
    0.87265719384507, 0.73043186548866, 0.846015009931598), OT_obsMarginal = c(0.868, 
    0.832, 0.74, 0.561, 0.466, 0.324, 0.699), OT_predMarginal = c(0.868, 
    0.819436340950056, 0.716009206228079, 0.541904462825413, 
    0.472896827861345, 0.345418912158432, 0.734341028620627)), class = "data.frame", row.names = c("A", 
                                                                                                   "B", "C", "D", "E", "F", "G"))

ex_ot_out4 <- structure(list(From = c("Root", "A", "E", "C", "B", "D"), To = c("A", 
"B", "C", "D", "E", "F"), edge = c("Root -> A", "A -> B", "E -> C", 
"C -> D", "B -> E", "D -> F"), OT_edgeBootFreq = c(NA, NA, NA, 
NA, NA, NA), OT_edgeWeight = c(0.798846911181523, 0.864081974864568, 
0.981087366597147, 0.674538658440992, 0.881678244444341, 0.430247912741961
), OT_obsMarginal = c(0.803, 0.745, 0.627, 0.428, 0.629, 0.268
), OT_predMarginal = c(0.803, 0.706952533306802, 0.624522271947019, 
0.452620235016971, 0.634704104103496, 0.249630357289804)), class = "data.frame", row.names = c("A", 
                                                                                               "B", "C", "D", "E", "F"))

ex_ot_out5 <- structure(list(From = c("Root", "E", "A", "A", "D", "A"), To = c("A", 
"B", "C", "D", "E", "F"), edge = c("Root -> A", "E -> B", "A -> C", 
"A -> D", "D -> E", "A -> F"), OT_edgeBootFreq = c(NA, NA, NA, 
NA, NA, NA), OT_edgeWeight = c(0.838483399072279, 0.317539113591577, 
0.915661053629564, 0.827342016912557, 0.802580929765542, 0.777992597809607
), OT_obsMarginal = c(0.8, 0.198, 0.771, 0.629, 0.52, 0.636), 
    OT_predMarginal = c(0.8, 0.181813790504971, 0.733932592469118, 
    0.664747357857712, 0.536799423955566, 0.62608920141162)), class = "data.frame", row.names = c("A", 
"B", "C", "D", "E", "F"))



cpm2tm(ex_ot_out2, max_f = NULL)$transition_matrix
cpm_access_genots_paths_w_simplified(list(edges = ex_ot_out2))$trans_mat_genots



## Same transitions, different fitness
cpm2tm(ex_cbn_out2, max_f = NULL)$transition_matrix
cpm2tm(ex_cbn_out2, 8)$transition_matrix

## Original code
cpm_access_genots_paths_w_simplified(list(edges = ex_cbn_out2))$trans_mat_genots


library(testthat)

## For testing
reorder_trans_mat <- function(x) {
    gg <- c(1, 1 + order(colnames(x)[-1]))
    return(as.matrix(x[gg, gg]))
}

## CBN
expect_equal(reorder_trans_mat(cpm2tm(ex_cbn_out1, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_cbn_out1))$trans_mat_genots),
             check.attributes = TRUE)

expect_equal(reorder_trans_mat(cpm2tm(ex_cbn_out2, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_cbn_out2))$trans_mat_genots),
             check.attributes = TRUE)

expect_equal(reorder_trans_mat(cpm2tm(ex_cbn_out3, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_cbn_out3))$trans_mat_genots),
             check.attributes = TRUE)

expect_equal(reorder_trans_mat(cpm2tm(ex_cbn_out4, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_cbn_out4))$trans_mat_genots),
             check.attributes = TRUE)

expect_equal(reorder_trans_mat(cpm2tm(ex_cbn_out5, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_cbn_out5))$trans_mat_genots),
             check.attributes = TRUE)

expect_equal(reorder_trans_mat(cpm2tm(ex_cbn_out6, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_cbn_out6))$trans_mat_genots),
             check.attributes = TRUE)


expect_equal(reorder_trans_mat(cpm2tm(ex_cbn_out7, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_cbn_out7))$trans_mat_genots),
             check.attributes = TRUE)

          
expect_equal(as.matrix(cpm2tm(ex_cbn_out2, max_f = NULL)$transition_matrix),
             as.matrix(cpm2tm(ex_cbn_out2, max_f = 2)$transition_matrix))
expect_equal(as.matrix(cpm2tm(ex_cbn_out1, max_f = NULL)$transition_matrix),
             as.matrix(cpm2tm(ex_cbn_out1, max_f = 2)$transition_matrix))
expect_equal(as.matrix(cpm2tm(ex_cbn_out4, max_f = NULL)$transition_matrix),
             as.matrix(cpm2tm(ex_cbn_out4, max_f = 2)$transition_matrix))
expect_equal(as.matrix(cpm2tm(ex_cbn_out4, max_f = NULL)$transition_matrix),
             as.matrix(cpm2tm(ex_cbn_out4, max_f = 4)$transition_matrix))
expect_equal(as.matrix(cpm2tm(ex_cbn_out5, max_f = NULL)$transition_matrix),
             as.matrix(cpm2tm(ex_cbn_out5, max_f = 2)$transition_matrix))
expect_equal(as.matrix(cpm2tm(ex_cbn_out6, max_f = NULL)$transition_matrix),
             as.matrix(cpm2tm(ex_cbn_out6, max_f = 2)$transition_matrix))
expect_equal(as.matrix(cpm2tm(ex_cbn_out7, max_f = NULL)$transition_matrix),
             as.matrix(cpm2tm(ex_cbn_out7, max_f = 4)$transition_matrix))

expect_equal(max(cpm2tm(ex_cbn_out2, max_f = 8)$accessible_genotypes),
             8)
expect_equal(max(cpm2tm(ex_cbn_out2, max_f = 3)$accessible_genotypes),
             3)


## OT 
expect_equal(reorder_trans_mat(cpm2tm(ex_ot_out1, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_ot_out1))$trans_mat_genots),
             check.attributes = TRUE)
expect_equal(reorder_trans_mat(cpm2tm(ex_ot_out2, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_ot_out2))$trans_mat_genots),
             check.attributes = TRUE)
expect_equal(reorder_trans_mat(cpm2tm(ex_ot_out3, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_ot_out3))$trans_mat_genots),
             check.attributes = TRUE)
expect_equal(reorder_trans_mat(cpm2tm(ex_ot_out4, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_ot_out4))$trans_mat_genots),
             check.attributes = TRUE)
expect_equal(reorder_trans_mat(cpm2tm(ex_ot_out5, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_ot_out5))$trans_mat_genots),
             check.attributes = TRUE)

## tripwire. Should fail
## expect_equal(as.matrix(cpm2tm(ex_ot_out1, max_f = NULL)$transition_matrix),
##              reorder_trans_mat(cpm_access_genots_paths_w_simplified(
##                  list(edges = ex_ot_out1))$trans_mat_genots),
##              check.attributes = FALSE)


expect_equal(as.matrix(cpm2tm(ex_ot_out2, max_f = NULL)$transition_matrix),
             as.matrix(cpm2tm(ex_ot_out2, max_f = 2)$transition_matrix))
expect_equal(as.matrix(cpm2tm(ex_ot_out1, max_f = NULL)$transition_matrix),
             as.matrix(cpm2tm(ex_ot_out1, max_f = 2)$transition_matrix))
expect_equal(as.matrix(cpm2tm(ex_ot_out4, max_f = NULL)$transition_matrix),
             as.matrix(cpm2tm(ex_ot_out4, max_f = 2)$transition_matrix))
expect_equal(as.matrix(cpm2tm(ex_ot_out4, max_f = NULL)$transition_matrix),
             as.matrix(cpm2tm(ex_ot_out4, max_f = 4)$transition_matrix))
expect_equal(as.matrix(cpm2tm(ex_ot_out5, max_f = NULL)$transition_matrix),
             as.matrix(cpm2tm(ex_ot_out5, max_f = 2)$transition_matrix))
expect_equal(as.matrix(cpm2tm(ex_ot_out3, max_f = NULL)$transition_matrix),
             as.matrix(cpm2tm(ex_ot_out3, max_f = 2)$transition_matrix))

expect_equal(max(cpm2tm(ex_ot_out2, max_f = 8)$accessible_genotypes),
             8)
expect_equal(max(cpm2tm(ex_ot_out2, max_f = 3)$accessible_genotypes),
             3)



## A minimal set of tests for genots_2_fgraph_and_trans_mat
## under general fitness landscapes

x1 <- c(WT = 1, A = 2.5, B = 1.5,
        "A, B" = 2, D = 4,
        "A, B, C" = 3, "A, B, C, D" = 3.5)

x1o <- genots_2_fgraph_and_trans_mat(x1)

expect_equal(unname(x1o$transition_matrix)[1, ],
             c(0, 1.5/(1.5 + .5 + 3), .5/(1.5 + .5 + 3), 3/(1.5 + .5 + 3), 0, 0, 0)
)
expect_equal(unname(x1o$transition_matrix["D", ]), rep(0, 7))
expect_equal(unname(x1o$transition_matrix["A", ]), rep(0, 7))

expect_equal(x1o$accessible_genotypes, c("A", "B", "D", "A, B", "A, B, C", "A, B, C, D"))


x2 <- c(WT = 2, A = 2.5, B = 1.5,
        "A, B" = 2, D = 4,
        "A, B, C" = 3, "A, B, C, D" = 3.5)

x2o <- suppressMessages(genots_2_fgraph_and_trans_mat(x2))

expect_equal(unname(x2o$transition_matrix)[1, ],
             c(0, 0.5/(0.5 + 2), 0, 2/(0.5 + 2), 0, 0, 0)
             )

expect_equal(unname(x2o$transition_matrix["D", ]), rep(0, 7))
expect_equal(unname(x2o$transition_matrix["A", ]), rep(0, 7))

expect_equal(x2o$accessible_genotypes, c("A", "D", "A, B", "A, B, C", "A, B, C, D"))


## tripwire
## expect_true(1 == 2)


## ## Running internal functions

## cpm_out_to_oncosimul(ex_pmce_out1)
## cpm_out_to_oncosimul(ex_cbn_out2)
## cpm_out_to_oncosimul(ex_cbn_out2)

## cpm_to_fitness_genots(ex_cbn_out1)
## cpm_to_fitness_genots(ex_cbn_out2)
## cpm_to_fitness_genots(ex_ot_out2)
## cpm_to_fitness_genots(ex_pmce_out1)

## cpm_to_fitness_genots(ex_cbn_out2, max_f = NULL)
## cpm_to_fitness_genots(ex_cbn_out2, max_f = 3)

## cpm2tm(ex_cbn_out2, max_f = NULL)$lambdas
## cpm2tm(ex_cbn_out2, 1.01, 8)$lambdas

})
