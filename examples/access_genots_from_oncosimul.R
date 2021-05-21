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


## From email 20-May-2021, at 17:32:51
ex_pmce_out1 <- read.table("ex_pmce_out1.txt", header = TRUE)


## fitness, target min fitness, target max fitness -> rescaled fitness
scale_fitness <- function(x, min_f, max_f) {
    min_x <- min(x)
    max_x <- max(x)
    return(min_f + (x - min_x) * (max_f - min_f)/(max_x - min_x))
}

## fitness, target max fitness. WT fitness always 1.
scale_fitness_2 <- function(x, max_f) {
    max_x <- max(x)
    return(1.0 +  (x - 1) * ( (max_f - 1)/(max_x - 1) ))
}


## output from CPMs, fitness scale -> input for OncoSimulR evalAllGenotypes
##   min_f: must pass this argument: minimum fitness
cpm_out_to_oncosimul <- function(x) {
    sh <- -Inf
    
    if("rerun_lambda" %in% names(x)) { ## CBN
        s <- exp(x$rerun_lambda) - 1
        typeDep <- "AND"
    } else if("Relation" %in% names(x)) { ## PMCE
        s <- exp(x$Lambda) - 1        
        typeDep <- x$Relation
        typeDep[typeDep == "Single"] <- "AND"
    } else if("OT_edgeWeight" %in% names(x) ) { ## OT
        s <- exp(x$OT_edgeWeight) - 1                
        typeDep <- "AND"
    } else if("whatever" %in% names(x)) { ## Something for DB
        
    } else if("otro" %in% names(x) ) { ## MCCBN
        
    } else {
        stop("Input not recognized")
    }
    x1 <- data.frame(parent = x$From,
                     child  = x$To,
                     s = s,
                     sh = sh,
                     typeDep = typeDep
                     )
    return(x1)
}


## output from CPMs, min and max final fitness, max num genots -> fitness of all genotypes
##    if max_f is NULL: no rescaling of fitness
##       max_f and min_f should have no effect in probs transition
##    max_genots: argument max of evalAllGenotypes
cpm_to_fitness_genots <- function(x, max_f = 3, max_genots = 2^15) {
    x1 <- cpm_out_to_oncosimul(x)
    x1 <- evalAllGenotypes(fitnessEffects = allFitnessEffects(rT = x1),
                           addwt = TRUE, max = max_genots)
    x1$Fitness[x1$Fitness > 0.0] <- log(x1$Fitness[x1$Fitness > 0.0]) + 1
    if(!is.null(max_f)) {
        if(max_f < 1) stop("max_f must be larger than min_f")
        x1$Fitness[x1$Fitness > 0.0] <-
            scale_fitness_2(x1$Fitness[x1$Fitness > 0.0], max_f)
    }
    return(x1)
}


## output from CPMs, min and max final fitness, max num genots ->
##                            list(accessible genotypes ,
##                                 fitness graph,
##                                 transition matrix between genotypes)
##    max_f and min_f should have no effect in probs transition
##    accessible_genotypes: gives the genotypes and their fitness,
##                          computed from the lambdas (thetas) and scaled
cpm_to_trans_mat_oncosimul <- function(x, max_f = 3,
                                       max_genots = 2^15) {
    fitness_all <- cpm_to_fitness_genots(x, max_f = max_f,
                                         max_genots = max_genots)

    access_genots_fitness <- fitness_all$Fitness[fitness_all$Fitness > 1.0]
    names(access_genots_fitness) <- fitness_all$Genotype[fitness_all$Fitness > 1.0]
    
    ## Silly? Inside the next, we now put them together. FIXME?
    access_genots_as_list <- lapply(names(access_genots_fitness),
                                    function(v) strsplit(v, ", ")[[1]])

    fgraph <- unrestricted_fitness_graph_sparseM(access_genots_as_list)
    
    ## Not very efficient, but works
    o_access_genots_fitness <- c("WT" = 1, access_genots_fitness)[colnames(fgraph)]
    mf <- matrix(rep(o_access_genots_fitness, nrow(fgraph)),
                 nrow = nrow(fgraph), byrow = TRUE)
    ## As difference w.r.t. original genotype
    mf <- mf - o_access_genots_fitness
    tm <- fgraph * mf
    tm_s <- tm / rowSums(tm)
    tm_s[nrow(tm_s), ] <- 0
    return(list(accessible_genotypes = access_genots_fitness,
                fitness_graph = fgraph,
                lambdas = tm, ## for checking
                transition_matrix = tm_s
                ))
}
## shorter
cpm2tm <- cpm_to_trans_mat_oncosimul

## Same transitions, different fitness
cpm2tm(ex_cbn_out2, max_f = NULL)$transition_matrix
cpm2tm(ex_cbn_out2, 8)$transition_matrix

## Original code
cpm_access_genots_paths_w_simplified(list(edges = ex_cbn_out2))$trans_mat_genots


library(testthat)

## For testing
reorder_trans_mat <- function(x) {
    gg <- c(1, 1 + order(colnames(x)[-1]))
    return(x[gg, gg])
}

expect_equal(reorder_trans_mat(cpm2tm(ex_cbn_out1, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_cbn_out1))$trans_mat_genots),
             check.attributes = FALSE)

expect_equal(reorder_trans_mat(cpm2tm(ex_cbn_out2, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_cbn_out2))$trans_mat_genots),
             check.attributes = FALSE)

expect_equal(reorder_trans_mat(cpm2tm(ex_cbn_out3, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_cbn_out3))$trans_mat_genots),
             check.attributes = FALSE)

expect_equal(reorder_trans_mat(cpm2tm(ex_cbn_out4, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_cbn_out4))$trans_mat_genots),
             check.attributes = FALSE)

expect_equal(reorder_trans_mat(cpm2tm(ex_cbn_out5, max_f = NULL)$transition_matrix),
             reorder_trans_mat(cpm_access_genots_paths_w_simplified(
                 list(edges = ex_cbn_out5))$trans_mat_genots),
             check.attributes = FALSE)


expect_equal(cpm2tm(ex_cbn_out2, max_f = NULL)$transition_matrix,
                  cpm2tm(ex_cbn_out2, max_f = 2)$transition_matrix)
expect_equal(max(cpm2tm(ex_cbn_out2, max_f = 8)$accessible_genotypes),
             8)
expect_equal(max(cpm2tm(ex_cbn_out2, max_f = 3)$accessible_genotypes),
             3)

expect_equal(cpm2tm(ex_cbn_out1, max_f = NULL)$transition_matrix,
             cpm2tm(ex_cbn_out1, max_f = 2)$transition_matrix)


expect_equal(cpm2tm(ex_cbn_out4, max_f = NULL)$transition_matrix,
             cpm2tm(ex_cbn_out4, max_f = 2)$transition_matrix)

expect_equal(cpm2tm(ex_cbn_out4, max_f = NULL)$transition_matrix,
             cpm2tm(ex_cbn_out4, max_f = 4)$transition_matrix)

expect_equal(cpm2tm(ex_cbn_out5, max_f = NULL)$transition_matrix,
             cpm2tm(ex_cbn_out5, max_f = 2)$transition_matrix)







## Running internal functions

cpm_out_to_oncosimul(ex_pmce_out1)
cpm_out_to_oncosimul(ex_cbn_out2)
cpm_out_to_oncosimul(ex_cbn_out2)

cpm_to_fitness_genots(ex_cbn_out1)
cpm_to_fitness_genots(ex_cbn_out2)
cpm_to_fitness_genots(ex_ot_out2)
cpm_to_fitness_genots(ex_pmce_out1)

cpm_to_fitness_genots(ex_cbn_out2, max_f = NULL)
cpm_to_fitness_genots(ex_cbn_out2, max_f = 3)

cpm2tm(ex_cbn_out2, max_f = NULL)$lambdas
cpm2tm(ex_cbn_out2, 1.01, 8)$lambdas













## ## output from CPMs, fitness scale -> input for OncoSimulR evalAllGenotypes
## ##   min_f: must pass this argument: minimum fitness
## cpm_out_to_oncosimul <- function(x, min_f = 1.01) {
##     sh <- -Inf
    
##     if("rerun_lambda" %in% names(x)) { ## CBN
##         s <- x$rerun_lambda * (min_f/min(x$rerun_lambda)) - 1
##         typeDep <- "AND"
##     } else if("Relation" %in% names(x)) { ## PMCE
##         s <- x$Lambda * (min_f/min(x$Lambda)) - 1
##         typeDep <- x$Relation
##         typeDep[typeDep == "Single"] <- "AND"
##     } else if("OT_edgeWeight" %in% names(x) ) { ## OT
##         s <- x$OT_edgeWeight * (min_f/min(x$OT_edgeWeight)) - 1
##         typeDep <- "AND"
##     } else if("whatever" %in% names(x)) { ## Something for DB
        
##     } else if("otro" %in% names(x) ) { ## MCCBN
        
##     } else {
##         stop("Input not recognized")
##     }
##     x1 <- data.frame(parent = x$From,
##                      child  = x$To,
##                      s = s,
##                      sh = sh,
##                      typeDep = typeDep
##                      )
##     return(x1)
## }


## ## output from CPMs, min and max final fitness, max num genots -> fitness of all genotypes
## ##    if max_f is NULL: no rescaling of fitness
## ##       max_f and min_f should have no effect in probs transition
## ##    max_genots: argument max of evalAllGenotypes
## cpm_to_fitness_genots <- function(x, min_f = 1.1, max_f = 3, max_genots = 2^15) {
##     x1 <- cpm_out_to_oncosimul(x, min_f)
##     x1 <- evalAllGenotypes(fitnessEffects = allFitnessEffects(rT = x1),
##                           addwt = TRUE, max = max_genots)
##     if(!is.null(max_f)) {
##         if(max_f < min_f) stop("max_f must be larger than min_f")
##         x1$Fitness[x1$Fitness > 1.000] <-
##             scale_fitness(x1$Fitness[x1$Fitness > 1.000], min_f, max_f)
##     }
##     return(x1)
## }
