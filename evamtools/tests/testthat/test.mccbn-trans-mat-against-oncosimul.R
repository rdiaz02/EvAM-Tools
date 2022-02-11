## Testing that cpm2tm gives same output as
## cpm_to_trans_mat_oncosimul, the function that uses OncoSimulR

test_that("Testing evamtools:::cpm2tm by comparing with
OncoSimulR's based cpm_to_trans_mat_oncosimul.", {
    ## Recall evamtools:::cpm2F2tm <- cpm_to_trans_mat_oncosimul

    ## The first set of tests use pre-run analyses. So we would not catch
    ## changes in the code that gives the analysis itself. See below for
    ## code that does that.

    print("Running these tests may take a while") ## Maybe we want to move it somewhere else?

    ## Functions for testing
    reorder_trans_mat <- function(x) {
        gg <- c(1, 1 + order(colnames(x)[-1]))
        return(as.matrix(x[gg, gg]))
    }
    
    run_test_for_dataset1 <- function(data) {
        expect_equal(
            reorder_trans_mat(evamtools:::cpm2F2tm(data, max_f = NULL)$transition_matrix),
            reorder_trans_mat(evamtools:::cpm2tm(
                                              list(edges = data))$trans_mat_genots),
            check.attributes = TRUE)

        ## Do not run next test if fitness is absurdly large as scaling will fail
        ## as it should.
        max_fitness <- max(evamtools:::cpm2F2tm(data, max_f = NULL)$accessible_genotypes)
        if(max_fitness < 1e10) {
            maxff <- sample(c(2, 3, 3.5, 3.8, 4, 5, 8), size = 1)
            expect_equal(as.matrix(evamtools:::cpm2F2tm(data, max_f = NULL)$transition_matrix),
                         as.matrix(evamtools:::cpm2F2tm(data, max_f = maxff)$transition_matrix))
        }
    }

    ## this takes the object, not the object$edges
    run_test_for_dataset <- function(data){
        expect_equal(reorder_trans_mat(evamtools:::cpm2F2tm(data$edges, max_f = NULL)$transition_matrix),
                    reorder_trans_mat(evamtools:::cpm2tm(
                        data)$trans_mat_genots),
                    check.attributes = TRUE)
        maxff <- sample(c(2, 3, 3.5, 3.8, 4, 5, 8), size = 1)
        expect_equal(as.matrix(evamtools:::cpm2F2tm(data$edges, max_f = NULL)$transition_matrix),
                    as.matrix(evamtools:::cpm2F2tm(data$edges, max_f = maxff)$transition_matrix))
    }

    ## First, load a bunch of data structures
    ## From ex1.R
    ex_mccbn_out1 <- structure(list(From = c("Root", "Root"),
                                  To = c("A", "B"),
                                  edge = c("Root -> A", "Root -> B"),
                                  lambda = c(88.234297, 268.921382)),
                             class = "data.frame", row.names = c("A", "B")) 

    ex_mccbn_out2 <- structure(list(From = c("Root", "A", "Root", "C"),
                                  To = c("A", "B", "C", "D"),
                                  edge = c("Root -> A", "A -> B", "Root -> C", "C -> D"),
                                  init_lambda = c(0.689845, 2.734304, 0.6988, 2.463583),
                                  final_lambda = c(0.689845, 2.734304, 0.6988, 2.463583),
                                  lambda = c(0.689845, 2.734304, 0.6988, 2.463583)),
                             class = "data.frame", row.names = c("A", "B", "C", "D"))

    ## From simulating-posets.R
    ex_mccbn_out3 <- structure(list(From = c("Root", "A", "B", "C", "Root", "E"), 
                                  To = c("A", "B", "C", "D", "E", "F"),
                                  edge = c("Root -> A", 
                                           "A -> B", "B -> C", "C -> D", "Root -> E", "E -> F"),
                                  lambda = c(0.874492, 
                                                   0.920866, 0.748741, 0.678371, 0.900958, 0.842805)),
                             class = "data.frame",
                             row.names = c("A", 
                                           "B", "C", "D", "E", "F"))

    ex_mccbn_out4 <- structure(list(From = c("Root", "A", "A", "Root", "Root", "D"
                                           ), To = c("A", "B", "C", "D", "E", "F"),
                                  edge = c("Root -> A", 
                                           "A -> B", "A -> C", "Root -> D", "Root -> E", "D -> F"),
                                  lambda = c(0.762984, 
                                                   0.49971, 0.748445, 0.712476, 0.979017, 0.830669)),
                             class = "data.frame",
                             row.names = c("A", 
                                           "B", "C", "D", "E", "F"))

    ex_mccbn_out5 <- structure(list(From = c("Root", "Root", "Root", "Root", "A", 
                                           "B", "D", "D", "B"), To = c("A", "B", "C", "D", "E", "E", "E", 
                                                                       "F", "G"),
                                  edge = c("Root -> A", "Root -> B", "Root -> C", "Root -> D", 
                                           "A -> E", "B -> E", "D -> E", "D -> F", "B -> G"),
                                  lambda = c(6.34855, 
                                                   5.749031, 3.436123, 2.721351, 0.249042, 0.249042, 0.249042, 3.046078, 
                                                   2.448528)),
                             class = "data.frame",
                             row.names = c(NA, -9L))

    ex_mccbn_out6 <- structure(list(From = c("Root", "Root", "Root", "Root", "D", 
                                           "A", "E", "A"), To = c("A", "B", "C", "D", "E", "F", "F", "G"
                                                                  ),
                                  edge = c("Root -> A", "Root -> B", "Root -> C", "Root -> D", 
                                           "D -> E", "A -> F", "E -> F", "A -> G"),
                                  lambda = c(6.598142, 4.978507, 
                                                   2.821383, 1.27322, 4.961959, 2.363096, 2.363096, 4.09447)),
                             class = "data.frame",
                             row.names = c(NA, 
                                           -8L))

    ex_mccbn_out7 <- structure(list(From = c("Root", "Root", "B", "B", "B", "Root"
                                           ), To = c("A", "B", "C", "D", "E", "F"),
                                  edge = c("Root -> A", 
                                           "Root -> B", "B -> C", "B -> D", "B -> E", "Root -> F"),
                                  lambda = c(4.107311, 
                                                   2.931829, 5.312215, 1.365132, 5.434725, 0.367228)),
                             class = "data.frame",
                             row.names = c("A", 
                                           "B", "C", "D", "E", "F"))

    precomputed_datasets <- list(ex_mccbn_out1, ex_mccbn_out2, ex_mccbn_out3,
                     ex_mccbn_out4, ex_mccbn_out5, ex_mccbn_out6,
                     ex_mccbn_out7)

    for(ex in precomputed_datasets){
        run_test_for_dataset1(ex)
    }

    data(examples_csd)

    mccbn_hcbn2_opts <- list()
    mccbn_hcbn2_opts$thrds <- 4
    mccbn_hcbn2_opts$max.iter.asa <- 10 #Just to run faster   
    mccbn_hcbn2_opts$addname <- "tmpo"

    do_MCCBN_HCBN<- function(x) suppressMessages(evamtools:::mccbn_hcbn_proc(x,
                                                                mccbn_hcbn2_opts))

    print("Running do_MCCBN_HCBN: this may take a while")
    print("Completed: 0/7")
    ex_mccbn_and <- do_MCCBN_HCBN(examples_csd$csd$AND$data)
    print("Completed: 1/7")
    ex_mccbn_linear <- do_MCCBN_HCBN(examples_csd$csd$Linear$data)
    print("Completed: 2/7")
    ex_mccbn_or <- do_MCCBN_HCBN(examples_csd$csd$OR$data)
    print("Completed: 3/7")
    ex_mccbn_xor <- do_MCCBN_HCBN(examples_csd$csd$XOR$data)
    print("Completed: 4/7")
    ex_mccbn_c1 <- do_MCCBN_HCBN(examples_csd$csd$c1$data)
    print("Completed: 5/7")
    ex_mccbn_c3 <- do_MCCBN_HCBN(examples_csd$csd$c3$data)
    print("Completed: 6/7")
    ex_mccbn_c4c2 <- do_MCCBN_HCBN(examples_csd$csd$c4c2$data)
    print("Completed: 7/7")
    
    all_mccbn_examples <- list(ex_mccbn_and, ex_mccbn_linear, ex_mccbn_or, ex_mccbn_xor,
        ex_mccbn_c1, ex_mccbn_c3, ex_mccbn_c4c2)
   
    for (ex in all_mccbn_examples) {
        run_test_for_dataset(ex)
    }

  ## Add the examples used for HESBCN with mixes of relationships
    set.seed(2)
    d1 <- data.frame(A = sample(c(1, 0), prob = c(0.7, 0.2), size = 200, replace = TRUE),
                     B = sample(c(1, 0), prob = c(0.85, 0.2), size = 200, replace = TRUE))
    d1$C <- 0
    d1$D <- 0
    d1$E <- 0
    d1$F <- 0
    d1$C[(d1$A == 1) & (d1$B == 1)] <- 1
    d1$D[(d1$A == 1) | (d1$B == 1)] <- 1
    d1$D[100:200] <- 0
    d1$E[xor((d1$A == 1), (d1$B == 1))] <- 1
    d1$F[(d1$C == 1) & (d1$D == 1)] <- 1
    d2 <- rbind(d1,
                data.frame(A = sample(c(1, 0), size = 25, prob = c(0.5, 0.2), replace = TRUE),
                           B = sample(c(1, 0), size = 25, prob = c(0.5, 0.2), replace = TRUE),
                           C = 0,
                           D = sample(c(1, 0), size = 25, replace = TRUE),
                           E = 0,
                           F = 0))
    d3 <- d2
    d3$C[(d3$A == 1) & (d3$B == 1)] <- 1
    
    set.seed(22)
    d111 <- data.frame(A = sample(c(1, 0), prob = c(0.7, 0.2), size = 2000, replace = TRUE),
                     B = sample(c(1, 0), prob = c(0.85, 0.2), size = 2000, replace = TRUE))
    d111$C <- 0
    d111$D <- 0
    d111$E <- 0
    d111$F <- 0
    d111$C[(d111$A == 1) & (d111$B == 1)] <- 1
    d111$D[(d111$A == 1) | (d111$B == 1)] <- 1
    d111$D[100:200] <- 0
    d111$E[xor((d111$A == 1), (d111$B == 1))] <- 1
    d111$F[(d111$C == 1) & (d111$D == 1)] <- 1
    d222 <- rbind(d111,
                data.frame(A = sample(c(1, 0), size = 250, prob = c(0.5, 0.2), replace = TRUE),
                           B = sample(c(1, 0), size = 250, prob = c(0.5, 0.2), replace = TRUE),
                           C = 0,
                           D = sample(c(1, 0), size = 25, replace = TRUE),
                           E = 0,
                           F = 0))
    d333 <- d222
    d333$C[(d333$A == 1) & (d333$B == 1)] <- 1

    set.seed(NULL)

    print("Running do_MCCBN_HCBN: this may take a while")
    print("Completed: 0/4")
    d3_1 <- do_MCCBN_HCBN(as.matrix(d3))
    print("Completed: 1/4")
    d2_1 <- do_MCCBN_HCBN(as.matrix(d2))
    print("Completed: 2/4")
    d222_1 <- do_MCCBN_HCBN(as.matrix(d222))
    print("Completed: 3/4")
    d333_1 <- do_MCCBN_HCBN(as.matrix(d333))
    print("Completed: 4/4")
    
    all_mixed_examples <- list(d3_1, d2_1, d222_1, d333_1)
    
    for(ex in all_mixed_examples){
        run_test_for_dataset(ex)
    }
    set.seed(NULL)

})

cat("\n Done test.mccbn-trans-mat-against-oncosimul.R \n")

