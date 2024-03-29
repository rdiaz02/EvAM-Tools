## Testing that cpm2tm gives same output as
## cpm_to_trans_mat_oncosimul, the function that uses OncoSimulR

t1 <- Sys.time()
test_that("Testing cpm2tm by comparing with
OncoSimulR's based cpm_to_trans_mat_oncosimul", {
    ## Recall cpm2F2tm <- cpm_to_trans_mat_oncosimul

    ## For testing
    reorder_trans_mat <- function(x) {
        gg <- c(1, 1 + evam_string_order(colnames(x)[-1]))
        return(as.matrix(x[gg, gg]))
    }

    run_test_for_dataset <- function(data) {
        ## Used in several places. Run just once
        cpm2_out <- cpm2F2tm(data$edges, max_f = NULL)
        if(min(data$edges$Lambda) < 1e-15) {
            warning("Skipping comparison with OncoSimul's ",
                    "transition rate matrix: ",
                    "smallest lambda < 1e-15")
            ## We compute products of numbers close to R's smallest limit.
            ## We are probably fine with numbers even as small as 2e-16, but to be safe.
        } else if (sum(data$edges$Lambda < 1e-9) >= 1) {
          warning("Skipping comparison with OncoSimul's ",
                  "numerical values of transition rate matrix: ",
                  "one or more lambdas < 1e-9.",
                  " Comparing only sets of accessible genotypes")
          cno <- colnames(reorder_trans_mat(cpm2_out$transition_matrix))
          cnd <- colnames(
              reorder_trans_mat(cpm2tm(data)$trans_mat_genots))
          expect_true(all(cno == cnd))
        } else {
            expect_equal(
            reorder_trans_mat(cpm2_out$transition_matrix),
            reorder_trans_mat(cpm2tm(
                                              data)$trans_mat_genots),
            check.attributes = TRUE)
        }

        ## Do not run next test if fitness is absurdly large as scaling will fail
        ## (as it should).  Do not run it either if ratio of largest to smallest
        ## lambda is > 1e9 as in products we accumulate tiny errors that give
        ## differences of the order of 1e-8.
        max_fitness <- max(cpm2_out$accessible_genotypes)
        ratio_lambdas <- max(data$edges$Lambdas)/min(data$edges$Lambdas)
        if((max_fitness < 1e10) && (ratio_lambdas < 1e9))  {
            maxff <- sample(c(2, 3, 3.5, 3.8, 4, 5, 8), size = 1)
            this_tolerance <- testthat_tolerance()
            if(min(data$edges$Lambdas) < 1e-5) {
                warning("Smallest lambda < 1e-5. Setting tolerance to 1e-6 ",
                        "for numerical comparison of transition rates ",
                        "between max_f = NULL and max_f = maxff")
                this_tolerance <- 1e-6
            }
            expect_equal(as.matrix(cpm2_out$transition_matrix),
                         as.matrix(cpm2F2tm(data$edges, max_f = maxff)$transition_matrix),
                         tolerance = this_tolerance
                         )
            rm(this_tolerance)
        } else {
            warning("Skipping comparison of transition matrices with ",
                    "different fitness scaling",
                    ifelse(max_fitness >= 1e10, ". max_fitness >= 1e10", ""),
                    ifelse(ratio_lambdas >= 1e9, ". Ratio of Lambdas >= 1e9", ""))
        }
    }
   
    
    ## First, load a bunch of data structures
    ## Adapted from test.OT-CBN-trans-mat-against-oncosimul.R
    ## ex_hesbcn_*: output from running evam:::do_HESBCN with examples_csd
    ## examples_csd: list of cross sectional data sets
    ## examples_csd: is located in file /data/examples_csd.RData
    ## examples_csd: is generated by script inst/miscell/toy_datasets.R

    ## Output from running do_HESBCN(examples_csd$csd$AND$data)

    data(examples_csd)
    ## Run HESBCN.
    ## The names declare intent; often, you just get "Single"
    ex_hesbcn_and <- do_HESBCN(examples_csd$csd$AND$data)
    ex_hesbcn_linear <- do_HESBCN(examples_csd$csd$Linear$data)
    ex_hesbcn_or <- do_HESBCN(examples_csd$csd$OR$data)
    ex_hesbcn_xor <- do_HESBCN(examples_csd$csd$XOR$data)
    ex_hesbcn_c1 <- do_HESBCN(examples_csd$csd$c1$data)
    ex_hesbcn_c3 <- do_HESBCN(examples_csd$csd$c3$data)
    ex_hesbcn_c4c2 <- do_HESBCN(examples_csd$csd$c4c2$data)
    ## This used to break with the strange AND from Root
    ddbroot <- do_HESBCN(examples_csd$csd$AND$data, seed = 8)
    ## Used to show inaccessible genotypes, like ABC
    ddex16 <- do_HESBCN(examples_csd$csd$c1$data, seed = 16)
    
    all_examples <- list(ex_hesbcn_and, ex_hesbcn_linear, ex_hesbcn_or, ex_hesbcn_xor,
        ex_hesbcn_c1, ex_hesbcn_c3, ex_hesbcn_c4c2, ddbroot)

    

    for (i in all_examples) {
        run_test_for_dataset(i)
    }


    ## Create two data sets that tend to lead to combination of AND, OR, XOR,
    ## Single.

    ## And if you thought CBN was unstable and sometimes counterintuitive ... try
    ## HESBCN :-)
    
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

    ## Examples that mix output
    d3_1 <- do_HESBCN(d3, seed = 26) ## AND, OR, XOR, Single
    set.seed(NULL)
    d3_3 <- do_HESBCN(d3)

    all_mixed_examples <- list(d3_1, d3_3)
    
    for(ex in all_examples){
        run_test_for_dataset(ex)
    }
    set.seed(NULL)
    
    ## Doesn't give us much, since most are Single
    ## data(every_which_way_data)
    ## see ../../other_miscell/test.run_hesbcn_continuously.R if you want that.
    
})

cat("\n Done test.HESBCN-trans-mat-against-oncosimul.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
