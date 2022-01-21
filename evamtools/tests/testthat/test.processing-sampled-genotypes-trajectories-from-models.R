## FIXME commented until rds issues solved
## use RData, not rds.
## test_that("Output are not generated with bad input", {
##     simGenotypes <- readRDS("../../data/simulation_output.rds")
##     x <- simGenotypes
##     x$trajectory <- NULL
##     expect_error(process_samples(x, 5), 
##         "trajectory is missing from your samples")
##     x <- simGenotypes
##     x$obs_events <- NULL
##     expect_error(process_samples(x, 5), 
##         "obs_events is missing from your samples")
## })

## test_that("Output is returned only with the requested fields", {
##     simGenotypes <- readRDS("../../data/simulation_output.rds")
##     out_params <- c("frequencies", "state_counts", "transitions")
##     out_sim <- process_samples(simGenotypes, 5, output = out_params[1])
##     expect_equal(sort(names(out_sim)), sort(out_params[1]))
    
##     out_sim <- process_samples(simGenotypes, 5, output = out_params[2:3])
##     expect_equal(sort(names(out_sim)), sort(out_params[2:3]))

##     out_sim <- process_samples(simGenotypes, 5, output = out_params)
##     expect_equal(sort(names(out_sim)), sort(out_params))

##     expect_error(process_samples(simGenotypes, 5, output = c()), "Specify valid output")

##     expect_error(process_samples(simGenotypes, 5, output = c("bad request", "bad request 2")), "Specify valid output")

##     expect_warning(process_samples(simGenotypes, 5, output = c(out_params, "bad request")), "The following parameters cannot be returned: bad request")

##     expect_warning(process_samples(simGenotypes, 5, output = c(out_params, "bad request", "bad request 2")), "The following parameters cannot be returned: bad request, bad request 2")
## })

## test_that("Output is correct", {
##     simGenotypes <- readRDS("../../data/simulation_output.rds")
##     out <- process_samples(simGenotypes, 5)
##     expect_equal(nrow(out$frequencies), 2**5)
##     expect_equal(dim(out$transitions), c(2**5, 2**5))
##     expect_equal(nrow(out$state_counts), 2**5)

##     sim_counts <- table(simGenotypes$obs_events)
##     for (i in 1:length(out$frequencies$Genotype)){
##         r <- out$frequencies[i, ]
##         if(r$Genotype %in% names(sim_counts)){
##             expect_equal(r$Counts, as.numeric(sim_counts[r$Genotype]))
##         }
##         else expect_equal(r$Counts, 0)

##     }

##     expect_equal(sum(out$frequencies$Counts), length(simGenotypes$obs_events))

##     state_counts <- table(unlist(simGenotypes$trajectory))
##     for (i in 1:length(out$state_counts$Genotype)){
##         r <- out$frequencies[i, ]
##         if(r$Genotype %in% names(sim_counts)){
##             expect_equal(r$Counts, as.numeric(sim_counts[r$Genotype]))
##         }
##         else expect_equal(r$Counts, 0)

##     }

##     expect_equal(sum(out$transitions)
##         , sum(vapply(simGenotypes$trajectory, length, numeric(1)) - 1))

##     expect_equal(out$transitions["WT", "A"], 2)
##     expect_equal(out$transitions["WT", "C"], 8)
##     expect_equal(out$transitions["C", "A, C"], 2)
##     expect_equal(out$transitions["A", "A, B"], 2)
##     expect_equal(out$transitions["C", "C, D"], 2)
##     expect_equal(out$transitions["C", "C, E"], 2)
##     expect_equal(out$transitions["C, E", "C, D, E"], 2)
##     expect_equal(out$transitions["A, B", "A, B, C"], 1)
##     expect_equal(out$transitions["A, B, C", "A, B, C, E"], 3)
##     expect_equal(out$transitions["C, D, E", "A, C, D, E"], 2)
## })
