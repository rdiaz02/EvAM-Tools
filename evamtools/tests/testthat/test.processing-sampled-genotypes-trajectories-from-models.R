t1 <- Sys.time()

test_that("Output is not generated with bad input", {
    trajectory <- list(
        c("WT", "C", "C, D"),
        c("WT", "C", "C, E", "C, D, E", "A, C, D, E"),
        c("WT", "A", "A, B"),
        c("WT", "C"),
        c("WT"),
        c("WT", "A", "A, B", "A, B, C", "A, B, C, E"),
        c("WT", "C", "C, E", "C, D, E", "A, C, D, E"),
        c("WT", "C", "C, D"),
        c("WT", "C", "A, C", "A, B, C", "A, B, C, E"),
        c("WT", "C"), 
        c("WT", "C", "A, C", "A, B, C", "A, B, C, E")
    )

    expected_transitions <- c(8, 2, 2, 2, 2, 2, 2, 1, 3, 2, 2)
    names(expected_transitions) <- c("WT -> C", "C -> C, D",
        "C -> C, E", "C, E -> C, D, E", "C, D, E -> A, C, D, E", 
        "WT -> A", "A -> A, B", "A, B -> A, B, C", "A, B, C -> A, B, C, E", 
        "C -> A, C", "A, C -> A, B, C")

    obs_events <- c(
        "C, D", "A, C, D, E", "A, B", "C", "WT", "A, B, C, E",
        "A, C, D, E", "C, D", "A, B, C, E", "C", "A, B, C, E"
    )

    sampling_output <- list(
        trajectory = trajectory, 
        obs_events = obs_events
    )

    x <- sampling_output
    x$trajectory <- NULL
    expect_error(process_samples(x, 5, gene_names = LETTERS[1:5]), 
                 "trajectory is missing from your samples")
    x <- sampling_output
    x$obs_events <- NULL
    expect_error(process_samples(x, 5, gene_names = LETTERS[1:5]), 
                 "obs_events is missing from your samples")
})

test_that("Output is returned only with the requested fields", {

    trajectory <- list(
        c("WT", "C", "C, D"),
        c("WT", "C", "C, E", "C, D, E", "A, C, D, E"),
        c("WT", "A", "A, B"),
        c("WT", "C"),
        c("WT"),
        c("WT", "A", "A, B", "A, B, C", "A, B, C, E"),
        c("WT", "C", "C, E", "C, D, E", "A, C, D, E"),
        c("WT", "C", "C, D"),
        c("WT", "C", "A, C", "A, B, C", "A, B, C, E"),
        c("WT", "C"), 
        c("WT", "C", "A, C", "A, B, C", "A, B, C, E")
    )

    expected_transitions <- c(8, 2, 2, 2, 2, 2, 2, 1, 3, 2, 2)
    names(expected_transitions) <- c("WT -> C", "C -> C, D",
        "C -> C, E", "C, E -> C, D, E", "C, D, E -> A, C, D, E", 
        "WT -> A", "A -> A, B", "A, B -> A, B, C", "A, B, C -> A, B, C, E", 
        "C -> A, C", "A, C -> A, B, C")

    obs_events <- c(
        "C, D", "A, C, D, E", "A, B", "C", "A, B, C, E",
        "A, C, D, E", "C, D", "A, B, C, E", "C", "A, B, C, E"
    )


    simGenotypes <- list(
        trajectory = trajectory,
        obs_events = obs_events
    )

    out_params <- c("sampled_genotype_counts",
                    "state_counts",
                    "obs_genotype_transitions")
    ## out_params <- c("frequencies", "state_counts", "transitions")
    out_sim <- process_samples(simGenotypes, 5,
                                           gene_names = LETTERS[1:5],
                                           output = out_params[1])
    expect_equal(sort(names(out_sim)), sort(out_params[1]))
    
    out_sim <- process_samples(simGenotypes, 5,
                                           gene_names = LETTERS[1:5],
                                           output = out_params[2:3])
    expect_equal(sort(names(out_sim)), sort(out_params[2:3]))

    out_sim <- process_samples(simGenotypes, 5,
                                           gene_names = LETTERS[1:5],
                                           output = out_params)
    expect_equal(sort(names(out_sim)), sort(out_params))


    expect_error(process_samples(simGenotypes, 5,
                                             gene_names = LETTERS[1:5],
                                             output = c()),
                 "No output specified", fixed = TRUE)

    expect_error(process_samples(simGenotypes, 5,
                                             gene_names = LETTERS[1:5],
                                             output = c("bad request", "bad request 2")),
                 "Incorrect output specified", fixed = TRUE)

    expect_error(process_samples(simGenotypes, 5,
                                               gene_names = LETTERS[1:5],
                                               output = c(out_params, "bad request")),
                    "Incorrect output specified", fixed = TRUE)


    expect_error(process_samples(simGenotypes, 5,
                                               gene_names = LETTERS[1:5],
                                               output = c(out_params, "bad request",
                                                          "bad request 2")),
                    "Incorrect output specified", fixed = TRUE)
    
})


test_that("Output is correct", {
     trajectory <- list(
        c("WT", "C", "C, D"),
        c("WT", "C", "C, E", "C, D, E", "A, C, D, E"),
        c("WT", "A", "A, B"),
        c("WT", "C"),
        c("WT"),
        c("WT", "A", "A, B", "A, B, C", "A, B, C, E"),
        c("WT", "C", "C, E", "C, D, E", "A, C, D, E"),
        c("WT", "C", "C, D"),
        c("WT", "C", "A, C", "A, B, C", "A, B, C, E"),
        c("WT", "C"), 
        c("WT", "C", "A, C", "A, B, C", "A, B, C, E")
    )

    expected_transitions <- c(8, 2, 2, 2, 2, 2, 2, 1, 3, 2, 2)
    names(expected_transitions) <- c("WT -> C", "C -> C, D",
        "C -> C, E", "C, E -> C, D, E", "C, D, E -> A, C, D, E", 
        "WT -> A", "A -> A, B", "A, B -> A, B, C", "A, B, C -> A, B, C, E", 
        "C -> A, C", "A, C -> A, B, C")

    obs_events <- c(
        "C, D", "A, C, D, E", "A, B", "C", "A, B, C, E",
        "A, C, D, E", "C, D", "A, B, C, E", "C", "A, B, C, E"
    )

    simGenotypes <- list(
        trajectory = trajectory,
        obs_events = obs_events
    )
    
   
    out <- process_samples(simGenotypes, 5, gene_names = LETTERS[1:5])
    expect_equal(length(out$sampled_genotype_counts), 2**5)
    expect_equal(nrow(out$state_counts), 2**5)

    sim_counts <- table(simGenotypes$obs_events)
    r <- out$sampled_genotype_counts
    for (i in names(r)){
        if(i %in% names(sim_counts)){
            expect_equal(r[i], sim_counts[i])
        }
        else{
            expect_equal(as.numeric(r[i]), 0)
        }
    }

    expect_equal(sum(out$sampled_genotype_counts), length(simGenotypes$obs_events))

    state_counts <- table(unlist(simGenotypes$trajectory))
    rownames(out$state_counts) <- out$state_counts$Genotype 
    for (i in out$state_counts$Genotype){
        tmp_row <- out$state_counts[i, ]
        if(tmp_row$Genotype %in% names(state_counts)){
            expect_equal(tmp_row$Counts, as.numeric(state_counts[tmp_row$Genotype]))
        }
        else expect_equal(tmp_row$Counts, 0)
    }

    expect_equal(sum(out$obs_genotype_transitions)
        , sum(vapply(simGenotypes$trajectory, length, numeric(1)) - 1))

     ## Note: this test has many identical entries. Not a strong one.
     ## This expects a vector
     ## expect_equal(out$transitions, expected_transitions)
     tmp <- igraph::as_data_frame(igraph::graph_from_adjacency_matrix(out$obs_genotype_transitions))
     tmp2 <- paste0(tmp[, 1], " -> ", tmp[, 2])
     ttmp2 <- table(tmp2)
     tmp2v <- as.vector(ttmp2)
     names(tmp2v) <- names(ttmp2)
     expect_equal(expected_transitions,
                  tmp2v[names(expected_transitions)])
     
})


test_that("New algorithm for state transitions", {
    sim1 <- list(trajectory = list(c("WT", "A"), c("WT", "A"), c("WT", "A")),
                 obs_events = c("A", "A", "A"))

    sim2 <- list(trajectory = list(c("WT", "A"), c("WT", "A"), c("WT", "A"), c("WT")),
                 obs_events = c("A", "A", "A", "WT"))

    sim3 <- list(trajectory = list(c("WT", "A"), c("WT", "A"), c("WT"),
                                   c("WT", "A", "B")),
                 obs_events = c("A", "A", "B", "WT"))

    sim4 <- list(trajectory = list(c("WT", "A"), c("WT", "A"), c("WT"),
                                   c("WT", "A", "B"), c("WT", "A", "B")),
                 obs_events = c("A", "A", "B", "WT", "B"))

    sim5 <- list(trajectory = list(c("WT", "A"), c("WT", "A"), c("WT"),
                                   c("WT", "A", "B"), c("WT", "A", "B", "C")),
                 obs_events = c("A", "A", "B", "WT", "C"))

    sim6 <- list(trajectory = list(c("WT", "A"), c("WT", "A"), c("WT"),
                                   c("WT", "A", "B", "C"), c("WT", "A", "B", "C")),
                 obs_events = c("A", "A", "WT", "C", "C"))

    t1 <- Matrix::sparseMatrix(i = NULL, j = NULL, dims = c(2, 2),
                               x = 0, dimnames = list(c("WT", "A"), c("WT", "A")))
    t1[1, 2] <- 3

    t2 <- Matrix::sparseMatrix(i = NULL, j = NULL, dims = c(2, 2),
                               x = 0, dimnames = list(c("WT", "A"), c("WT", "A")))
    t2[1, 2] <- 3
    
    t3 <- Matrix::sparseMatrix(i = NULL, j = NULL, dims = c(3, 3),
                               x = 0, dimnames = list(c("WT", "A", "B"),
                                                      c("WT", "A", "B")))
    t3[1, 2] <- 3
    t3[2, 3] <- 1

    t4 <- Matrix::sparseMatrix(i = NULL, j = NULL, dims = c(3, 3),
                               x = 0, dimnames = list(c("WT", "A", "B"),
                                                      c("WT", "A", "B")))
    t4[1, 2] <- 4
    t4[2, 3] <- 2

    t5 <- Matrix::sparseMatrix(i = NULL, j = NULL, dims = c(4, 4),
                               x = 0, dimnames = list(c("WT", "A", "B", "C"),
                                                      c("WT", "A", "B", "C")))
    t5[1, 2] <- 4
    t5[2, 3] <- 2
    t5[3, 4] <- 1

    t6 <- Matrix::sparseMatrix(i = NULL, j = NULL, dims = c(4, 4),
                               x = 0, dimnames = list(c("WT", "A", "B", "C"),
                                                      c("WT", "A", "B", "C")))
    t6[1, 2] <- 4
    t6[2, 3] <- 2
    t6[3, 4] <- 2
    
    os1 <- process_samples(sim1, 5, gene_names = LETTERS[1:5])
    os2 <- process_samples(sim2, 5, gene_names = LETTERS[1:5])
    os3 <- process_samples(sim3, 5, gene_names = LETTERS[1:5])
    os4 <- process_samples(sim4, 5, gene_names = LETTERS[1:5])
    os5 <- process_samples(sim5, 5, gene_names = LETTERS[1:5])
    os6 <- process_samples(sim6, 5, gene_names = LETTERS[1:5])
    
    expect_equal(os1$obs_genotype_transitions, t1)
    expect_equal(os2$obs_genotype_transitions, t2)
    expect_equal(os3$obs_genotype_transitions, t3)
    expect_equal(os4$obs_genotype_transitions, t4)
    expect_equal(os5$obs_genotype_transitions, t5)
    expect_equal(os6$obs_genotype_transitions, t6)
})

cat("\n Done test.processing-sampled-genotypes-trajectories-from-models.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
