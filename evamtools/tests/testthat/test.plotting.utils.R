t1 <- Sys.time()

test_that("Processing samples & Plotting of CPMs is correct", {
  sample_evam_output <- evam(examples_csd$csd$AND$data)
  models <- c("OT", "CBN", "OncoBN", "HESBCN", "MHN")

  # Does not process
  expect_error(process_data(NULL, "OT", "something_not_supported"))
  expect_error(process_data(sample_evam_output, "OT", "something_not_supported"))

  # Checking without sampling
  for(model in models){
      model_transitions <- process_data(sample_evam_output, model,
                                                    "obs_genotype_transitions")
      model_trm <- process_data(sample_evam_output, model,
                                            "trans_rate_mat")
      model_probabilities <- process_data(sample_evam_output, model,
                                                      "trans_mat")

    expect_equal(model_transitions$data2plot, NULL) ## We have not sampled
    expect_equal(model_probabilities$data2plot, sample_evam_output[[paste0(model, "_trans_mat")]])

    if(model %in% c("OT", "OncoBN")){
      expect_equal(model_trm$data2plot, NULL)
    } else {
      expect_equal(model_trm$data2plot, sample_evam_output[[paste0(model, "_trans_rate_mat")]]) 
    }

    for(i in c("parent_set")){
      expect_equal(model_trm[[i]], sample_evam_output[[paste0(model, "_", i)]]) 
      expect_equal(model_transitions[[i]], sample_evam_output[[paste0(model, "_", i)]]) 
      expect_equal(model_probabilities[[i]], sample_evam_output[[paste0(model, "_", i)]]) 
    }
  }

  # Checking with sampling
  samples <- sample_evam(sample_evam_output, 100)

  for(model in models) {
      cat("\n model ", model, "\n")
      model_transitions <- process_data(sample_evam_output, model,
                                                    "obs_genotype_transitions", samples)
      model_trm <- process_data(sample_evam_output, model,
                                            "trans_rate_mat", samples)
      model_probabilities <- process_data(sample_evam_output, model,
                                                      "trans_mat", samples)

      expect_equal(model_probabilities$data2plot,
                   sample_evam_output[[paste0(model, "_trans_mat")]])

    if(model %in% c("OT", "OncoBN")){
      expect_equal(model_transitions$data2plot, NULL) 
      expect_equal(model_trm$data2plot,
                   sample_evam_output[[paste0(model, "_trans_rate_mat")]])
    } else {
      tmp <- model_transitions$data2plot
      expect_equal(tmp, samples[[paste0(model, "_obs_genotype_transitions")]])
      expect_equal(model_trm$data2plot,
                   sample_evam_output[[paste0(model, "_trans_rate_mat")]]) 
    }

    for(i in c("parent_set", "predicted_genotype_freqs")){
      expect_equal(model_trm[[i]], sample_evam_output[[paste0(model, "_", i)]]) 
      expect_equal(model_transitions[[i]], sample_evam_output[[paste0(model, "_", i)]]) 
      expect_equal(model_probabilities[[i]], sample_evam_output[[paste0(model, "_", i)]]) 
    }

    for(i in c("sampled_genotype_counts")){
      expect_equal(model_trm[[i]], samples[[paste0(model, "_", i)]]) 
      expect_equal(model_transitions[[i]], samples[[paste0(model, "_", i)]]) 
      expect_equal(model_probabilities[[i]], samples[[paste0(model, "_", i)]]) 
    }
  }
})

test_that("how it handles processing missing data", {
  methods <- c("OT", "CBN", "HESBCN", "MHN")
  
  sample_evam_output <- evam(examples_csd$csd$AND$data, methods = methods)
  out <- process_data(sample_evam_output, "OncoBN", "trans_mat", NULL)
  expect_equal(out$method_info, NULL)
  expect_equal(out$data2plot, NA)
  expect_equal(out$predicted_genotype_freqs, NA)
  expect_equal(out$parent_set, NA)
  expect_equal(out$sampled_genotype_counts, NULL)
  expect_equal(out$edges, NA)
  
  out <- process_data(sample_evam_output, "MCCBN", "trans_mat", NULL)
  expect_equal(out$method_info, NULL)
  expect_equal(out$data2plot, NA)
  expect_equal(out$predicted_genotype_freqs, NA)
  expect_equal(out$parent_set, NULL)
  expect_equal(out$sampled_genotype_counts, NULL)
  expect_equal(out$edges, NA)
})

test_that("plot_evam handles argument correctly", {
  sample_evam_output <- evam(examples_csd$csd$AND$data)
  expect_error(plot_evam(sample_evam_output, orientation="horizontal",
                         plot_type = "not_supported"))
  expect_error(plot_evam(sample_evam_output, orientation="horizontal",
                         plot_type = "transitions"))
})


test_that("plotting works, minimal, with mixed edges", {
    ## Yes, in the help file but make sure we run this
    data("ex_mixed_and_or_xor")

    out_AND_OR_XOR <- evam(ex_mixed_and_or_xor,
                       methods = c("OT", "HESBCN", "MHN", "OncoBN"),
                       hesbcn_opts = list(seed = 26))

    plot_evam(out_AND_OR_XOR, plot_type = "trans_mat", top_paths = 4)

    ## Now, only some, and change options
    plot_evam(out_AND_OR_XOR, methods = "OT", top_paths = 2,
              orientation = "vertical")
    expect_warning(plot_evam(out_AND_OR_XOR, methods = c("OT", "CBN"),
                             top_paths = 2,
                             orientation = "vertical"),
                   "At least one method you asked",
                   fixed = TRUE)
    
    expect_error(plot_evam(out_AND_OR_XOR, methods = c("CBN"),
                           top_paths = 2,
                           orientation = "vertical"),
                 "No valid methods", fixed = TRUE)
    
    expect_error(plot_evam(NULL),
                 "No valid methods", fixed = TRUE)
    
})

test_that("Exercise other plotting options", {
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
    out <- evam(dB_c1[, 1:3],
                methods = c("CBN", "OT", "MHN"))

    plot_evam(out, plot_type = "trans_mat")
    plot_evam(out, plot_type = "trans_rate_mat")
    out_samp <- sample_evam(out, 100,
                            output = c("sampled_genotype_counts",
                                       "obs_genotype_transitions"))
    plot_evam(out, out_samp, plot_type = "obs_genotype_transitions")
    plot_evam(out, out_samp, plot_type = "obs_genotype_transitions", 
              label_type = "acquisition")

    expect_error(plot_evam(out, samples = NULL,
                           plot_type = "obs_genotype_transitions"),
                 "obs_genotype_transitions needs", fixed = TRUE)
})


test_that("Check implementation of graphAM for DAG", {
     dfe <- data.frame(From = c("Root", "A"),
                                   To = c("A", "B"),
                                   lambda = c(3, 4),
                                   OT_edgeWeight = c(.2, .7))
     expect_error(DAG_plot_graphAM(dfe, NULL),
                  "more than one column with weights",
                  fixed = TRUE)

     dfe2 <- data.frame(From = c("Root", "A"),
                                     To = c("A", "B"),
                                     lambda = c(0, 0.00001))
     expect_silent(DAG_plot_graphAM(dfe2, "something"))

     dfe3 <- data.frame(From = c("Root", "A"),
                        To = c("A", "B"))
     expect_silent(DAG_plot_graphAM(dfe3, "something"))
})

test_that("plot_genotype_counts", {
    ## This cannot be tested with the usual code
    ## as only called from shiny

    dx1 <- data.frame(Genotype = c("A", "A, B, C, D", "E"),
                      Freq = c(0.1, 0.2, 0.7))
    expect_silent(plot_genotype_counts(dx1))

    dx2 <- data.frame(Genotype = c("A", "A, B, C, D", "E"),
                      Counts = c(8, 100, 50))
    expect_silent(plot_genotype_counts(dx2))

    dx3 <- data.frame()
    expect_silent(plot_genotype_counts(dx3))
})




cat("\n Done test.plotting.utils.R.Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
