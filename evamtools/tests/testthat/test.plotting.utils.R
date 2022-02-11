

test_that("Processing samples & Plotting of CPMs is correct", {
  sample_evam_output <- evam(examples_csd$csd$AND$data)
  models <- c("OT", "CBN", "OncoBN", "HESBCN", "MHN")

  # Does not process
  expect_error(evamtools:::process_data(NULL, "OT", "something_not_supported"))
  expect_error(evamtools:::process_data(sample_evam_output, "OT", "something_not_supported"))

  # Checking without sampling
  for(model in models){
    model_transitions <- evamtools:::process_data(sample_evam_output, model, "transitions")
    model_trm <- evamtools:::process_data(sample_evam_output, model, "trm")
    model_probabilities <- evamtools:::process_data(sample_evam_output, model, "probabilities")

    expect_equal(model_transitions$data2plot, NULL) ## We have not sampled
    expect_equal(model_probabilities$data2plot, sample_evam_output[[sprintf("%s_trans_mat", model)]])

    if(model %in% c("OT", "OncoBN")){
      expect_equal(model_trm$data2plot, NULL)
    } else {
      expect_equal(model_trm$data2plot, sample_evam_output[[sprintf("%s_trans_rate_mat", model)]]) 
    }

    for(i in c("theta", "parent_set")){
      expect_equal(model_trm[[i]], sample_evam_output[[sprintf("%s_%s", model, i)]]) 
      expect_equal(model_transitions[[i]], sample_evam_output[[sprintf("%s_%s", model, i)]]) 
      expect_equal(model_probabilities[[i]], sample_evam_output[[sprintf("%s_%s", model, i)]]) 
    }
  }

  # Checking with sampling
  samples <- evamtools:::sample_all_CPMs(sample_evam_output, 100)

  for(model in models){
    model_transitions <- evamtools:::process_data(sample_evam_output, model, "transitions", samples)
    model_trm <- evamtools:::process_data(sample_evam_output, model, "trm", samples)
    model_probabilities <- evamtools:::process_data(sample_evam_output, model, "probabilities", samples)

    expect_equal(model_probabilities$data2plot, sample_evam_output[[sprintf("%s_trans_mat", model)]])

    if(model %in% c("OT", "OncoBN")){
      expect_equal(model_transitions$data2plot, NULL) 
      expect_equal(model_trm$data2plot,sample_evam_output[[sprintf("%s_trans_rate_mat", model)]])
    } else {
      tmp <- model_transitions$data2plot$Counts
      names(tmp) <- rownames(model_transitions$data2plot)
      expect_equal(tmp, samples[[sprintf("%s_genotype_transitions", model)]])
      expect_equal(model_trm$data2plot, sample_evam_output[[sprintf("%s_trans_rate_mat", model)]]) 
    }

    for(i in c("theta", "parent_set")){
      expect_equal(model_trm[[i]], sample_evam_output[[sprintf("%s_%s", model, i)]]) 
      expect_equal(model_transitions[[i]], sample_evam_output[[sprintf("%s_%s", model, i)]]) 
      expect_equal(model_probabilities[[i]], sample_evam_output[[sprintf("%s_%s", model, i)]]) 
    }

    for(i in c("genotype_freqs")){
      expect_equal(model_trm[[i]], samples[[sprintf("%s_%s", model, i)]]) 
      expect_equal(model_transitions[[i]], samples[[sprintf("%s_%s", model, i)]]) 
      expect_equal(model_probabilities[[i]], samples[[sprintf("%s_%s", model, i)]]) 
    }
  }
})

test_that("plot_CPMs handles argument correctly", {
  sample_evam_output <- evam(examples_csd$csd$AND$data)
  expect_error(plot_CPMs(sample_evam_output, orientation="horizontal", plot_type = "not_supported"))
  expect_error(plot_CPMs(sample_evam_output, orientation="horizontal", plot_type = "transitions"))
})