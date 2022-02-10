

test_that("Processing samples & Plotting of CPMs is correct", {
  sample_evam_output <- evam(examples_csd$csd$AND$data)


  models <- c("OT", "CBN", "OncoBN", "HESBCN", "MHN")

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
      expect_equal(model_transitions$data2plot, samples[[sprintf("%s_genotype_transitions", model)]])
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

  #Checking orientation
  # plot_CPMs(sample_evam_output, orientation="horizontal")
  # plot_CPMs(sample_evam_output, orientation="vertical")


  #Checking models to plot
  # plot_CPMs(sample_evam_output, models = c("OT", "CBN", "OncoBN"), orientation="horizontal")
  # plot_CPMs(sample_evam_output, models = c("OT", "CBN", "OncoBN"),orientation="vertical")

  # #Checking type of plot
  # plot_CPMs(sample_evam_output, models = c("OT", "CBN", "OncoBN"), orientation="horizontal", plot_type="probabilities")
  # plot_CPMs(sample_evam_output, models = c("OT", "CBN", "OncoBN"),orientation="vertical", plot_type="trm")
  # plot_CPMs(sample_evam_output, models = c("OT", "CBN", "OncoBN"),orientation="vertical", plot_type="transitions")

  # expected_not_null_outputs <- list(
  #   OT = c("dag_tree", "")
  # )
  # processed_output <- evamtools:::process_data(sample_evam_output, "OT")
  # processed_output <- evamtools:::process_data(sample_evam_output, "OncoBN")
  # processed_output <- evamtools:::process_data(sample_evam_output, "MHN")
  # processed_output <- evamtools:::process_data(sample_evam_output, "HyperTraps")

  # processed_output <- evamtools:::process_data(sample_evam_output, "bad_name")
})


# test_that("CPM layout", {

# })