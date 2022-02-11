library(evamtools)

and_cpm <- evam(examples_csd$csd$AND$data)
and_cpm_with_simulations <- sample_CPMs(and_cpm, 10000, 5)
orig_data <- list(data = examples_csd$csd$AND$data, name = "AND_new"
            , type = "csd", gene_names = colnames(examples_csd$csd$AND$data)
            , thetas = NULL, lambdas = NULL
            , dag = NULL, dag_parent_set = NULL)
and_cpm_with_simulations$source_data <- orig_data
and_cpm_with_simulations$name <- "AND_new"
save(and_cpm_with_simulations, file = "AND_test_cpm.RData")
saveRDS(and_cpm_with_simulations, file = "AND_test_cpm.RDS")

