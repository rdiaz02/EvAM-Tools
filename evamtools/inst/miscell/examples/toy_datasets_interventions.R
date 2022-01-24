#' Toy dataset to play with interventions
#' Matrix with 8 gentoypes derived from a 3 genes DAG
#' Use it as template to modify the probability, and fluxes,
#' of transitioning between genotypes
template_path <- matrix(0, ncol = 8, nrow = 8)
colnames(template_path) <- rownames(template_path) <- 
  c("WT", "A", "B", "C", "A, B", "A, C", "B, C", "A, B, C")

low_transited_probability <- 0.1
template_path["WT", "WT"] <- low_transited_probability
template_path["WT", "A"] <- low_transited_probability
template_path["WT", "B"] <- low_transited_probability
template_path["WT", "C"] <- low_transited_probability
template_path["A", "A"] <- low_transited_probability
template_path["A", "A, B"] <- low_transited_probability
template_path["A", "A, C"] <- low_transited_probability
template_path["B", "B"] <- low_transited_probability
template_path["B", "A, B"] <- low_transited_probability
template_path["B", "B, C"] <- low_transited_probability
template_path["C", "C"] <- low_transited_probability
template_path["C", "A, C"] <- low_transited_probability
template_path["C", "B, C"] <- low_transited_probability
template_path["A, B", "A, B"] <- low_transited_probability
template_path["A, B", "A, B, C"] <- low_transited_probability
template_path["A, C", "A, C"] <- low_transited_probability
template_path["A, C", "A, B, C"] <- low_transited_probability
template_path["B, C", "B, C"] <- low_transited_probability
template_path["B, C", "A, B, C"] <- low_transited_probability
template_path["A, B, C", "A, B, C"] <- 1

# Going from WT -> A -> AB -> ABC
path_simple_tpm <- template_path
path_simple_tpm["WT", "A"] <- 0.7 
path_simple_tpm["A", "A, B"] <- 0.8 
path_simple_tpm["A, B", "A, B, C"] <- 0.9 
path_simple_tpm["B", "B"] <- 0.8 
path_simple_tpm["C", "C"] <- 0.8 
path_simple_tpm["A, C", "A, C"] <- 0.9 
path_simple_tpm["B, C", "B, C"] <- 0.9 

# Going from WT -> A -> AC -> ABC
# Going from WT -> C -> AC -> ABC
path_converging_tpm <- template_path
path_converging_tpm["WT", "A"] <- 0.4 
path_converging_tpm["WT", "C"] <- 0.4 
path_converging_tpm["A", "A, C"] <- 0.8 
path_converging_tpm["C", "A, C"] <- 0.8 
path_converging_tpm["A, C", "A, B, C"] <- 0.9 
path_converging_tpm["B", "B"] <- 0.8 
path_converging_tpm["A, B", "A, B"] <- 0.9 
path_converging_tpm["B, C", "B, C"] <- 0.9 


# Going from WT -> B -> AB -> ABC
# Going from WT -> B -> BC -> ABC
path_diverging_tpm <- template_path
path_diverging_tpm["WT", "B"] <- 0.7
path_diverging_tpm["A", "A"] <- 0.8 
path_diverging_tpm["C", "C"] <- 0.8 
path_diverging_tpm["B", "A, B"] <- 0.45 
path_diverging_tpm["B", "B, C"] <- 0.45 
path_diverging_tpm["A, C", "A, C"] <- 0.9 
path_diverging_tpm["A, B", "A, B, C"] <- 0.9 
path_diverging_tpm["B, C", "A, B, C"] <- 0.9 


# Going from WT -> A -> AB -> ABC
# Going from WT -> C -> BC -> ABC
path_parallel_tpm <- template_path
path_parallel_tpm["WT", "A"] <- 0.3
path_parallel_tpm["WT", "B"] <- 0.3
path_parallel_tpm["WT", "C"] <- 0.3
path_parallel_tpm["B", "B"] <- 0.8 
path_parallel_tpm["A", "A, B"] <- 0.8
path_parallel_tpm["C", "B, C"] <- 0.8 
path_parallel_tpm["A, C", "A, C"] <- 0.9 
path_parallel_tpm["A, B", "A, B, C"] <- 0.9 
path_parallel_tpm["B, C", "A, B, C"] <- 0.9 


transition_probability_matrix2count_matrix <- function(tpm, start_size = 10000){
  if(!all(rowSums(tpm) == 1)) stop("Not all rows sum one")
  count_matrix <- tpm * 0
  count_matrix[1, ] <- tpm[1, ] * start_size
  for (i in 2:ncol(tpm)){
    count_matrix[i, ] <- sum(count_matrix[, i]) * tpm[i, ]
  }

  return(count_matrix)
}

# Going from WT -> A -> AB -> ABC
path_simple <- transition_probability_matrix2count_matrix(path_simple_tpm)
# Going from WT -> B -> AB -> ABC
# Going from WT -> B -> BC -> ABC
path_diverging <- transition_probability_matrix2count_matrix(path_diverging_tpm)
# Going from WT -> A -> AC -> ABC
# Going from WT -> C -> AC -> ABC
path_converging <- transition_probability_matrix2count_matrix(path_converging_tpm)
# Going from WT -> A -> AB -> ABC
# Going from WT -> C -> BC -> ABC
path_parallel <- transition_probability_matrix2count_matrix(path_parallel_tpm)

#' List with everything
all_examples_interventions <- list(path_simple = path_simple
  , path_converging = path_converging
  , path_diverging = path_diverging
  , path_parallel = path_parallel
  )

# save(path_simple, file = "toy_datasets_interventions/path_simple.RData")
# save(path_converging, file = "toy_datasets_interventions/path_converging.RData")
# save(path_diverging, file = "toy_datasets_interventions/path_diverging.RData")
# save(path_parallel, file = "toy_datasets_interventions/path_parallel.RData")