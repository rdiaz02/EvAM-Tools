library(evamtools)

#Creating inputs
## CSD - good
dB_AND <- matrix(
  c(
    rep(c(1, 0, 0, 0), 200) #A
    , rep(c(1, 0, 1, 0), 100) #AC
    , rep(c(1, 1, 0, 0), 100) #AB
    , rep(c(1, 1, 1, 0), 50) #ABC
    , rep(c(1, 1, 1, 1), 10) #ABCD
    , rep(c(0, 0, 0, 0), 10) #WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_AND) <- c("A1", "B2", "C3", "D4")
write.csv(dB_AND, file = "good_csd.csv", row.names = FALSE)

csd_data <- list(data = dB_AND)
csd2save <- evamtools:::standarize_dataset(csd_data)
csd2save$name <- "CSD_custom"
csd2save$type <- "csd"
saveRDS(csd2save, file="csd_good.RDS")

## CSD - corrupt
db_AND_corrupt <- dB_AND
db_AND_corrupt[42, 2] <- 42
write.csv(db_AND_corrupt, file = "corrupt_csd.csv", row.names = FALSE)
bad_csd2save <- csd2save
bad_csd2save$data[42, 2] <- 42
saveRDS(bad_csd2save, file="csd_bad.RDS")


## DAG - good
test_dag <- matrix(c(
    c(0, 1, 0, 0, 0),
    c(0, 0, 1, 1, 0),
    c(0, 0, 0, 0, 1),
    c(0, 0, 0, 0, 1),
    c(0, 0, 0, 0, 0)
    ), byrow = TRUE, ncol=5
)
colnames(test_dag) <- rownames(test_dag) <- c("WT", "A1", "B2", "C3", "D4")

parent_set <- c("Single", "OR", "AND", "AND")
names(parent_set) <- c("x1", "x2", "x3", "x4")
dag_data <- list(dag = test_dag
  , lambdas = c(1, 2, 3)
  , dag_parent_set = parent_set
)

dag2save <- evamtools:::standarize_dataset(dag_data)
dag2save$name <- "DAG_custom"
dag2save$type <- "dag"
saveRDS(dag2save, file="dag_good.RDS")

## DAG - corrupt
bad_dag2save <- dag2save
bad_dag2save$dag["A1", "B2"] <- 2
saveRDS(bad_dag2save, file="dag_bad.RDS")

## MHN - good
test_thetas <- matrix(c(
  c(-1, 1, 1, 3),
  c(1, 3, 2, 1),
  c(0.5, 0, 0, 1),
  c(0.7, 0, 0, 0)
  ), byrow = TRUE, ncol=4
)
colnames(test_thetas) <- rownames(test_thetas) <- c("A1", "B2", "C3", "D4")

matrix_data <- list(thetas = test_thetas)
matrix2save <- evamtools:::standarize_dataset(matrix_data)
matrix2save$name <- "MHN_custom"
matrix2save$type <- "matrix"
saveRDS(matrix2save, file="matrix_good.RDS")
## MHN - corrupt
bad_matrix2save <-  matrix2save
bad_matrix2save$thetas["A1", "B2"] <- "abc"
saveRDS(bad_matrix2save, file="matrix_bad.RDS")

## EvAM - output
methods <- c("OT", "CBN", "MCCBN", "HESBCN", "OncoBN", "MHN")
and_cpm <- evam(examples_csd$csd$AND$data, methods = methods)
and_cpm_with_simulations <- sample_CPMs(and_cpm, 10000, methods, c("sampled_genotype_counts", "obs_genotype_transitions"))
# browser()
orig_data <- list(data = examples_csd$csd$AND$data, name = "AND_new"
            , type = "csd", gene_names = colnames(examples_csd$csd$AND$data)
            , thetas = NULL, lambdas = NULL
            , dag = NULL, dag_parent_set = NULL)
tabular_data <- evamtools:::create_tabular_data(c(and_cpm, and_cpm_with_simulations))
sample_evam_output <- list("cpm_output" = c(and_cpm, and_cpm_with_simulations)
      , "orig_data" = orig_data, "tabular_data" = tabular_data
    )
# and_cpm_with_simulations$source_data <- orig_data
# and_cpm_with_simulations$name <- "AND_new"
save(sample_evam_output, file = "../../../../data/sample_evam_output.RData")


# methods <- c("OT", "CBN", "MCCBN", "HESBCN", "OncoBN", "MHN")
# dB_AND <- matrix(
#   c(
#       rep(c(1, 0, 0), 200) #A
#     ,  rep(c(0, 1, 0), 200) #B
#     , rep(c(1, 1, 0), 100) #AB
#     , rep(c(0, 1, 1), 100) #BC
#     , rep(c(1, 1, 1), 50) #ABC
#     , rep(c(0, 0, 0), 10) #WT
#   ), ncol = 3, byrow = TRUE
# )
# colnames(dB_AND) <- LETTERS[1:3]

# and_cpm <- evam(dB_AND, methods = methods)
# and_cpm_with_simulations <- sample_CPMs(and_cpm, 10000, methods, c("sampled_genotype_counts", "obs_genotype_transitions"))
# orig_data <- list(data = examples_csd$csd$AND$data, name = "AND_new"
#             , type = "csd", gene_names = colnames(examples_csd$csd$AND$data)
#             , thetas = NULL, lambdas = NULL
#             , dag = NULL, dag_parent_set = NULL)
# cpm_output <- c(and_cpm, and_cpm_with_simulations)
# sample_evam_output <- list("cpm_output" = c(and_cpm, and_cpm_with_simulations)
#       , "orig_data" = orig_data
#     )
# and_cpm_with_simulations$source_data <- orig_data
# and_cpm_with_simulations$name <- "AND_new"
# save(sample_evam_output, file = "AND_test_cpm.RData")
# save(sample_evam_output, file = "../../../../data/AND_test_cpm.RData")
# saveRDS(sample_evam_output, file = "AND_test_cpm.RDS")
