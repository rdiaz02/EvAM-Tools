# test_that("Test standarize datasets", {
#   check_fields <- function(orig_data, expected_data){
#     standard_data <- standarize_dataset(orig_data)

#     expected_attr <- names(expected_data)
#     for(i in expected_attr){
#       expect_equal(standard_data$i, expected_data$i)
#     }

#     orig_attr <- names(orig_data)
#     for (i in setdiff(orig_attr, expected_attr)){
#       expect_equal(standard_data$i, orig_data$i)
#     }

#     already_checked_attr <- unique(c(orig_attr, expected_attr))
#     for(i in setdiff(names(standard_data), already_checked_attr)){
#       ## The rest has to be equal to SHINY_DEFAULTS
#       if(i == "gene_names" & "data" %in% names(orig_data)) next
#       expect_equal(standard_data$i, SHINY_DEFAULTS$i)
#     }
#     expect_equal(standard_data$gene_names, colnames(standard_data$thetas))
#     expect_equal(c("WT", standard_data$gene_names), colnames(standard_data$dag))
#     expect_equal(standard_data$gene_names, names(standard_data$dag_parent_set))
#   }

#   test_data <- examples_csd$csd$AND$data

#   data1 <- NULL

#   # Try with a CSD
#   data2 <- list(data = test_data)

#   tmp_data3 <- test_data
#   colnames(tmp_data3) <- c("A1", "C2", "D3", "X4")
#   data3 <- list(data = tmp_data3)
#   expected_gene_names_3 <- c("A1", "C2", "D3", "X4",  "E", "F", "G", "H", "I", "J")

#   tmp_lambdas <- c(1,2,3,4)
#   names(tmp_lambdas) <- c("z", "x", "y", "w")
#   data4 <- list(lambdas = tmp_lambdas, data = tmp_data3) ## Different gene names

#   expected_lambdas_4 <- rep(1, 10)
#   expected_lambdas_4[1:4] <- tmp_lambdas
#   names(expected_lambdas_4) <- c("A1", "C2", "D3", "X4", "E", "F", "G", "H", "I", "J")

#   # Try with a DAG
#   test_dag <- matrix(c(
#     c(0, 1, 1, 0),
#     c(0, 0, 0, 1),
#     c(0, 0, 0, 1),
#     c(0, 0, 0, 0)
#     ), byrow = TRUE, ncol=4
#   )
#   colnames(test_dag) <- rownames(test_dag) <- c("A1", "B2", "C3", "D4")
  
#   data5 <- list(dag = test_dag
#     , lambdas = c(1, 2, 3)
#     , parent_set = c("Single", "OR", "AND")
#   )

#   gene_names_5 <- c("A1", "B2", "C3", "D4", "E", "F", "G", "H", "I", "J")
#   parent_set_5 <- rep("Single", 10)
#   names(parent_set_5) <- gene_names_5 

#   expected_data_5 <- list(
#     dag_parent_set = parent_set_5,
#     gene_names = gene_names_5
#   )

#   ## Try with a Matrix
#   test_thetas <- matrix(c(
#     c(-1, 1, 1, 3),
#     c(1, 3, 2, 1),
#     c(0.5, 0, 0, 1),
#     c(0.7, 0, 0, 0)
#     ), byrow = TRUE, ncol=4
#   )
#   colnames(test_thetas) <- rownames(test_thetas) <- c("A1", "B2", "C3", "D4")
#   gene_names_6 <- c("A1", "B2", "C3", "D4", "E", "F", "G", "H", "I", "J")

#   expected_thetas_6 <- matrix(0, ncol=10, nrow = 10)
#   colnames(expected_thetas_6) <- rownames(expected_thetas_6) <- gene_names_6

#   expected_thetas_6[1:4, 1:4] <- test_thetas
#   expected_data_6 <- list(thetas = expected_thetas_6, gene_names = gene_names_6)
#   # test_thetas <- examples_csd$matrix$test1$thetas
#   # data5 <- list(thetas = test_thetas)

#   data2test <- list(
#     # list(orig = data1, expected = list()), 
#     # list(orig = data2, expected = list()),
#     # list(orig = data3, expected = list(gene_names = expected_gene_names_3))
#     # list(orig = data4, expected = list(
#     #   gene_names = expected_gene_names_3, 
#     #   lambdas = expected_lambdas_4))
#     list(orig = data5, expected = expected_data_5)
#     # list(orig = data6, expected = expected_data_6)
#   )

#   for(i in data2test){
#     check_fields(i$orig, i$expected)
#   }
# })