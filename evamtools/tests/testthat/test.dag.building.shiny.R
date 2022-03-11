test_that("Modify dags works correctly",{
  dag <- SHINY_DEFAULTS$template_data$dag
  dag_parent_set <- SHINY_DEFAULTS$template_data$dag_parent_set

  mod1 <- dag
  mod1["WT", "A"] <- 1
  expect_equal(evamtools:::modify_dag(dag, "Root", "A", "add", dag_parent_set)$dag, mod1)
  expect_equal(evamtools:::modify_dag(mod1, "Root", "A", "remove", dag_parent_set)$dag, dag)
  expect_error(evamtools:::modify_dag(NULL, "Root", "A", "remove", dag_parent_set)$dag, "From and To options and DAG have to be defined")
  expect_error(evamtools:::modify_dag(dag, NULL, "A", "remove", dag_parent_set)$dag, "From and To options and DAG have to be defined")
  expect_error(evamtools:::modify_dag(dag, "WT", NULL, "remove", dag_parent_set)$dag,"From and To options and DAG have to be defined")
  expect_error(evamtools:::modify_dag(dag, "WT2", "A", "add", dag_parent_set)$dag, "Both From and To options have to be valid gene names")
  expect_error(evamtools:::modify_dag(dag, "WT", "A2", "add", dag_parent_set)$dag, "Both From and To options have to be valid gene names")

  mod2 <- mod1
  mod2["A", "B"] <- 1
  expect_error(evamtools:::modify_dag(mod2, "B", "A", "add", dag_parent_set)$dag, "Relationships cannot be bidirectional")
  expect_error(evamtools:::modify_dag(mod2, "B", "B", "add", dag_parent_set)$dag, "Both From and To options must be different")
  expect_error(evamtools:::modify_dag(mod2, "A", "B", "add", dag_parent_set)$dag,"That edge is already present")
  expect_error(evamtools:::modify_dag(mod2, "A", "B", "add", dag_parent_set)$dag,"That edge is already present")
  expect_error(evamtools:::modify_dag(mod2, "D", "A", "add", dag_parent_set)$dag, "A direct children of Root cannot have multiple parents")
  expect_equal(evamtools:::modify_dag(mod2, "WT", "A", "clear", dag_parent_set)$dag, SHINY_DEFAULTS$template_data$dag)
  expect_equal(evamtools:::modify_dag(mod2, "WT", "A", "remove", dag_parent_set)$dag, SHINY_DEFAULTS$template_data$dag)

  ## More complex scenarios where we also test parent_set
  mod3 <- dag
  mod3["WT", "A"] <- 1
  mod3["A", "B"] <- 1
  mod3["A", "C"] <- 1
  mod3["B", "D"] <- 1
  mod3["C", "D"] <- 1
  tmp_parent_set3 <- dag_parent_set
  tmp_parent_set3["D"] <- "AND"
  tmp_res1 <- evamtools:::modify_dag(mod3, "C", "D", "remove", tmp_parent_set3)
  expect_equal(tmp_res1$parent_set, dag_parent_set)
  res_dag1 <- mod3
  res_dag1["C", "D"] <- 0
  expect_equal(tmp_res1$dag, res_dag1)

  mod4 <- dag
  mod4["WT", "A"] <- 1
  mod4["A", "B"] <- 1
  mod4["A", "C"] <- 1
  mod4["B", "D"] <- 1
  mod4["C", "D"] <- 1
  mod5 <-mod4
  mod5["D", "E"] <- 1
  mod5["E", "F"] <- 1
  mod5["E", "G"] <- 1
  mod5["F", "H"] <- 1
  mod5["G", "H"] <- 1
  tmp_parent_set4 <- dag_parent_set
  tmp_parent_set4["D"] <- "AND"
  tmp_parent_set5 <- tmp_parent_set4
  tmp_parent_set5["H"] <- "OR"
  tmp_parent_set6 <- tmp_parent_set5
  tmp_parent_set6["A"] <- "Bad_rel"
  tmp_parent_set6["D"] <- "NOTAND"

  expect_error(evamtools:::modify_dag(mod5, "G", "B", "add", tmp_parent_set6), "This relationship breaks the DAG. Revise it.")

  mod5_results <- evamtools:::modify_dag(mod5, "D", "E", "remove", tmp_parent_set5)
  expect_equal(mod5_results$dag, mod4)
  expect_equal(mod5_results$parent_set, tmp_parent_set4)

  mod5_results_2 <- evamtools:::modify_dag(mod5, "D", "E", "remove", tmp_parent_set5)
  expect_equal(mod5_results_2$dag, mod4)
  expect_equal(mod5_results_2$parent_set, tmp_parent_set4)

  mod6 <- mod4
  mod6["WT", "G"] <- 1
  mod6["G", "B"] <- 1
  mod6_results <- evamtools:::modify_dag(mod4, "G", "B", "add", tmp_parent_set6)
  expect_equal(mod6_results$dag, mod6)
})

test_that("Test that modify lambdas and parent set is correct", {
  dag <- SHINY_DEFAULTS$template_data$dag
  dag_parent_set <- SHINY_DEFAULTS$template_data$dag_parent_set
  lambdas <- SHINY_DEFAULTS$template_data$lambdas

  dag_table1 <- data.frame(
    From = c("Root", "A", "A", "B", "C"),
    To = c("A", "B", "C", "D", "D"),
    Relation = c("Single", "Single", "Single", "AND", "AND"),
    Lambdas = c(1,2,3,4,4)
  )


  mod1 <- dag
  mod1["WT", "A"] <- 1
  mod1["A", "B"] <- 1
  mod1["A", "C"] <- 1
  mod1["B", "D"] <- 1
  mod1["C", "D"] <- 1

  parent_set1 <- dag_parent_set
  parent_set1["D"] <- "AND"
  lambdas1 <- lambdas
  lambdas1[1:4] <- c(1:4)

  new_data1 <- data.frame(
    row = c(rep(1,4), rep(2,4), rep(3, 4), rep(4,4), rep(5,4)),
    col = c(rep(c(0,1,2,3),5)),
    value = c(
      "Root", "A", "Single", 4, 
      "A", "B", "badRel", 2,
      "A", "C", "Single", 3,
      "B", "D", "OR", 4,
      "C", "D", "AND", 4
    )
  )

  expected_parent_set1 <- parent_set1
  expected_parent_set1["D"] <- "OR"
  expected_lambdas1 <- lambdas1
  expected_lambdas1["A"] <- 4

  results1 <- evamtools:::modify_lambdas_and_parent_set_from_table(
    dag_table1, new_data1, lambdas1, mod1, parent_set1, "HESBCN"
  )

  expect_equal(results1$lambdas, expected_lambdas1)
  expect_equal(results1$parent_set, expected_parent_set1)

  new_data2 <- data.frame(
    row = c(rep(1,4), rep(2,4), rep(3, 4), rep(4,4), rep(5,4)),
    col = c(rep(c(0,1,2,3),5)),
    value = c(
      "Root", "A", "Single", "asf", 
      "A", "B", "badRel", 2,
      "A", "C", "Single", 3,
      "B", "D", "OR", 4,
      "C", "D", "AND", 4
    )
  )
  
  expect_error(evamtools:::modify_lambdas_and_parent_set_from_table(
    dag_table1, new_data2, lambdas1, mod1, parent_set1, "HESBCN"
  ), "There are missing lambdas")

  new_data3 <- data.frame(
    row = c(rep(1,4), rep(2,4), rep(3, 4), rep(4,4), rep(5,4)),
    col = c(rep(c(0,1,2,3),5)),
    value = c(
      "Root", "A", "Single", 1, 
      "A2", "B", "badRel", 2,
      "A", "C", "Single", 3,
      "B", "D", "OR", 4,
      "C", "D", "AND", 4
    )
  )
  
  expect_error(evamtools:::modify_lambdas_and_parent_set_from_table(
    dag_table1, new_data3, lambdas1, mod1, parent_set1, "HESBCN"
  ), "There are unkown genes")


  new_data2 <- data.frame(
    row = c(rep(1,4), rep(2,4), rep(3, 4), rep(4,4), rep(5,4)),
    col = c(rep(c(0,1,2,3),5)),
    value = c(
      "Root", "A", "Single", 6, 
      "A", "B", "badRel", 2,
      "A", "C", "Single", 3,
      "B", "D", "OR", 4,
      "C", "D", "XOR", 2
    )
  )

  results2 <- evamtools:::modify_lambdas_and_parent_set_from_table(
    dag_table1, new_data2, lambdas1, mod1, parent_set1, "HESBCN"
  )
  expected_lambdas2 <- lambdas
  expected_lambdas2["A"] <- 6
  expected_lambdas2["B"] <- 2
  expected_lambdas2["C"] <- 3
  expected_lambdas2["D"] <- 2
  expect_equal(results2$lambdas, expected_lambdas2)
  expect_equal(results2$parent_set, expected_parent_set1)
})

test_that("Modify dags works correctly on a more comples example",{
  dag <- SHINY_DEFAULTS$template_data$dag
  dag_parent_set <- SHINY_DEFAULTS$template_data$dag_parent_set

  x <- evamtools:::modify_dag(dag, "WT", "A", "add", dag_parent_set)
  x <- evamtools:::modify_dag(x$dag, "WT", "B", "add", x$parent_set)
  x <- evamtools:::modify_dag(x$dag, "A", "C", "add", x$parent_set)
  x <- evamtools:::modify_dag(x$dag, "B", "C", "add", x$parent_set)
  x <- evamtools:::modify_dag(x$dag, "B", "D", "add", x$parent_set)
  x <- evamtools:::modify_dag(x$dag, "A", "D", "add", x$parent_set)
  
  lambdas <- SHINY_DEFAULTS$template_data$lambdas
  dag_parent_set <- SHINY_DEFAULTS$template_data$dag_parent_set
  dag_parent_set["C"] <- "AND"
  dag_parent_set["D"] <- "AND"
  dag_table <- data.frame(
    From = c("Root", "Root", "A", "B", "A", "B"),
    To = c("A", "B", "C", "C", "D", "D"),
    Relation = c("Single", "Single", "AND", "AND", "AND", "AND"),
    Lambdas = c(1,1,1,1,1,1)
  )

  new_data <- data.frame(
    row = c(rep(1,4), rep(2,4), rep(3, 4), rep(4,4), rep(5,4), rep(6,4)),
    col = c(rep(c(0,1,2,3),6)),
    value = c(
      "Root", "A", "Single", 1, 
      "Root", "B", "Single", 1, 
      "A", "C", "AND", 1,
      "B", "C", "AND", 1,
      "A", "D", "OR", 1,
      "B", "D", "AND", 1
    )
  )

  results <- evamtools:::modify_lambdas_and_parent_set_from_table(
    dag_table, new_data, lambdas, x$dag, x$parent_set, "HESBCN"
  )
  new_parent_set <- x$parent_set
  new_parent_set["D"] <- "OR"
  expect_equal(results$parent_set, new_parent_set)
})
