## This test was very kludgy as it used as reference
## a set of defaults that one might one to change.
## And that object (.ev_SHINY_dflt) is loaded implicitly.

## So now we create a function that is only used to
## generate that old default, for the pieces that are needed.
## And yes, this is OK. All these tests are transformations, that changing
## something works. So we just use the same standard that is expected.
## What changes,now, is that we are explicit about the standard (whatever
## is produced by generate_old) instead of using, in a non-explicit
## way, whatever is in a different object, which we might want to
## modify anyway.

t1 <- Sys.time()
generate_old <- function() {
    ## default_genes <- 3
    max_genes <- 10
    all_gene_names <- LETTERS[1: max_genes]
    template_dag <- matrix(0, ncol = max_genes + 1, nrow = max_genes + 1)
    rownames(template_dag) <- colnames(template_dag) <- c("Root", all_gene_names)
    template_parent_set <- rep("Single", max_genes)
    names(template_parent_set) <- all_gene_names
    template_lambdas <- rep(1, max_genes)
    names(template_lambdas) <- all_gene_names
    template_thetas <- matrix(0, ncol = max_genes, nrow = max_genes)
    rownames(template_thetas) <- colnames(template_thetas) <- all_gene_names
    template_csd_counts <- data.frame(Genotype = character(), Counts = integer())
    template_csd_data <- matrix(0, ncol=3, nrow=0)


    old_.ev_SHINY_dflt <- list(
        ## max_genes = 10,
        ## min_genes = 2,
        ## cpm_samples = 10000,
        ## ngenes = 3,
        ## csd_samples = 1000,
        ## dag_model = "HESBCN",
        ## all_cpms = c("OT", "CBN", "OncoBN", "MHN", "MCCBN", "HESBCN"),
        ## cpms2run = c("OT", "CBN", "OncoBN", "MHN"), ## , "HESBCN"),
        template_data = list(
            csd_counts =  template_csd_counts
          , data = NULL
          , dag = template_dag
          , dag_parent_set = template_parent_set
          , lambdas = template_lambdas
          , thetas = template_thetas
          , gene_names = LETTERS[1: max_genes]
          , name = "New_CSD"
        )
    )
    return(old_.ev_SHINY_dflt)
}


test_that("Modify dags works correctly",{
    ## dag <- .ev_SHINY_dflt$template_data$dag
    ## dag_parent_set <- .ev_SHINY_dflt$template_data$dag_parent_set
    old <- generate_old()
    dag <- old$template_data$dag
    dag_parent_set <- old$template_data$dag_parent_set

    
  mod1 <- dag
  mod1["Root", "A"] <- 1
  expect_equal(modify_dag(dag, "Root", "A", "add", dag_parent_set,
                          default_dag = old$template_data$dag)$dag,
               mod1)
  expect_equal(modify_dag(mod1, "Root", "A", "remove", dag_parent_set, default_dag = old$template_data$dag)$dag, dag)
  expect_error(modify_dag(NULL, "Root", "A", "remove", dag_parent_set, default_dag = old$template_data$dag)$dag, "From and To options and DAG have to be defined")
  expect_error(modify_dag(dag, NULL, "A", "remove", dag_parent_set, default_dag = old$template_data$dag)$dag, "From and To options and DAG have to be defined")
  expect_error(modify_dag(dag, "Root", NULL, "remove", dag_parent_set, default_dag = old$template_data$dag)$dag,"From and To options and DAG have to be defined")
  expect_error(modify_dag(dag, "Root2", "A", "add", dag_parent_set, default_dag = old$template_data$dag)$dag, "Both From and To options have to be valid gene names")
  expect_error(modify_dag(dag, "Root", "A2", "add", dag_parent_set, default_dag = old$template_data$dag)$dag, "Both From and To options have to be valid gene names")

  mod2 <- mod1
  mod2["A", "B"] <- 1
  expect_error(modify_dag(mod2, "B", "A", "add", dag_parent_set, default_dag = old$template_data$dag)$dag, "Relationships cannot be bidirectional")
  expect_error(modify_dag(mod2, "B", "B", "add", dag_parent_set, default_dag = old$template_data$dag)$dag, "Both From and To options must be different")
  expect_error(modify_dag(mod2, "A", "B", "add", dag_parent_set, default_dag = old$template_data$dag)$dag,"That edge is already present")
  expect_error(modify_dag(mod2, "A", "B", "add", dag_parent_set, default_dag = old$template_data$dag)$dag,"That edge is already present")
  expect_error(modify_dag(mod2, "D", "A", "add", dag_parent_set, default_dag = old$template_data$dag)$dag, "A direct children of Root cannot have multiple parents")
  expect_equal(modify_dag(mod2, "Root", "A", "clear", dag_parent_set, default_dag = old$template_data$dag)$dag,
               ## .ev_SHINY_dflt$template_data$dag
               old$template_data$dag
               )
  expect_equal(modify_dag(mod2, "Root", "A", "remove", dag_parent_set, default_dag = old$template_data$dag)$dag,
               ##.ev_SHINY_dflt$template_data$dag
               old$template_data$dag
               )

  ## More complex scenarios where we also test parent_set
  mod3 <- dag
  mod3["Root", "A"] <- 1
  mod3["A", "B"] <- 1
  mod3["A", "C"] <- 1
  mod3["B", "D"] <- 1
  mod3["C", "D"] <- 1
  tmp_parent_set3 <- dag_parent_set
  tmp_parent_set3["D"] <- "AND"
  tmp_res1 <- modify_dag(mod3, "C", "D", "remove", tmp_parent_set3, default_dag = old$template_data$dag)
  expect_equal(tmp_res1$parent_set, dag_parent_set)
  res_dag1 <- mod3
  res_dag1["C", "D"] <- 0
  expect_equal(tmp_res1$dag, res_dag1)

  mod4 <- dag
  mod4["Root", "A"] <- 1
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

  expect_error(modify_dag(mod5, "G", "B", "add", tmp_parent_set6, default_dag = old$template_data$dag), "This relationship breaks the DAG. Revise it.")

  mod5_results <- modify_dag(mod5, "D", "E", "remove", tmp_parent_set5, default_dag = old$template_data$dag)
  expect_equal(mod5_results$dag, mod4)
  expect_equal(mod5_results$parent_set, tmp_parent_set4)

  mod5_results_2 <- modify_dag(mod5, "D", "E", "remove", tmp_parent_set5, default_dag = old$template_data$dag)
  expect_equal(mod5_results_2$dag, mod4)
  expect_equal(mod5_results_2$parent_set, tmp_parent_set4)

  mod6 <- mod4
  mod6["Root", "G"] <- 1
  mod6["G", "B"] <- 1
  mod6_results <- modify_dag(mod4, "G", "B", "add", tmp_parent_set6, default_dag = old$template_data$dag)
  expect_equal(mod6_results$dag, mod6)
})

test_that("Test that modify lambdas and parent set is correct", {
  ## dag <- .ev_SHINY_dflt$template_data$dag
  ## dag_parent_set <- .ev_SHINY_dflt$template_data$dag_parent_set
  ## lambdas <- .ev_SHINY_dflt$template_data$lambdas

  old <- generate_old()
  dag <- old$template_data$dag
  dag_parent_set <- old$template_data$dag_parent_set
  lambdas <- old$template_data$lambdas

  
  dag_table1 <- data.frame(
    From = c("Root", "A", "A", "B", "C"),
    To = c("A", "B", "C", "D", "D"),
    Relation = c("Single", "Single", "Single", "AND", "AND"),
    Lambdas = c(1,2,3,4,4)
  )


  mod1 <- dag
  mod1["Root", "A"] <- 1
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

  results1 <- modify_lambdas_and_parent_set_from_table(
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

  
  expect_error(modify_lambdas_and_parent_set_from_table(
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
  
  expect_error(modify_lambdas_and_parent_set_from_table(
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

  results2 <- modify_lambdas_and_parent_set_from_table(
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
    old <- generate_old()
    dag <- old$template_data$dag
    dag_parent_set <- old$template_data$dag_parent_set


  ##   dag <- .ev_SHINY_dflt$template_data$dag
  ## dag_parent_set <- .ev_SHINY_dflt$template_data$dag_parent_set

  x <- modify_dag(dag, "Root", "A", "add", dag_parent_set, default_dag = old$template_data$dag)
  x <- modify_dag(x$dag, "Root", "B", "add", x$parent_set, default_dag = old$template_data$dag)
  x <- modify_dag(x$dag, "A", "C", "add", x$parent_set, default_dag = old$template_data$dag)
  x <- modify_dag(x$dag, "B", "C", "add", x$parent_set)
  x <- modify_dag(x$dag, "B", "D", "add", x$parent_set)
  x <- modify_dag(x$dag, "A", "D", "add", x$parent_set)

    lambdas <- old$template_data$lambdas
    ## lambdas <- .ev_SHINY_dflt$template_data$lambdas
    ## dag_parent_set <- .ev_SHINY_dflt$template_data$dag_parent_set
    
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

  results <- modify_lambdas_and_parent_set_from_table(
    dag_table, new_data, lambdas, x$dag, x$parent_set, "HESBCN"
  )
  new_parent_set <- x$parent_set
  new_parent_set["D"] <- "OR"
  expect_equal(results$parent_set, new_parent_set)
})

rm(generate_old)
cat("\n Done test.dag.building.shiny. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
