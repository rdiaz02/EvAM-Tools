pwd2 <- getwd()
setwd("../../data/")
source("toy_datasets.R")
setwd("../examples/")
source("access_genots_from_oncosimul.R")
setwd(pwd2)
rm(pwd2)

compare_DBN_cpm2tm <- function(codename){
    out <- readRDS(sprintf("../../data/toy_datasets_cpms/%s.rds", codename))
    out_onco <- cpm2tm(out$DBN_model)

    order1 <- sort(rownames(out$DBN_trans_mat), index.return = TRUE)$ix
    order2 <- sort(rownames(out_onco$transition_matrix), index.return = TRUE)$ix

    ordered_computed_trm <- out$DBN_trans_mat[order1, order1]
    ordered_trm_onco <- out_onco$transition_matrix[order2, order2]

    expect_equal(ordered_computed_trm, ordered_trm_onco)
}

test_that("DBN gives the same results as OncoSimul", {
   for (i in names(all_examples)){
        print(sprintf("Dataset %s", i))
        compare_methods <- compare_DBN_cpm2tm(i)
    }
})