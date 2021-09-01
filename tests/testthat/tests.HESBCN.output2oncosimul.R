pwd2 <- getwd()
setwd("../../data/")
source("toy_datasets.R")
setwd("../examples/")
source("access_genots_from_oncosimul.R")
setwd("../code_from_what_genotype_next/")
source("code-all-methods-minimal.R")
setwd(pwd2)
rm(pwd2)

compare_HESBCN_cpm2tm <- function(data){
    out <- all_methods_2_trans_mat(data)
    out$HESBCN_model$Relation <- sapply(out$HESBCN_model$To
        , function(x) out$HESBCN_parent_set[x])
    out_onco <- cpm2tm(out$HESBCN_model)

    order1 <- sort(rownames(out$HESBCN_trans_mat), index.return = TRUE)$ix
    order2 <- sort(rownames(out_onco$transition_matrix), index.return = TRUE)$ix

    out1 <- out$HESBCN_trans_mat[order1, order1]
    out2 <- out_onco$transition_matrix[order2, order2]

    sum(round(out1 - out2, 6)) == 0
}

test_that("HESBCN gives the same results as OncoSimul", {
   for (i in names(all_examples)){
        tmp_dataset <- all_examples[[i]]
        compare_methods <- compare_HESBCN_cpm2tm(tmp_dataset)
        expect_equal(compare_methods, TRUE)
    }
})