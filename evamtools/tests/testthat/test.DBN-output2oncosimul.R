test_that("DBN gives the same results as OncoSimul", {

    compare_DBN_cpm2tm <- function(codename) {
        ## FIXME: don't do this. 
        out <- cpm_output[[codename]]
        out_onco <- evamtools:::cpm2tm(out$DBN_model)

        order1 <- sort(rownames(out$DBN_trans_mat), index.return = TRUE)$ix
        order2 <- sort(rownames(out_onco$transition_matrix), index.return = TRUE)$ix

        ordered_computed_trm <- out$DBN_trans_mat[order1, order1]
        ordered_trm_onco <- out_onco$transition_matrix[order2, order2]

        expect_equal(ordered_computed_trm, ordered_trm_onco)
    }

    ## FIXME: don't do this. 
    for (i in names(examples_csd[["csd"]])[2:11]) {
        compare_methods <- compare_DBN_cpm2tm(i)
    }
})
cat("\n Done test.DBN-output2oncosimul.R \n")


