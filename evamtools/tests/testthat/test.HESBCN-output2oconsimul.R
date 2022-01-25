## FIXME: commented until fixed
test_that("HESBCN gives the same results as OncoSimul", {

    compare_HESBCN_cpm2tm <- function(codename) {
        out <- cpm_output[[codename]]
        out$HESBCN_model$Relation <- sapply(
            out$HESBCN_model$To,
            function(x) out$HESBCN_parent_set[x]
        )
        out_onco <- evamtools:::cpm2tm(out$HESBCN_model)

        order1 <- sort(rownames(out$HESBCN_trans_mat), index.return = TRUE)$ix
        order2 <- sort(rownames(out_onco$transition_matrix), index.return = TRUE)$ix

        ordered_computed_trm <- out$HESBCN_trans_mat[order1, order1]
        ordered_trm_onco <- out_onco$transition_matrix[order2, order2]

        expect_equal(ordered_computed_trm, ordered_trm_onco)
    }

    for (i in names(examples_csd[["csd"]])[2:11]) {
        compare_methods <- compare_HESBCN_cpm2tm(i)
    }
})
cat("\n Done test.HESBCN-output2oconsimul.R \n")
