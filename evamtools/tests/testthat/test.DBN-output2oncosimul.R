test_that("DBN gives the same results as OncoSimul", {


    ## FIXME here, explain, EXPLICITLY and FULLY, where the data
    ## are  (i.e., where the RData is) and how and where the R file
    ## that generate the RData is.
    ## We need to make sure we can examine what are the scenarios tested. Right
    ## now, this is opaque.

    ## FIXME: what is this really testing? What function of evamtools is being
    ## tested? I only see a call to cpm2tm. But not other calls to other
    ## functions. For instance, to cpm_access_genots_paths_w_simplified or to
    ## evam or whatever. In other words, if we were to make changes in functions
    ## in evamtools, and introduce bugs, would this detect anything?

    ## FIXME: the key is NOT testing cpm2tm, but testing the code used in
    ## evamtools for the analysis of data.

    ## FIXME: I do not see this function really testing anything.
    
    ## FIXME: isn't this using global variables? What is cpm_output?
    
    compare_DBN_cpm2tm <- function(codename) {
        out <- cpm_output[[codename]]
        out_onco <- evamtools:::cpm2tm(out$DBN_model)

        order1 <- sort(rownames(out$DBN_trans_mat), index.return = TRUE)$ix
        order2 <- sort(rownames(out_onco$transition_matrix), index.return = TRUE)$ix

        ordered_computed_trm <- out$DBN_trans_mat[order1, order1]
        ordered_trm_onco <- out_onco$transition_matrix[order2, order2]

        expect_equal(ordered_computed_trm, ordered_trm_onco)
    }

    ## FIXME: referring to a list by index instead of name is often confusing in
    ## these cases. Why 2? Why not 1? Why not beyond 11? Etc. And what are 2 to
    ## 11?

    for (i in names(examples_csd[["csd"]])[2:11]) {
        compare_methods <- compare_DBN_cpm2tm(i)
    }
})
cat("\n Done test.DBN-output2oncosimul.R \n")


