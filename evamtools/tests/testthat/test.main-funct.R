if(FALSE) {
    test_that("evam gives identical results to pre-refactoring function", {
        ## Compare against a few pre-computed examples to allow
        ## for refactoring of the main function.

        ## No need to run on every R CMD check.

        input_data <-
            structure(
                c(0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 
                  1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 0L, 
                  0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 
                  0L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 
                  0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
                  0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                  1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 
                  1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 
                  0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
                  0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                  0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 
                  1L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 
                  0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 
                  0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L),
                .Dim = c(50L, 4L),
                .Dimnames = list(
                    c("1", "2", "3", "4", "5", "6", "7", 
                      "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", 
                      "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", 
                      "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", 
                      "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", 
                      "48", "49", "50"),
                    c("A", "B", "C", "D")
                )
            )


        data(every_which_way_data)

        Dat1 <- every_which_way_data[[2]][1:50, 1:4]
        Dat2 <- every_which_way_data[[4]][1:50, 1:5]
        Dat3 <- every_which_way_data[[20]][1:50, 1:5]
        Dat4 <- every_which_way_data[[16]][1:40, 2:6]

        ## The above were run using a former version of the function, and output was
        ## stored. See details in /misc/data_creation_for_test_main-funct.R

        ## Beware that CBN has a random component that cannot be fixed using
        ## set.seed. And the same thing happens with MCCBN.

        identical_parts <- c(
            "OT_model",                 "OT_f_graph"              ,
            "OT_trans_mat",             "OT_genots_predicted"     ,
            ## "CBN_model",                "CBN_f_graph"             ,
            ## "CBN_trans_mat",            "CBN_td_trans_mat"        ,
            ## "MCCBN_model",              "MCCBN_f_graph"           ,
            ## "MCCBN_trans_mat",          "MCCBN_td_trans_mat"      ,
            "MHN_theta",                "MHN_exp_theta"           ,
            "MHN_transitionRateMatrix", "MHN_trans_mat"           ,
            "MHN_td_trans_mat",         "DBN_model"               ,
            "DBN_likelihood",           "DBN_f_graph"             ,
            "DBN_trans_mat",            "HESBCN_model"            ,
            "HESBCN_parent_set",        "HESBCN_f_graph"          ,
            "HESBCN_trans_mat",         "HESBCN_td_trans_mat"     ,
            "HyperTraPS_model",         "HyperTraPS_f_graph"      ,
            "HyperTraPS_trans_mat",     "HyperTraPS_td_trans_mat" ,
            "csd_data")


        ## Why do we run it if we don't test it? because it should run
        ## And you can still examine output yourself. Often, they will be
        ## equal.
        if(requireNamespace("mccbn", quietly = TRUE)) {
            warning("To use mccbn we need to unset the environment ",
                    "variable _R_CHECK_LENGTH_1_LOGIC2_")
            warning("To use mccbn we need to load library relations")
            Sys.unsetenv("_R_CHECK_LENGTH_1_LOGIC2_")
            library(relations)
            mccbn_M <- "MCCBN"
        } else {
            mccbn_M <- NULL
        }

        set.seed(1)
        test_main_input_data <- evam(input_data,
                                     methods = c("CBN", "OT",
                                                 "MHN", "HESBCN",
                                                 "DBN", mccbn_M))

        set.seed(111)
        test_main_Dat1 <- evam(Dat1,
                               methods = c("CBN", "OT",
                                           "MHN", "HESBCN",
                                           "DBN", mccbn_M))

        set.seed(653)
        test_main_Dat2 <- evam(Dat2,
                               methods = c("CBN", "OT",
                                           "MHN", "HESBCN",
                                           "DBN", mccbn_M))

        set.seed(98)
        test_main_Dat3 <- evam(Dat3,
                               methods = c("CBN", "OT",
                                           "MHN", "HESBCN",
                                           "DBN", mccbn_M))
        set.seed(44)
        test_main_Dat4 <- evam(Dat4,
                               methods = c("CBN", "OT",
                                           "MHN", "HESBCN",
                                           "DBN", mccbn_M))


        data(runs_for_refactor_main_funct)

        expect_equal(test_main_input_data[identical_parts],
                     runs_for_refactor_main_funct[["test_main_input_data"]][identical_parts])

        expect_equal(test_main_Dat1[identical_parts],
                     runs_for_refactor_main_funct[["test_main_Dat1"]][identical_parts])

        expect_equal(test_main_Dat2[identical_parts],
                     runs_for_refactor_main_funct[["test_main_Dat2"]][identical_parts])

        expect_equal(test_main_Dat3[identical_parts],
                     runs_for_refactor_main_funct[["test_main_Dat3"]][identical_parts])

        expect_equal(test_main_Dat4[identical_parts],
                     runs_for_refactor_main_funct[["test_main_Dat4"]][identical_parts])

        
        cat("\n Done test.main-funct.R \n")
    })
}
