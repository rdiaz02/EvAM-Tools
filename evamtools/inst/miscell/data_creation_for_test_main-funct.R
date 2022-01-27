## How we run former all_methods_2_trans_mat to generate output so that we can
## refactor the function. Note that we have to unset an env. var. to allow MCCBN
## to run. Also need to load relations (MCCBN). Beware that CBN has a random
## component that cannot be fixed using set.seed. And the same thing happens with
## MCCBN.


library(relations)
Sys.unsetenv("_R_CHECK_LENGTH_1_LOGIC2_")
set.seed(1)
test_main_input_data <- all_methods_2_trans_mat(input_data,
                                                methods = c("CBN", "OT",
                                                            "MHN", "HESBCN",
                                                            "DBN", "MCCBN"))

set.seed(111)
test_main_Dat1 <- all_methods_2_trans_mat(Dat1,
                                          methods = c("CBN", "OT",
                                                      "MHN", "HESBCN",
                                                      "DBN", "MCCBN"))

set.seed(653)
test_main_Dat2 <- all_methods_2_trans_mat(Dat2,
                                          methods = c("CBN", "OT",
                                                      "MHN", "HESBCN",
                                                      "DBN", "MCCBN"))

set.seed(98)
test_main_Dat3 <- all_methods_2_trans_mat(Dat3,
                                          methods = c("CBN", "OT",
                                                      "MHN", "HESBCN",
                                                      "DBN", "MCCBN"))
set.seed(44)
test_main_Dat4 <- all_methods_2_trans_mat(Dat4,
                                          methods = c("CBN", "OT",
                                                      "MHN", "HESBCN",
                                                      "DBN", "MCCBN"))


runs_for_refactor_main_funct<- list(test_main_input_data = test_main_input_data,
                                    test_main_Dat1 = test_main_Dat1,
                                    test_main_Dat2 = test_main_Dat2,
                                    test_main_Dat3 = test_main_Dat3,
                                    test_main_Dat4 = test_main_Dat4
                                    )

save(file = "runs_for_refactor_main_funct.RData",
     runs_for_refactor_main_funct)


## If you run them again, they should be identical to the above
## in some of the components. Not necessarily CBN or MCCBN.

set.seed(1)
t0 <- all_methods_2_trans_mat(input_data,
                              methods = c("CBN", "OT",
                                          "MHN", "HESBCN",
                                          "DBN", "MCCBN"))
expect_identical(test_main_input_data, t0)

set.seed(111)
t1 <- all_methods_2_trans_mat(Dat1,
                              methods = c("CBN", "OT",
                                          "MHN", "HESBCN",
                                          "DBN", "MCCBN"))
expect_identical(test_main_Dat1, t1)

set.seed(653)
t2 <- all_methods_2_trans_mat(Dat2,
                              methods = c("CBN", "OT",
                                          "MHN", "HESBCN",
                                          "DBN", "MCCBN"))
expect_identical(test_main_Dat2, t2)

set.seed(98)
t3 <- all_methods_2_trans_mat(Dat3,
                              methods = c("CBN", "OT",
                                          "MHN", "HESBCN",
                                          "DBN", "MCCBN"))
expect_identical(test_main_Dat3, t3)

set.seed(44)
t4 <- all_methods_2_trans_mat(Dat4,
                              methods = c("CBN", "OT",
                                          "MHN", "HESBCN",
                                          "DBN", "MCCBN"))
expect_identical(test_main_Dat4, t4)

