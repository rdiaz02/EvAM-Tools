data(every_which_way_data)
Dat1 <- every_which_way_data[[16]][1:40, 2:6]
## null <-
##     all_methods_2_trans_mat(Dat1,
##                             do_MCCBN = requireNamespace("mccbn", quietly = TRUE))
null <- all_methods_2_trans_mat(Dat1, do_MCCBN = FALSE)
cat("\n Done test.exercise-main-funct.R \n")
