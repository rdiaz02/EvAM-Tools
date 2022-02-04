test_that("Minimal test: we can run", {
    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:6]
    null <- expect_output(suppressMessages(evam(Dat1,
                                  methods = c("CBN", "OT",
                                              "MHN", "HESBCN",
                                              "DBN"))))
})

cat("\n Done test.exercise-main-funct.R \n")
