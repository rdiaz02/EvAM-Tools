## Testing HyperHMM

t1 <- Sys.time()
test_that("Testing HyperHMM integration in EvAM-Tools", {
    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:6]
    out <- suppressMessages(evam(Dat1, methods = 'HyperHMM'))
    expect_true(inherits(out$HyperHMM_trans_mat, 'Matrix'))
    expect_true(all(slotNames(out$HyperHMM_trans_mat)==c('i','p','Dim','Dimnames','x','factors')))        
})

cat("\n Done test.HyperHMM.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")