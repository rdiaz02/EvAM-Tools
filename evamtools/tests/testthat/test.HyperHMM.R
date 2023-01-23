## Testing HyperHMM

t1 <- Sys.time()
test_that("Testing HyperHMM integration in EvAM-Tools", {
    data <- cbind(c(1,1,0,0), c(1,0,0,0),c(0,0,1,0),c(0,1,1,0))
    colnames(data) <- c('A','B','C','D')
    rownames(data) <- c('A','B','C','D')
    out <- suppressMessages(evam(data, methods = 'HyperHMM'))
    expect_true(inherits(out$HyperHMM_trans_mat, 'Matrix'))
    expect_true(all(slotNames(out$HyperHMM_trans_mat)==c('i','p','Dim','Dimnames','x','factors')))
    lista<-lapply(1:length(out$HyperHMM_features),
                  function(y) combn(out$HyperHMM_features,y,simplify=FALSE))
    lista<-c('WT',unlist(lapply(lista, 
                                function(y) (lapply(y, function(z) toString(z))))))
    expect_true(all(out$HyperHMM_trans_mat@Dim==c(length(lista), length(lista))))
})

cat("\n Done test.HyperHMM.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")