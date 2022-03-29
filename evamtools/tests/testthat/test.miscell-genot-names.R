## Test miscellaneous functions for dealing with genotype names
## and ordering vectors of genotypes

t1 <- Sys.time()

test_that("canonicalize", {
    www <- c("B" = 8,
            "A,B,C" = 10,
            "WT" = 2,
            "A" = 3,
            "B,E,C,A" = 103,
            "E,    D  ,  A   ,B" = 98,
            " U#/!=?¿%&8  ,  F_#$@  ,  A  " = 5)

    expect_identical(canonicalize_genotype_names(names(www)),
                     c("B", "A, B, C", "WT", "A", "A, B, C, E",
                       "A, B, D, E", "A, F_#$@, U#/!=?¿%&8")
                     )
    
})



test_that("reorder_to_pD", {
    www <- c("B" = 8,
            "A,B,C" = 10,
            "WT" = 2,
            "A" = 3,
            "B,E,C,A" = 103,
            "E,    D  ,  A   ,B" = 98,
            " U#/!=?¿%&8  ,  F_#$@  ,  A  " = 5)


    sgu <- generate_pD_sorted_genotypes(n_genes = 7,
                                     gene_names = c("A", "B", "C",
                                                    "D", "E",
                                                    "U#/!=?¿%&8",
                                                    "F_#$@"))

    sguv <- rep(NA, length(sgu))
    names(sguv) <- sgu
    ## Yes, pass by position, not matching names on purpose for test
    sguv[c(3, 8, 1, 2, 24, 28, 98)] <- unname(www)
    
    expect_identical(reorder_to_pD(www),
                     sguv)


     www2 <- c("B" = .8,
            "A,B,C" = .10,
            "WT" = .22,
            "A" = .33,
            "B,E,C,A" = .103,
            "E,    D  ,  A   ,B" = .98,
            " U#/!=?¿%&8  ,  F_#$@  ,  A  " = .5,
            "    B  ,  F_#$@   , D " = 0.12
            )

    sguv2 <- rep(NA, length(sgu))
    names(sguv2) <- sgu
    sguv2[c(3, 8, 1, 2, 24, 28, 98, 43)] <- unname(www2)
    
})

cat("\n Done test.miscell-genot-names.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
