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



test_that("reorder_to_standard_order minimal tests", {
    ## This actually happened to break code
    v1 <- c(WT = "WT", ATP2B2 = "ATP2B2", PIK3CA = "PIK3CA", PNPLA3 = "PNPLA3",
            RB1 = "RB1", TP53 = "TP53", TRIM6 = "TRIM6", `ATP2B2, PIK3CA` = "ATP2B2, PIK3CA",
            `ATP2B2, PNPLA3` = "ATP2B2, PNPLA3", `ATP2B2, RB1` = "ATP2B2, RB1",
            `ATP2B2, TP53` = "ATP2B2, TP53", `ATP2B2, TRIM6` = "ATP2B2, TRIM6",
            `PIK3CA, PNPLA3` = "PIK3CA, PNPLA3", `PIK3CA, RB1` = "PIK3CA, RB1",
            `PIK3CA, TP53` = "PIK3CA, TP53", `PIK3CA, TRIM6` = "PIK3CA, TRIM6",
            `ATP2B2, PNPLA3` = "PNPLA3, ATP2B2", `PIK3CA, PNPLA3` = "PNPLA3, PIK3CA",
            `PNPLA3, RB1` = "PNPLA3, RB1", `PNPLA3, TP53` = "PNPLA3, TP53",
            `PNPLA3, TRIM6` = "PNPLA3, TRIM6", `RB1, TP53` = "RB1, TP53",
            `RB1, TRIM6` = "RB1, TRIM6", `ATP2B2, TP53` = "TP53, ATP2B2",
            `PIK3CA, TP53` = "TP53, PIK3CA", `RB1, TP53` = "TP53, RB1", `TP53, TRIM6` = "TP53, TRIM6",
            `ATP2B2, TRIM6` = "TRIM6, ATP2B2", `ATP2B2, PIK3CA, PNPLA3` = "ATP2B2, PIK3CA, PNPLA3",
            `ATP2B2, PIK3CA, RB1` = "ATP2B2, PIK3CA, RB1", `ATP2B2, PIK3CA, TP53` = "ATP2B2, PIK3CA, TP53",
            `ATP2B2, PIK3CA, TRIM6` = "ATP2B2, PIK3CA, TRIM6", `ATP2B2, PNPLA3, RB1` = "ATP2B2, PNPLA3, RB1",
            `ATP2B2, PNPLA3, TP53` = "ATP2B2, PNPLA3, TP53", `ATP2B2, PNPLA3, TRIM6` = "ATP2B2, PNPLA3, TRIM6",
            `ATP2B2, RB1, TP53` = "ATP2B2, RB1, TP53", `ATP2B2, RB1, TRIM6` = "ATP2B2, RB1, TRIM6",
            `ATP2B2, TP53, TRIM6` = "ATP2B2, TP53, TRIM6", `PIK3CA, PNPLA3, RB1` = "PIK3CA, PNPLA3, RB1",
            `PIK3CA, PNPLA3, TP53` = "PIK3CA, PNPLA3, TP53", `PIK3CA, PNPLA3, TRIM6` = "PIK3CA, PNPLA3, TRIM6",
            `PIK3CA, RB1, TP53` = "PIK3CA, RB1, TP53", `PIK3CA, RB1, TRIM6` = "PIK3CA, RB1, TRIM6",
            `PIK3CA, TP53, TRIM6` = "PIK3CA, TP53, TRIM6", `PIK3CA, PNPLA3, RB1` = "PNPLA3, PIK3CA, RB1",
            `PNPLA3, RB1, TP53` = "PNPLA3, RB1, TP53", `PNPLA3, RB1, TRIM6` = "PNPLA3, RB1, TRIM6",
            `ATP2B2, PNPLA3, TP53` = "PNPLA3, TP53, ATP2B2", `PIK3CA, PNPLA3, TP53` = "PNPLA3, TP53, PIK3CA",
            `PNPLA3, RB1, TP53` = "PNPLA3, TP53, RB1", `PNPLA3, TP53, TRIM6` = "PNPLA3, TP53, TRIM6",
            `ATP2B2, PNPLA3, TRIM6` = "PNPLA3, TRIM6, ATP2B2", `PIK3CA, PNPLA3, TRIM6` = "PNPLA3, TRIM6, PIK3CA",
            `RB1, TP53, TRIM6` = "RB1, TP53, TRIM6", `ATP2B2, PIK3CA, TP53` = "TP53, ATP2B2, PIK3CA",
            `ATP2B2, RB1, TP53` = "TP53, ATP2B2, RB1", `PIK3CA, RB1, TP53` = "TP53, PIK3CA, RB1",
            `ATP2B2, TP53, TRIM6` = "TP53, TRIM6, ATP2B2", `PIK3CA, TP53, TRIM6` = "TP53, TRIM6, PIK3CA",
            `RB1, TP53, TRIM6` = "TP53, TRIM6, RB1", `ATP2B2, PIK3CA, PNPLA3, RB1` = "ATP2B2, PIK3CA, PNPLA3, RB1",
            `ATP2B2, PIK3CA, PNPLA3, TP53` = "ATP2B2, PIK3CA, PNPLA3, TP53",
            `ATP2B2, PIK3CA, PNPLA3, TRIM6` = "ATP2B2, PIK3CA, PNPLA3, TRIM6",
            `ATP2B2, PIK3CA, RB1, TP53` = "ATP2B2, PIK3CA, RB1, TP53", `ATP2B2, PIK3CA, RB1, TRIM6` = "ATP2B2, PIK3CA, RB1, TRIM6",
            `ATP2B2, PIK3CA, TP53, TRIM6` = "ATP2B2, PIK3CA, TP53, TRIM6",
            `ATP2B2, PNPLA3, RB1, TP53` = "ATP2B2, PNPLA3, RB1, TP53", `ATP2B2, PNPLA3, RB1, TRIM6` = "ATP2B2, PNPLA3, RB1, TRIM6",
            `ATP2B2, PNPLA3, TP53, TRIM6` = "ATP2B2, PNPLA3, TP53, TRIM6",
            `ATP2B2, RB1, TP53, TRIM6` = "ATP2B2, RB1, TP53, TRIM6", `PIK3CA, PNPLA3, RB1, TP53` = "PIK3CA, PNPLA3, RB1, TP53",
            `PIK3CA, PNPLA3, RB1, TRIM6` = "PIK3CA, PNPLA3, RB1, TRIM6",
            `PIK3CA, PNPLA3, TP53, TRIM6` = "PIK3CA, PNPLA3, TP53, TRIM6",
            `PIK3CA, RB1, TP53, TRIM6` = "PIK3CA, RB1, TP53, TRIM6", `ATP2B2, PIK3CA, PNPLA3, RB1` = "PNPLA3, ATP2B2, PIK3CA, RB1",
            `PNPLA3, RB1, TP53, TRIM6` = "PNPLA3, RB1, TP53, TRIM6", `ATP2B2, PIK3CA, PNPLA3, TP53` = "PNPLA3, TP53, ATP2B2, PIK3CA",
            `ATP2B2, PNPLA3, RB1, TP53` = "PNPLA3, TP53, ATP2B2, RB1", `PIK3CA, PNPLA3, RB1, TP53` = "PNPLA3, TP53, PIK3CA, RB1",
            `PIK3CA, PNPLA3, TP53, TRIM6` = "PNPLA3, TP53, TRIM6, PIK3CA",
            `PNPLA3, RB1, TP53, TRIM6` = "PNPLA3, TP53, TRIM6, RB1", `ATP2B2, PIK3CA, PNPLA3, TRIM6` = "PNPLA3, TRIM6, ATP2B2, PIK3CA",
            `PIK3CA, PNPLA3, RB1, TRIM6` = "PNPLA3, TRIM6, PIK3CA, RB1",
            `ATP2B2, PIK3CA, RB1, TP53` = "TP53, ATP2B2, PIK3CA, RB1", `ATP2B2, PIK3CA, TP53, TRIM6` = "TP53, TRIM6, ATP2B2, PIK3CA",
            `ATP2B2, RB1, TP53, TRIM6` = "TP53, TRIM6, ATP2B2, RB1", `PIK3CA, RB1, TP53, TRIM6` = "TP53, TRIM6, PIK3CA, RB1",
            `ATP2B2, PIK3CA, PNPLA3, RB1, TP53` = "ATP2B2, PIK3CA, PNPLA3, RB1, TP53",
            `ATP2B2, PIK3CA, PNPLA3, RB1, TRIM6` = "ATP2B2, PIK3CA, PNPLA3, RB1, TRIM6",
            `ATP2B2, PIK3CA, PNPLA3, TP53, TRIM6` = "ATP2B2, PIK3CA, PNPLA3, TP53, TRIM6",
            `ATP2B2, PIK3CA, RB1, TP53, TRIM6` = "ATP2B2, PIK3CA, RB1, TP53, TRIM6",
            `ATP2B2, PNPLA3, RB1, TP53, TRIM6` = "ATP2B2, PNPLA3, RB1, TP53, TRIM6",
            `PIK3CA, PNPLA3, RB1, TP53, TRIM6` = "PIK3CA, PNPLA3, RB1, TP53, TRIM6",
            `ATP2B2, PIK3CA, PNPLA3, RB1, TP53` = "PNPLA3, TP53, ATP2B2, PIK3CA, RB1",
            `ATP2B2, PIK3CA, PNPLA3, TP53, TRIM6` = "PNPLA3, TP53, TRIM6, ATP2B2, PIK3CA",
            `ATP2B2, PNPLA3, RB1, TP53, TRIM6` = "PNPLA3, TP53, TRIM6, ATP2B2, RB1",
            `PIK3CA, PNPLA3, RB1, TP53, TRIM6` = "PNPLA3, TP53, TRIM6, PIK3CA, RB1",
            `ATP2B2, PIK3CA, RB1, TP53, TRIM6` = "TP53, TRIM6, ATP2B2, PIK3CA, RB1",
            `ATP2B2, PIK3CA, PNPLA3, RB1, TP53, TRIM6` = "ATP2B2, PIK3CA, PNPLA3, RB1, TP53, TRIM6",
            `ATP2B2, PIK3CA, PNPLA3, RB1, TP53, TRIM6` = "PNPLA3, TP53, TRIM6, ATP2B2, PIK3CA, RB1"
            )

    expect_error(reorder_to_standard_order(v1),
                 "x has duplicated genotype names; at least these positions 17, 18")

    ## Same thing, before ordering. This came from a bug in
    ## probs_from_HT
    v1o <- c(WT = 0.0045379746835443, ATP2B2 = 0.00182784810126582, PIK3CA = 4.30379746835443e-05,
             PNPLA3 = 0.000330379746835443, RB1 = 8.86075949367089e-05, TP53 = 0.328784810126582,
             TRIM6 = 2.78481012658228e-05, `ATP2B2, PIK3CA` = 0, `ATP2B2, PNPLA3` = 0,
             `ATP2B2, RB1` = 0, `ATP2B2, TP53` = 0, `ATP2B2, TRIM6` = 0, `PIK3CA, PNPLA3` = 0,
             `PIK3CA, RB1` = 0, `PIK3CA, TP53` = 0, `PIK3CA, TRIM6` = 0, `PNPLA3, RB1` = 0.000230379746835443,
             `PNPLA3, TP53` = 0.00421392405063291, `PNPLA3, TRIM6` = 3.16455696202532e-05,
             `RB1, TP53` = 0, `RB1, TRIM6` = 0, `TP53, TRIM6` = 0.0413405063291139,
             `ATP2B2, PIK3CA, PNPLA3` = 0, `ATP2B2, PIK3CA, RB1` = 0, `ATP2B2, PIK3CA, TP53` = 0,
             `ATP2B2, PIK3CA, TRIM6` = 0, `ATP2B2, PNPLA3, RB1` = 0, `ATP2B2, PNPLA3, TP53` = 0,
             `ATP2B2, PNPLA3, TRIM6` = 0, `ATP2B2, RB1, TP53` = 0, `ATP2B2, RB1, TRIM6` = 0,
             `ATP2B2, TP53, TRIM6` = 0, `PIK3CA, PNPLA3, RB1` = 0, `PIK3CA, PNPLA3, TP53` = 0,
             `PIK3CA, PNPLA3, TRIM6` = 0, `PIK3CA, RB1, TP53` = 0, `PIK3CA, RB1, TRIM6` = 0,
             `PIK3CA, TP53, TRIM6` = 0, `PNPLA3, RB1, TP53` = 0, `PNPLA3, RB1, TRIM6` = 0,
             `PNPLA3, TP53, TRIM6` = 0.000330379746835443, `RB1, TP53, TRIM6` = 0,
             `ATP2B2, PIK3CA, PNPLA3, RB1` = 0, `ATP2B2, PIK3CA, PNPLA3, TP53` = 0,
             `ATP2B2, PIK3CA, PNPLA3, TRIM6` = 0, `ATP2B2, PIK3CA, RB1, TP53` = 0,
             `ATP2B2, PIK3CA, RB1, TRIM6` = 0, `ATP2B2, PIK3CA, TP53, TRIM6` = 0,
             `ATP2B2, PNPLA3, RB1, TP53` = 0, `ATP2B2, PNPLA3, RB1, TRIM6` = 0,
             `ATP2B2, PNPLA3, TP53, TRIM6` = 0, `ATP2B2, RB1, TP53, TRIM6` = 0,
             `PIK3CA, PNPLA3, RB1, TP53` = 0, `PIK3CA, PNPLA3, RB1, TRIM6` = 0,
             `PIK3CA, PNPLA3, TP53, TRIM6` = 0, `PIK3CA, RB1, TP53, TRIM6` = 0,
             `PNPLA3, RB1, TP53, TRIM6` = 0, `ATP2B2, PIK3CA, PNPLA3, RB1, TP53` = 0,
             `ATP2B2, PIK3CA, PNPLA3, RB1, TRIM6` = 0, `ATP2B2, PIK3CA, PNPLA3, TP53, TRIM6` = 0,
             `ATP2B2, PIK3CA, RB1, TP53, TRIM6` = 0, `ATP2B2, PNPLA3, RB1, TP53, TRIM6` = 0,
             `PIK3CA, PNPLA3, RB1, TP53, TRIM6` = 0, `ATP2B2, PIK3CA, PNPLA3, RB1, TP53, TRIM6` = 0,
             `TP53, RB1` = 0.0247569620253165, `TP53, PIK3CA` = 0.222087341772152,
             `TP53, PIK3CA, RB1` = 0.0137746835443038, `TP53, ATP2B2` = 0.00796582278481013,
             `TP53, ATP2B2, RB1` = 0.0121632911392405, `TP53, ATP2B2, PIK3CA` = 0.00261772151898734,
             `TP53, TRIM6, RB1` = 0.00459873417721519, `TP53, TRIM6, PIK3CA` = 0.0911405063291139,
             `TP53, TRIM6, PIK3CA, RB1` = 0.0166481012658228, `TP53, TRIM6, ATP2B2` = 0.00050506329113924,
             `TP53, TRIM6, ATP2B2, PIK3CA` = 0.00382278481012658, `TP53, TRIM6, ATP2B2, PIK3CA, RB1` = 0.00259240506329114,
             `PNPLA3, PIK3CA` = 0.000684810126582279, `PNPLA3, PIK3CA, RB1` = 0.00306708860759494,
             `PNPLA3, ATP2B2` = 0.000691139240506329, `PNPLA3, TP53, RB1` = 0.0111063291139241,
             `PNPLA3, TP53, PIK3CA` = 0.0166341772151899, `PNPLA3, TP53, PIK3CA, RB1` = 0.0394721518987342,
             `PNPLA3, TP53, ATP2B2` = 0.000460759493670886, `PNPLA3, TP53, ATP2B2, PIK3CA, RB1` = 0.00512025316455696,
             `PNPLA3, TP53, TRIM6, RB1` = 0.00230506329113924, `PNPLA3, TP53, TRIM6, PIK3CA` = 0.0169101265822785,
             `PNPLA3, TP53, TRIM6, PIK3CA, RB1` = 0.0984810126582279, `PNPLA3, TP53, TRIM6, ATP2B2, PIK3CA, RB1` = 0.00848227848101266,
             `TP53, ATP2B2, PIK3CA, RB1` = 0.00534050632911392, `TP53, TRIM6, ATP2B2, RB1` = 0.000813924050632911,
             `PNPLA3, ATP2B2, PIK3CA, RB1` = 0.00103164556962025, `PNPLA3, TP53, ATP2B2, RB1` = 0.00200759493670886,
             `PNPLA3, TP53, ATP2B2, PIK3CA` = 0.00183291139240506, `PNPLA3, TRIM6, PIK3CA, RB1` = 0.000346835443037975,
             `PNPLA3, TRIM6, PIK3CA` = 4.81012658227848e-05, `TRIM6, ATP2B2` = 0.000167088607594937,
             `PNPLA3, TP53, TRIM6, ATP2B2, RB1` = 0.000279746835443038, `PNPLA3, TP53, TRIM6, ATP2B2, PIK3CA` = 0.000231645569620253,
             `PNPLA3, TRIM6, ATP2B2` = 2.15189873417722e-05, `PNPLA3, TRIM6, ATP2B2, PIK3CA` = 2.53164556962025e-06
             )

    expect_error(reorder_to_standard_order(v1o),
                 "x has duplicated genotype names; at least these positions 65, 66,")

    ## Minimal test cases
    v2 <- c("A" = 3, "A, B" = 8, "H, C, C" = 14, "WT" = 9, "H" = 1)
    expect_error(reorder_to_standard_order(v2),
                 "At least one genotype name not in sorted_genots")

    v3 <- c("A" = 3, "A, B" = 8, "H, C, G" = 14, "WT" = 9, "H" = 1)
    v3o <- reorder_to_standard_order(v3)
    expect_true(length(v3o) == 32)
    ## The third element of v3 is not in canonical form, so there should be an NA
    ## when indexing with those names in the ordered and canonicalized v3o
    expect_true(identical(unname(v3o[names(v3)]), unname(c(v3[c(1, 2)], NA, v3[c(4, 5)]))))
    expect_true(v3o["C, G, H"] == 14)
})


cat("\n Done test.miscell-genot-names.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
