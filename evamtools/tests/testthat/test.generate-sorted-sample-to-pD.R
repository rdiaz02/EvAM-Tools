test_that("generate_sorted_genotypes and sample_to_pD_order same order", {
    ## 6 WT, 3 A, 4 AC, 5 D, 2 AD
    gg <- c("A, C", "C, A", "C, A", "A, C",
            "D, A", "A, D",
            "D", "D", "D", "D", "D",
            "A", "A", "A",
            rep("WT", 6))
    expect_out_gg_counts <- c(6, 3, 0, 4, 5, 2, 0, 0)
    expect_out_gsg <- c("WT", "A", "C", "A, C", "D", "A, D", "C, D", "A, C, D")
    for(i in 1:50) {
        out_pd <- evamtools:::sample_to_pD_order(sample(gg),
                                                 3,
                                                 sample(c("A", "C", "D")))
        out_gsg <- evamtools:::generate_sorted_genotypes(3,
                                                         sample(c("A", "C", "D")))
        expect_equal(out_pd, expect_out_gg_counts)
        expect_equal(out_gsg, expect_out_gsg)
    }


    gg2 <- c("WT", "A", "A",
             "B", "B", "B",
             rep("C", 4),
             rep("A, B", 5),
             rep("A, C", 6),
             rep("B, C", 7),
             rep("A, B, C", 8))

    expect_out_gg2_counts <- c(1, 2, 3, 5, 4, 6, 7, 8)
    expect_out_gsg2 <- c("WT", "A", "B", "A, B", "C", "A, C", "B, C", "A, B, C")
    for(i in 1:50) {
        out_pd2 <- evamtools:::sample_to_pD_order(sample(gg2),
                                                  3,
                                                  sample(c("A", "C", "B")))
        out_gsg2 <- evamtools:::generate_sorted_genotypes(3,
                                                          sample(c("A", "C", "B")))
        expect_equal(out_pd2, expect_out_gg2_counts)
        expect_equal(out_gsg2, expect_out_gsg2)
    }



    gg3 <- c(rep("WT", 23),
             rep("ZE", 15),
             rep("b_mut", 12),
             rep("F_NRAS", 4),
             rep("ZE, b_mut", 3), rep("b_mut, ZE", 2), ## 5
             rep("ZE, F_NRAS", 2), rep("F_NRAS, ZE", 4), ## 6
             rep("b_mut, F_NRAS", 5), rep("F_NRAS, b_mut", 2), ## 7
             rep("ZE, b_mut, F_NRAS", 2),
             rep("ZE, F_NRAS, b_mut", 2),
             "b_mut, ZE, F_NRAS",
             "b_mut, F_NRAS, ZE",
             "F_NRAS, b_mut, ZE", ## 8 of these
             "F_NRAS, ZE, b_mut",
             "colorin3, b_mut",
             "colorin4, F_NRAS, ZE, colorin3",
             "F_NRAS, colorin4",
             "colorin4, F_NRAS")

    expect_out_gsg3 <- c(
        "WT"
      , "b_mut"
      , "colorin3"
      , "b_mut, colorin3"
      , "colorin4"
      , "b_mut, colorin4"
      , "colorin3, colorin4"
      , "b_mut, colorin3, colorin4"
      , "F_NRAS"
      , "b_mut, F_NRAS"
      , "colorin3, F_NRAS"
      , "b_mut, colorin3, F_NRAS"
      , "colorin4, F_NRAS"
      , "b_mut, colorin4, F_NRAS"
      , "colorin3, colorin4, F_NRAS"
      , "b_mut, colorin3, colorin4, F_NRAS"
      , "no_ta"
      , "b_mut, no_ta"
      , "colorin3, no_ta"
      , "b_mut, colorin3, no_ta"
      , "colorin4, no_ta"
      , "b_mut, colorin4, no_ta"
      , "colorin3, colorin4, no_ta"
      , "b_mut, colorin3, colorin4, no_ta"
      , "F_NRAS, no_ta"
      , "b_mut, F_NRAS, no_ta"
      , "colorin3, F_NRAS, no_ta"
      , "b_mut, colorin3, F_NRAS, no_ta"
      , "colorin4, F_NRAS, no_ta"
      , "b_mut, colorin4, F_NRAS, no_ta"
      , "colorin3, colorin4, F_NRAS, no_ta"
      , "b_mut, colorin3, colorin4, F_NRAS, no_ta"
      , "ZE"
      , "b_mut, ZE"
      , "colorin3, ZE"
      , "b_mut, colorin3, ZE"
      , "colorin4, ZE"
      , "b_mut, colorin4, ZE"
      , "colorin3, colorin4, ZE"
      , "b_mut, colorin3, colorin4, ZE"
      , "F_NRAS, ZE"
      , "b_mut, F_NRAS, ZE"
      , "colorin3, F_NRAS, ZE"
      , "b_mut, colorin3, F_NRAS, ZE"
      , "colorin4, F_NRAS, ZE"
      , "b_mut, colorin4, F_NRAS, ZE"
      , "colorin3, colorin4, F_NRAS, ZE"
      , "b_mut, colorin3, colorin4, F_NRAS, ZE"
      , "no_ta, ZE"
      , "b_mut, no_ta, ZE"
      , "colorin3, no_ta, ZE"
      , "b_mut, colorin3, no_ta, ZE"
      , "colorin4, no_ta, ZE"
      , "b_mut, colorin4, no_ta, ZE"
      , "colorin3, colorin4, no_ta, ZE"
      , "b_mut, colorin3, colorin4, no_ta, ZE"
      , "F_NRAS, no_ta, ZE"
      , "b_mut, F_NRAS, no_ta, ZE"
      , "colorin3, F_NRAS, no_ta, ZE"
      , "b_mut, colorin3, F_NRAS, no_ta, ZE"
      , "colorin4, F_NRAS, no_ta, ZE"
      , "b_mut, colorin4, F_NRAS, no_ta, ZE"
      , "colorin3, colorin4, F_NRAS, no_ta, ZE"
      , "b_mut, colorin3, colorin4, F_NRAS, no_ta, ZE")




    expect_out_gg3_counts <- c(
        23, #"WT"
        12, # "b_mut"
        0, # "colorin3"
        1, # "b_mut, colorin3"
        0, # "colorin4"
        0, # "b_mut, colorin4"
        0, # "colorin3, colorin4"
        0, # "b_mut, colorin3, colorin4"
        4, # "F_NRAS"
        7, # "b_mut, F_NRAS"
        0, # "colorin3, F_NRAS"
        0, # "b_mut, colorin3, F_NRAS"
        2, # "colorin4, F_NRAS"
        0, # "b_mut, colorin4, F_NRAS"
        0, # "colorin3, colorin4, F_NRAS"
        0, # "b_mut, colorin3, colorin4, F_NRAS"
        0, # "no_ta"
        0, # "b_mut, no_ta"
        0, # "colorin3, no_ta"
        0, # "b_mut, colorin3, no_ta"
        0, # "colorin4, no_ta"
        0, # "b_mut, colorin4, no_ta"
        0, # "colorin3, colorin4, no_ta"
        0, # "b_mut, colorin3, colorin4, no_ta"
        0, # "F_NRAS, no_ta"
        0, # "b_mut, F_NRAS, no_ta"
        0, # "colorin3, F_NRAS, no_ta"
        0, # "b_mut, colorin3, F_NRAS, no_ta"
        0, # "colorin4, F_NRAS, no_ta"
        0, # "b_mut, colorin4, F_NRAS, no_ta"
        0, # "colorin3, colorin4, F_NRAS, no_ta"
        0, # "b_mut, colorin3, colorin4, F_NRAS, no_ta"
        15, # "ZE"
        5, # "b_mut, ZE"
        0, # "colorin3, ZE"
        0, # "b_mut, colorin3, ZE"
        0, # "colorin4, ZE"
        0, # "b_mut, colorin4, ZE"
        0, # "colorin3, colorin4, ZE"
        0, # "b_mut, colorin3, colorin4, ZE"
        6, # "F_NRAS, ZE"
        8, # "b_mut,  F_NRAS,  ZE"
        0, # "colorin3,  F_NRAS,  ZE"
        0, # "b_mut,  colorin3,  F_NRAS,  ZE"
        0, # "colorin4,  F_NRAS,  ZE"
        0, # "b_mut,  colorin4,  F_NRAS,  ZE"
        1, # "colorin3,  colorin4,  F_NRAS,  ZE"
        0, # "b_mut,  colorin3,  colorin4,  F_NRAS,  ZE"
        0, # "no_ta,  ZE"
        0, # "b_mut,  no_ta,  ZE"
        0, # "colorin3,  no_ta,  ZE"
        0, # "b_mut,  colorin3,  no_ta,  ZE"
        0, # "colorin4,  no_ta,  ZE"
        0, # "b_mut,  colorin4,  no_ta,  ZE"
        0, # "colorin3,  colorin4,  no_ta,  ZE"
        0, # "b_mut,  colorin3,  colorin4,  no_ta,  ZE"
        0, # "F_NRAS,  no_ta,  ZE"
        0, # "b_mut,  F_NRAS,  no_ta,  ZE"
        0, # "colorin3,  F_NRAS,  no_ta,  ZE"
        0, # "b_mut,  colorin3,  F_NRAS,  no_ta,  ZE"
        0, # "colorin4,  F_NRAS,  no_ta,  ZE"
        0, # "b_mut,  colorin4,  F_NRAS,  no_ta,  ZE"
        0, # "colorin3,  colorin4,  F_NRAS,  no_ta,  ZE"
        0 # "b_mut,  colorin3,  colorin4,  F_NRAS,  no_ta,  ZE"
    )



    
    for(i in 1:50) {
        genes <- c("ZE", "b_mut", "F_NRAS", "colorin3", "no_ta", "colorin4")
        out_pd3 <- evamtools:::sample_to_pD_order(sample(gg3),
                                                  6,
                                                  sample(genes))
        out_gsg3 <- evamtools:::generate_sorted_genotypes(6,
                                                          sample(genes))
        expect_equal(out_pd3, expect_out_gg3_counts)
        expect_equal(out_gsg3, expect_out_gsg3)
    }
})
cat("\n Done test.generate-sorted-sample-to-pD. \n")
