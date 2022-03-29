t1 <- Sys.time()



## The possible problems are, of course, the sorting
## See these examples
if(FALSE) {
    lcc <- Sys.setlocale("LC_COLLATE", "C")
    sort(c("a", "A", "z", "Z"))
    sort(c("a", "A", "ba", "ZX", "C-2", "C-100"))

    lcc <- Sys.setlocale("LC_COLLATE", "en_US.UTF-8")
    sort(c("a", "A", "z", "Z"))
    sort(c("a", "A", "ba", "ZX", "C-2", "C-100"))


    stri_sort(c("a", "A", "ba", "ZX", "C-2", "C-100"),
              locale = "en", uppercase_first = FALSE,
              numeric = TRUE)
}



test_that("sorted_minimal_0", {
    expect_out_gsg0 <- c(
        "WT"
      , "ba"
      , "ZX"
      , "ba, ZX")

    out_gsg0 <- generate_pD_sorted_genotypes(2,
                                             c("ZX", "ba"))
    expect_equal(out_gsg0, expect_out_gsg0)
    out_gsg0b <- generate_pD_sorted_genotypes(2,
                                             c("ba", "ZX"))
    expect_equal(out_gsg0b, expect_out_gsg0)
})



test_that("generate_pD_sorted_genotypes and sample_to_pD_order same order", {
    ## 6 WT, 3 A, 4 AC, 5 D, 2 AD
    gg <- c("A, C", "C, A", "C, A", "A, C",
            "D, A", "A, D",
            "D", "D", "D", "D", "D",
            "A", "A", "A",
            rep("WT", 6))
    expect_out_gg_counts <- c(6, 3, 0, 4, 5, 2, 0, 0)
    expect_out_gsg <- c("WT", "A", "C", "A, C", "D", "A, D", "C, D", "A, C, D")
    for(i in 1:50) {
        out_pd <- sample_to_pD_order(sample(gg),
                                                 3,
                                                 sample(c("A", "C", "D")))
        out_gsg <- generate_pD_sorted_genotypes(3,
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
        out_pd2 <- sample_to_pD_order(sample(gg2),
                                                  3,
                                                  sample(c("A", "C", "B")))
        out_gsg2 <- generate_pD_sorted_genotypes(3,
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

    for (i in 1:50) {
        genes <- c("ZE", "b_mut", "F_NRAS", "colorin3", "no_ta", "colorin4")
        out_pd3 <- sample_to_pD_order(sample(gg3),
                                      6,
                                      sample(genes))
        out_gsg3 <- generate_pD_sorted_genotypes(6,
                                                 sample(genes))
        expect_equal(out_pd3, expect_out_gg3_counts)
        expect_equal(out_gsg3, expect_out_gsg3)
    }
})


test_that("Paranoia: checking an example table", {
    s1 <- c("EZH2_msmut, C", "C", "WT", "EZH2_msmut, F_NRAS_msmut", "EZH2_msmut, D", 
"C", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "EZH2_msmut", 
"WT", "EZH2_msmut", "WT", "EZH2_msmut", "D", "WT", "EZH2_msmut", 
"WT", "WT", "F_NRAS_msmut", "WT", "WT", "WT", "C", "WT", "WT", 
"C", "WT", "WT", "WT", "WT", "WT", "C", "WT", "WT", "WT", "WT", 
"WT", "WT", "WT", "EZH2_msmut", "WT", "WT", "C, F_NRAS_msmut", 
"WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "EZH2_msmut, F_NRAS_msmut", 
"F_NRAS_msmut", "EZH2_msmut", "WT", "WT", "WT", "WT", "WT", "EZH2_msmut", 
"WT", "WT", "WT", "WT", "C, D", "WT", "WT", "EZH2_msmut", "WT", 
"WT", "C, F_NRAS_msmut", "WT", "EZH2_msmut", "WT", "WT", "WT", 
"EZH2_msmut", "C", "WT", "WT", "F_NRAS_msmut", "D", "WT", "WT", 
"WT", "WT", "WT", "EZH2_msmut", "F_NRAS_msmut", "WT", "F_NRAS_msmut", 
"EZH2_msmut, D", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", 
"WT", "WT", "EZH2_msmut, D", "C", "EZH2_msmut, D", "D", "WT", 
"WT", "WT", "WT", "WT", "WT", "EZH2_msmut", "EZH2_msmut, C", 
"WT", "EZH2_msmut", "EZH2_msmut", "C, F_NRAS_msmut", "WT", "WT", 
"WT", "WT", "WT", "WT", "WT", "WT", "D, F_NRAS_msmut", "EZH2_msmut", 
"EZH2_msmut", "WT", "WT", "EZH2_msmut", "WT", "C", "WT", "WT", 
"WT", "WT", "WT", "WT", "WT", "D", "D, F_NRAS_msmut", "EZH2_msmut, D, F_NRAS_msmut", 
"WT", "WT", "F_NRAS_msmut", "WT", "C", "WT", "WT", "EZH2_msmut, F_NRAS_msmut", 
"EZH2_msmut", "WT", "WT", "WT", "WT", "C", "WT", "WT", "WT", 
"WT", "WT", "WT", "EZH2_msmut, D", "EZH2_msmut", "C", "WT", "WT", 
"WT", "WT", "WT", "WT", "WT", "WT", "WT", "C", "F_NRAS_msmut", 
"WT", "WT", "EZH2_msmut", "WT", "C", "WT", "WT", "C", "WT", "EZH2_msmut", 
"EZH2_msmut", "EZH2_msmut", "WT", "EZH2_msmut, F_NRAS_msmut", 
"WT", "WT", "WT", "D", "WT", "WT", "WT", "WT", "WT", "WT", "C", 
"WT", "WT", "WT", "C", "WT", "WT", "C", "WT", "F_NRAS_msmut", 
"WT", "C", "WT", "WT", "WT", "WT", "WT", "C", "D", "WT", "D", 
"EZH2_msmut", "WT", "WT", "C", "WT", "WT", "WT", "WT", "WT", 
"WT", "WT", "WT", "EZH2_msmut", "D", "WT", "EZH2_msmut", "WT", 
"WT", "EZH2_msmut", "C", "C, D, F_NRAS_msmut", "D", "C, D", "D", 
"WT", "WT", "EZH2_msmut, D", "WT", "F_NRAS_msmut", "WT", "WT", 
"WT", "EZH2_msmut", "F_NRAS_msmut", "C", "WT", "WT", "C", "WT", 
"WT", "WT", "WT", "WT", "WT", "EZH2_msmut", "EZH2_msmut, C, F_NRAS_msmut", 
"D", "EZH2_msmut, F_NRAS_msmut", "WT", "WT", "WT", "C", "WT", 
"WT", "WT", "F_NRAS_msmut", "EZH2_msmut, C", "WT", "WT", "WT", 
"WT", "WT", "EZH2_msmut, D", "WT", "WT", "EZH2_msmut, D", "EZH2_msmut, F_NRAS_msmut", 
"WT", "F_NRAS_msmut", "WT", "EZH2_msmut", "EZH2_msmut", "WT", 
"WT", "WT", "WT", "EZH2_msmut, D", "EZH2_msmut", "WT", "WT", 
"WT", "WT", "WT", "WT", "WT", "WT", "EZH2_msmut", "WT", "EZH2_msmut, D", 
"EZH2_msmut, C, D", "WT", "WT", "WT", "EZH2_msmut, F_NRAS_msmut", 
"C", "C", "EZH2_msmut, C", "EZH2_msmut, C", "WT", "C", "WT", 
"EZH2_msmut", "WT", "EZH2_msmut", "WT", "F_NRAS_msmut", "EZH2_msmut", 
"WT", "WT", "WT", "D", "WT", "WT", "WT", "F_NRAS_msmut", "WT", 
"WT", "C", "WT", "WT", "EZH2_msmut", "WT", "WT", "D", "WT", "WT", 
"WT", "EZH2_msmut", "EZH2_msmut", "WT", "C", "WT", "WT", "WT", 
"WT", "WT", "WT", "WT", "C", "WT", "EZH2_msmut", "WT", "EZH2_msmut", 
"WT", "WT", "WT", "WT", "WT", "EZH2_msmut, C, F_NRAS_msmut", 
"WT", "WT", "WT", "WT", "WT", "WT", "D", "F_NRAS_msmut", "WT", 
"C", "D", "WT", "C", "C", "WT", "EZH2_msmut", "WT", "F_NRAS_msmut", 
"WT", "C", "WT", "WT", "F_NRAS_msmut", "C", "EZH2_msmut", "WT", 
"WT", "WT", "EZH2_msmut", "F_NRAS_msmut", "WT", "C, D", "WT", 
"WT", "WT", "WT", "EZH2_msmut", "EZH2_msmut", "EZH2_msmut, D", 
"WT", "WT", "WT", "EZH2_msmut", "WT", "EZH2_msmut", "WT", "EZH2_msmut", 
"D", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", 
"C", "EZH2_msmut", "C", "C", "WT", "EZH2_msmut", "WT", "F_NRAS_msmut", 
"C", "C", "WT", "WT", "WT", "EZH2_msmut", "WT", "C", "WT", "WT", 
"EZH2_msmut", "WT", "WT", "D", "F_NRAS_msmut", "WT", "WT", "WT", 
"WT", "WT", "WT", "WT", "EZH2_msmut, C", "WT", "D", "WT", "EZH2_msmut", 
"WT", "EZH2_msmut", "WT", "WT", "WT", "WT", "WT", "EZH2_msmut, F_NRAS_msmut", 
"WT", "WT", "WT", "WT", "WT", "WT", "C", "WT", "EZH2_msmut, D", 
"EZH2_msmut, D", "WT", "WT", "WT", "C", "C", "C", "WT", "EZH2_msmut, C", 
"F_NRAS_msmut", "D", "F_NRAS_msmut", "WT", "WT", "WT", "WT", 
"WT", "WT", "WT", "EZH2_msmut", "WT", "WT", "C, F_NRAS_msmut", 
"WT", "WT", "WT", "EZH2_msmut", "WT", "WT", "WT", "WT", "EZH2_msmut", 
"WT", "WT", "WT", "EZH2_msmut, F_NRAS_msmut", "WT", "WT", "EZH2_msmut", 
"WT", "EZH2_msmut, F_NRAS_msmut", "WT", "WT", "WT", "WT", "C", 
"WT", "WT", "C", "WT", "WT", "WT", "D", "WT", "WT", "EZH2_msmut, F_NRAS_msmut", 
"WT", "WT", "WT", "WT", "WT", "WT", "EZH2_msmut", "WT", "WT", 
"WT", "EZH2_msmut", "WT", "C", "WT", "EZH2_msmut, D, F_NRAS_msmut", 
"WT", "WT", "WT", "WT", "D", "WT", "WT", "WT", "WT", "WT", "EZH2_msmut", 
"WT", "WT", "WT", "WT", "WT", "D", "WT", "F_NRAS_msmut", "EZH2_msmut, D", 
"EZH2_msmut", "F_NRAS_msmut", "EZH2_msmut", "WT", "EZH2_msmut", 
"EZH2_msmut, D, F_NRAS_msmut", "C", "WT", "F_NRAS_msmut", "WT", 
"WT", "WT", "WT", "F_NRAS_msmut", "WT", "WT", "WT", "WT", "WT", 
"WT", "WT", "WT", "EZH2_msmut, C", "WT", "WT", "WT", "WT", "C", 
"WT", "WT", "WT", "WT", "EZH2_msmut, D", "F_NRAS_msmut", "WT", 
"F_NRAS_msmut", "WT", "WT", "WT", "EZH2_msmut", "WT", "WT", "WT", 
"EZH2_msmut", "EZH2_msmut", "WT", "D", "WT", "WT", "WT", "EZH2_msmut", 
"WT", "WT", "WT", "WT", "EZH2_msmut", "D", "WT", "WT", "WT", 
"EZH2_msmut", "WT", "F_NRAS_msmut", "EZH2_msmut, F_NRAS_msmut", 
"WT", "WT", "EZH2_msmut", "C", "WT", "WT", "EZH2_msmut", "WT", 
"EZH2_msmut", "WT", "WT", "C", "WT", "WT", "WT", "WT", "WT", 
"EZH2_msmut", "WT", "WT", "WT", "WT", "WT", "WT", "EZH2_msmut", 
"WT", "WT", "D", "WT", "WT", "D", "WT", "D", "WT", "C", "WT", 
"WT", "WT", "EZH2_msmut, F_NRAS_msmut", "WT", "C", "WT", "WT", 
"D", "WT", "EZH2_msmut, F_NRAS_msmut", "WT", "WT", "WT", "WT", 
"WT", "WT", "WT", "EZH2_msmut", "WT", "WT", "WT", "WT", "WT", 
"C", "WT", "WT", "EZH2_msmut, C", "F_NRAS_msmut", "EZH2_msmut, D", 
"WT", "C", "WT", "WT", "D", "WT", "WT", "WT", "WT", "EZH2_msmut, F_NRAS_msmut", 
"EZH2_msmut, F_NRAS_msmut", "EZH2_msmut, D", "EZH2_msmut", "EZH2_msmut", 
"EZH2_msmut, D", "WT", "EZH2_msmut", "C", "F_NRAS_msmut", "WT", 
"WT", "WT", "WT", "WT", "WT", "D", "EZH2_msmut, F_NRAS_msmut", 
"WT", "WT", "WT", "C", "WT", "EZH2_msmut", "WT", "WT", "C", "WT", 
"C", "EZH2_msmut, C", "D", "F_NRAS_msmut", "C", "WT", "EZH2_msmut", 
"WT", "WT", "WT", "EZH2_msmut, D, F_NRAS_msmut", "WT", "EZH2_msmut, D", 
"WT", "F_NRAS_msmut", "WT", "C", "WT", "D", "WT", "D", "WT", 
"EZH2_msmut", "WT", "WT", "C", "D", "WT", "WT", "WT", "WT", "WT", 
"WT", "C", "WT", "WT", "WT", "WT", "WT", "EZH2_msmut, C", "WT", 
"WT", "WT", "WT", "WT", "EZH2_msmut", "F_NRAS_msmut", "EZH2_msmut", 
"EZH2_msmut, D, F_NRAS_msmut", "WT", "WT", "F_NRAS_msmut", "WT", 
"WT", "WT", "WT", "WT", "EZH2_msmut, C, F_NRAS_msmut", "WT", 
"WT", "WT", "WT", "WT", "D, F_NRAS_msmut", "WT", "WT", "WT", 
"WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "EZH2_msmut, C", 
"WT", "WT", "WT", "WT", "C", "WT", "WT", "WT", "WT", "WT", "C", 
"D, F_NRAS_msmut", "WT", "WT", "WT", "WT", "EZH2_msmut, C", "WT", 
"WT", "WT", "WT", "WT", "WT", "WT", "WT", "D, F_NRAS_msmut", 
"WT", "EZH2_msmut", "WT", "WT", "WT", "C", "WT", "EZH2_msmut", 
"EZH2_msmut, D", "WT", "WT", "C", "EZH2_msmut", "WT", "EZH2_msmut, F_NRAS_msmut", 
"EZH2_msmut, D", "WT", "WT", "WT", "D", "WT", "C", "WT", "WT", 
"WT", "WT", "EZH2_msmut", "WT", "WT", "EZH2_msmut", "WT", "WT", 
"WT", "WT", "WT", "EZH2_msmut, D", "WT", "WT", "WT", "D, F_NRAS_msmut", 
"C", "WT", "WT", "WT", "WT", "EZH2_msmut, C", "WT", "WT", "C", 
"WT", "F_NRAS_msmut", "WT", "WT", "WT", "WT", "D", "WT", "WT", 
"WT", "C", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", 
"WT", "F_NRAS_msmut", "WT", "WT", "WT", "WT", "WT", "EZH2_msmut, C, D", 
"C", "WT", "WT", "WT", "WT", "WT", "WT", "D", "WT", "WT", "D", 
"WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", 
"WT", "WT", "EZH2_msmut", "EZH2_msmut", "WT", "WT", "WT", "D", 
"C", "WT", "C", "WT", "WT", "WT", "D", "WT", "D", "EZH2_msmut", 
"WT", "WT", "WT", "EZH2_msmut", "WT", "EZH2_msmut, C", "WT", 
"WT", "WT", "WT", "WT", "EZH2_msmut, F_NRAS_msmut", "F_NRAS_msmut", 
"WT", "WT", "D", "WT", "WT", "C", "WT", "WT", "WT", "WT", "WT", 
"WT", "WT", "EZH2_msmut, D", "EZH2_msmut", "EZH2_msmut, F_NRAS_msmut", 
"WT", "WT", "WT", "WT", "WT", "EZH2_msmut", "WT")

    table(s1)
    
    gene_names <- c("EZH2_msmut", "C", "D", "F_NRAS_msmut")
    n_genes <- length(gene_names)

    outgs <- generate_pD_sorted_genotypes(n_genes, gene_names)
    outgpd <- sample_to_pD_order(s1, n_genes, gene_names)

    ## I can't do  table(s1)[outgs]
    ## as the order of the names of genes is different in some genotypes
    ## and that is precisely what I want to verify: it works OK if I resort them

    s2 <- s1
    s2 <- gsub("EZH2_msmut, C, F_NRAS_msmut", "C, EZH2_msmut, F_NRAS_msmut",
               s2, fixed = TRUE)
    s2 <- gsub("EZH2_msmut, D, F_NRAS_msmut", "D, EZH2_msmut, F_NRAS_msmut",
               s2, fixed = TRUE)
    ## s2 <- gsub("C, EZH2_msmut, D", "C, D, EZH2_msmut",
    ##            s2, fixed = TRUE)
    s2 <- gsub("EZH2_msmut, C, D", "C, D, EZH2_msmut",
               s2, fixed = TRUE)
    s2 <- gsub("EZH2_msmut, D", "D, EZH2_msmut", s2, fixed = TRUE)
    s2 <- gsub("EZH2_msmut, C", "C, EZH2_msmut", s2, fixed = TRUE)

    nums <- table(s2)[outgs]
    nums[is.na(nums)] <- 0
    
    expect_equal(as.vector(nums), outgpd, ignore_attr = TRUE)
})




test_that("Test with same order of genes", {
    ## This could fail when things are OK
    ## because in OncoSimulR I use mixedsort
    ## which sorts slightly differently and ignores case
    N <- 10000
    iters <- 20
    for(i in 1:iters) {
        ngenes <- 7
        gn <- vector(mode = "character", length = ngenes)

        ## create weird gene names
        for(g in seq_len(ngenes)) {
            l1 <- sample(c(LETTERS, letters), 1)
            rest <-  c(sample(
                c(
                    sample(c("_",  "=", "?", "#", "@", "%", "&", "!"), 4, replace = TRUE)
                    , sample(0:9, 4, replace = TRUE)
                    , sample(c(letters, LETTERS), 4, replace = TRUE)
                )))
            
            rest <- paste(rest, sep = "", collapse = "")
            gn[g] <- paste0(l1, rest)
        }
        
        gn <- evam_string_sort(gn)
        rf <- OncoSimulR::rfitness(ngenes)
        colnames(rf)[1:ngenes] <- gn
        ag <- OncoSimulR::evalAllGenotypes(OncoSimulR::allFitnessEffects(genotFitness = rf),
                               addwt = TRUE)

        ## Undo what gtools::mixedsort has done or this test
        ## can fail because the table has, as entries, the original
        ## genotypes, those with the genes sorted by mixedsort.
        genotypes_resorted <-
            unname(vapply(ag[, "Genotype"],
                   function(u)
                       paste(evam_string_sort(strsplit(u, split = ", ", fixed = TRUE)[[1]]),
                             collapse = ", "), "something"))
        
        sx <- sample(genotypes_resorted, N, prob = ag[, "Fitness"], replace = TRUE)

        outgs <- generate_pD_sorted_genotypes(ngenes, gn)
        outgpd <-  sample_to_pD_order(sx, ngenes, gn)

        ot <- table(sx)
        oto <- ot[outgs]
        oto[is.na(oto)] <- 0
        expect_equal(as.vector(oto), outgpd, ignore_attr = TRUE)
    }
})


cat("\n Done test.generate-sorted-sample-to-pD. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
