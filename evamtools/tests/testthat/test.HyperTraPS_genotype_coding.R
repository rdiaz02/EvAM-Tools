## This isn't really a test, but code to show the coding and decoding
## of genotypes is correct.

## The binary state shown as 0,1,0,1, for example,
## when labels are A, B, C, D, is B, D.

## We use this conversion in probs_from_HyperTraPS_continuous, using binary2str_label,
## and a different conversion (from the integer coding of HyperTraPS
## to the states) is used in run_HyperTraPS, when obtaining the
## transition matrix.

## The code below shows:
## a) The names in the transition matrix are correct.
##    (shows code to obtain trans mat is right)
## b) The genotypes along correspond to the decoding.
##    (shows conversion from binary to genotype name
##     in probs_from_HyperTraPS is correct).


if (FALSE) {
    d1 <- matrix(
        c(
            rep(c(0, 1, 0), 300) #B
          , rep(c(0, 1, 1), 200) #B,C
          , rep(c(1, 1, 1), 300) #A,B,C
          , rep(c(0, 0, 0), 300) # WT
        ), ncol = 3, byrow = TRUE
    )
    colnames(d1) <- LETTERS[1:3]

    o1 <- evam(d1,
               methods = c("HyperTraPS"),
               hyper_traps_opts = list(length = 3,
                                       model = 2,
                                       walkers = 400,
                                       seed = 1))
    ## See how main transitions are
    ## WT -> B -> BC -> ABC
    o1$HyperTraPS_trans_mat
    o1$HyperTraPS_trans_mat["WT", ]
    o1$HyperTraPS_trans_mat["B", ]
    o1$HyperTraPS_trans_mat["B, C", ]



    d2 <- matrix(
        c(
            rep(c(0, 0, 1, 0), 300) #C
          , rep(c(0, 0, 1, 1), 200) #CD
          , rep(c(1, 0, 1, 1), 100) #ACD
          , rep(c(1, 1, 1, 1), 50) # ABCD
          , rep(c(0, 0, 0, 0), 10) # WT
        ), ncol = 4, byrow = TRUE
    )
    colnames(d2) <- LETTERS[1:4]

    o2 <- evam(d2,
               methods = c("HyperTraPS"),
               hyper_traps_opts = list(length = 3,
                                       model = 2,
                                       walkers = 400,
                                       seed = 23))
    ## See how main transitions are
    ## WT -> C -> CD -> ACD
    o2$HyperTraPS_trans_mat
    o2$HyperTraPS_trans_mat["C", ]
    o2$HyperTraPS_trans_mat["C, D", ]
    o2$HyperTraPS_trans_mat["A, C, D", ]


    o3 <- evam(d2,
               methods = c("HyperTraPS"),
               hyper_traps_opts = list(length = 3,
                                       model = 2,
                                       walkers = 200,
                                       seed = 23))

    o3$HyperTraPS_trans_mat["WT", ]
    o3$HyperTraPS_trans_mat["C", ]
    o3$HyperTraPS_trans_mat["C, D", ]
    o3$HyperTraPS_trans_mat["A, C, D", ]
}
