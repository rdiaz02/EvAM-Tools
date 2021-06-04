pwd0 <- getwd()
setwd("../code_from_what_genotype_next/")
source("code-all-methods-minimal.R")
setwd(pwd0)
rm(pwd0)

dB_AND <- matrix(
  c(
    rep(c(1, 0, 0, 0), 100)
    , rep(c(1, 1, 0, 0), 100)
    , rep(c(1, 1, 1, 0), 100)
    , rep(c(1, 1, 1, 1), 100)
    , rep(c(0, 0, 0, 0), 10)
  ), ncol = 4, byrow = TRUE
)
colnames(dB_AND) <- LETTERS[1:4]

out <- all_methods_2_trans_mat(dB_AND, HT_folder = "HyperTraPS_examples/HP_AND")
# do_HyperTraPS(dB_AND, "HyperTraPS_examples/HP_AND", runs = 500, bi=200, dry_run = TRUE)