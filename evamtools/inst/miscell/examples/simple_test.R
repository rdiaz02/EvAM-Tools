library(devtools)
load_all()

dB_c1 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), 300) #A
    , rep(c(0, 0, 1, 0, 0), 300) #C
    , rep(c(1, 1, 0, 0, 0), 200) #AB
    , rep(c(0, 0, 1, 1, 0), 200) #CD
    , rep(c(1, 1, 1, 0, 0), 100) #ABC
    , rep(c(1, 0, 1, 1, 0), 100) #ACD
    , rep(c(1, 1, 0, 0, 1), 100) #ABE
    , rep(c(0, 0, 1, 1, 1), 100) #CDE
    , rep(c(1, 1, 1, 0, 1), 100) #ABCE
    , rep(c(1, 0, 1, 1, 1), 100) #ACDE
    , rep(c(1, 1, 1, 1, 0), 50) # ABCD
    , rep(c(0, 0, 0, 0, 0), 10) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c1) <- LETTERS[1:5]


# out_with_simulations <- readRDS("cpm_out_with_simulations.rds")
out <- all_methods_2_trans_mat(dB_c1)
# do_HyperTraPS(dB_AND, "HyperTraPS_examples/HP_AND", runs = 500, bi=200, dry_run = TRUE)
# plot_DAG_fg(out, dB_c1, plot_type = "matrix")

out_with_simulations <- sample_all_CPMs(out, 100000, 5)
saveRDS(out_with_simulations, "cpm_out_with_simulations.rds")
# png("graph.png", width = 1000, height = 600, units = "px")
# plot_DAG_fg(out_with_simulations, dB_c1, plot_type = "genotypes")
# dev.off()
# png("transitions.png", width = 1000, height = 600, units = "px")
# plot_DAG_fg(out_with_simulations, dB_c1, plot_type = "transitions", top_paths = 4)
# dev.off()