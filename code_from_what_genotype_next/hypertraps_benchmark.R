source("hypertraps-process.R")

N <- 100
na <- N
nc <- N + round( 10 * runif(1))
nab <- 1.6 * N + round( 10 * runif(1))
ncd <- 1.5 * N + round( 10 * runif(1))
n00 <- round( 10 * runif(1))
dB <- matrix(
  c(
    rep(c(1, 0, 0, 0), na) 
    , rep(c(0, 0, 1, 0), nc)
    , rep(c(1, 1, 0, 0), nab)
    , rep(c(0, 0, 1, 1), ncd)        
    , rep(c(0, 0, 0, 0), n00)
  ), ncol = 4, byrow = TRUE
)
colnames(dB) <- LETTERS[1:4]

iterations <- c(5, 50, 500, 1000, 2500)
bi <- c(500, 5000, 20000, 45000)
rs <- c(10, 100, 200, 500)
out <- matrix(ncol=4, nrow=length(iterations)*length(bi)*length(rs))
counter <- 0 
for (it in iterations){
  for (b in bi){
    for (r in rs){
      tmp_folder = paste(c("runs", it, "bi", b, "r", r), collapse="_")
      time_sampling_posterior <- do_HyperTraPS(dB, tmp_folder=tmp_folder, runs=it, bi=b, r=r)
      counter <- counter + 1
      out[counter, ] <- c(it, b, r, time_sampling_posterior)
      data_out <- data.frame(out)
      colnames(data_out) <- c("Iteration", "bi", "r", "Time elapsed")
      write.csv(data_out, "/home/pablo/CPM-SSWM-Sampling/code_from_what_genotype_next/elapsed_time_hypertraps.csv")
    }
  }
}

