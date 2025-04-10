## Example data generation
rc <- random_evam(model = "CBN", ngenes = 8)
rd <- sample_evam(rc, N = 20000, obs_noise = 0.05)
ddd <- rd[[2]]

rc2 <- random_evam(model = "MHN", ngenes = 8)
rd2 <- sample_evam(rc2, N = 20000, obs_noise = 0.05)
ddd2 <- rd2[[2]]



## Session-wide
setup_python_MHN()


opts0_MHN_python <- list(
  omp_threads = 1, ##
  Type = "cMHN", ## "oMHN" won't be used for now
  Penalty = "L1", ## SYM_SPARSE, L1 or L2. We'll use SYM and L1
  seed = NA,
  steps = 10,
  nfolds = 5,
  maxit = 5000,
  reltol = 1e-07,
  round_result = "True", ## Rounding was used in the 2020 paper
  ## and it is the default in the Python examples
  show_progressbar = "True",
  verbose = TRUE
)


## Interesting it gives a warning if lambda_max is 1
## but not with other values at the boundary
r1 <- run_MHN_python(ddd[1:1000, ],
                     opts = c(opts0_MHN_python,
                              lambda_min = 0.1/100,
                              lambda_max = 0.5))
r1$out$lambda_used



r11 <- run_MHN_python(ddd2[1:100, ],
                      opts = c(opts0_MHN_python,
                              lambda_min = 0.1/100,
                              lambda_max = 0.1))
r11$out$lambda_used
r11$out$score

r11 <- run_MHN_python(ddd2[1:100, ],
                      opts = c(opts0_MHN_python,
                               lambda_min = 0.1/100,
                               lambda_max = 0.05))
r11$out$lambda_used
r11$out$score
r11$out$theta

r11 <- run_MHN_python(ddd2[1:100, ],
                      opts = c(opts0_MHN_python,
                               lambda_min = 0.1/100,
                               lambda_max = 1))
r11$out$lambda_used
r11$out$score
r11$out$theta

r11 <- run_MHN_python(ddd2[1:100, ],
                      opts = c(opts0_MHN_python,
                               lambda_min = 0.1/100,
                               lambda_max = 10))
r11$out$lambda_used
r11$out$score
r11$out$theta

r111 <- run_MHN_python(ddd2[1:100, ],
                       opts = c(opts0_MHN_python,
                               lambda_min = 1,
                               lambda_max = 1))
r111$out$lambda_used
r111$out$score
r111$out$theta

r111 <- run_MHN_python(ddd2[1:100, ],
                       opts = c(opts0_MHN_python,
                                lambda_min = 10,
                                lambda_max = 10))
r111$out$lambda_used
r111$out$score
r111$out$theta



## This can be slow! 1100 seconds
system.time(r2 <- run_MHN_python(ddd, opts = c(opts0_MHN_python,
                                               lambda_min = 0.1/nrow(ddd),
                                               lambda_max = 100/nrow(ddd))))




### Running the LUAD data
##  As in the example here https://github.com/spang-lab/LearnMHN/blob/main/demo/demo.ipynb
## But note I am using cMHN, not oMHN here

## Download and do
luad <- read.csv("~/Downloads/LUAD_n12.csv")
## Or if you already have it, load it from an RData

## 350 seconds
system.time(
  rluad <- run_MHN_python(luad, opts = c(opts0_MHN_python,
                                         lambda_min = 0.1/nrow(luad),
                                         lambda_max = 100/nrow(luad))))



### Minimal checks

#### evamtools and Python calls

## Compare the Python and the evamtools implementations
## with a fixed lambda. Results are not exactly identical
## but close enough.

rc <- random_evam(model = "CBN", ngenes = 8)
rd <- sample_evam(rc, N = 10000, obs_noise = 0.05)
x <- rd[[2]]

opts_fixed_l_MHN_python <- list(
  omp_threads = 1, ##
  Type = "cMHN", ## "oMHN" won't be used for now
  Penalty = "L1",
  ## the rest are irrelevant if no CV, except maxit and further
  seed = NA,
  steps = 4, # 10,
  nfolds = 3, # 5,
  ## Default is 5000
  maxit = 5000, ## 5000
  reltol = 1e-07,
  round_result = "True",
  show_progressbar = "False",
  verbose = TRUE
)

out_mhn_evam <- evam(x, methods = "MHN",
                     mhn_opts = list(lambda = 1/nrow(x),
                                     omp_threads = 1))


## Much slower for large number of samples
out_mhn_python <- run_MHN_python(x,
                                 opts = c(opts_fixed_l_MHN_python,
                                          lambda_min = 1/nrow(x),
                                          lambda_max = 1/nrow(x)))
## Check lambda used

stopifnot(identical(out_mhn_python$out$lambda_used,
                    1/nrow(x)))

## Are we maxing the iterations?
out_mhn_python$out$nit
out_mhn_python$out$message

out_mhn_python$out$predicted_genotype_freqs
out_mhn_evam$MHN_predicted_genotype_freqs

out_mhn_python$out$theta
out_mhn_evam$MHN_theta




#### Using the correct theta matrix

## Actually, a moot point as we know we are doing things
## correctly based on the comparison of evam and Python codes

## A simple check that we are using the correct theta matrix
## Compare genotype frequencies of a simulated sample using the
## Python with the predicted frequencies we get after
## reading the theta matrix, obtaining the  transition rate matrix
## and from it the genotype frequencies


f_compare_pred_sampled <- function(N = 5e6, k = 6) {
  rc <- random_evam(model = "CBN", ngenes = k)
  rd <- sample_evam(rc, N = 1000, obs_noise = 0.05)
  x <- rd[[2]]

  t1 <- do_MHN_python(x, opts = c(opts0_MHN_python,
                                  lambda_min = 0.1/nrow(x),
                                  lambda_max = 100/nrow(x)))
  t1_data <- t1$sample_artificial_data(N, as_dataframe = "True")
  t1_freq <- OncoSimulR::sampledGenotypes(t1_data)
  t1_freqv <- t1_freq$Freq
  names(t1_freqv) <- t1_freq$Genotype
  t1_freqv <- evamtools:::reorder_to_standard_order(t1_freqv)
  t1_freqvp <- t1_freqv/sum(t1_freqv)

  ## Exactly as in run_MHN_python
  t1_theta <- t1$log_theta
  colnames(t1_theta) <- rownames(t1_theta) <- colnames(x)
  t1_trm <- evamtools:::theta_to_trans_rate_3(t1_theta)
  t1_predgenots <- evamtools:::probs_from_trm(t1_trm)

  if (!identical(names(t1_predgenots), names(t1_freqv))) {
    message("Genotype names not identical")
    browser()
  }

  m_freqs <- cbind(Simul = t1_freqvp, Pred = t1_predgenots)
  rownames(m_freqs) <- names(t1_freqvp)

  print(m_freqs)

  obs_pred <- round(N * t1_predgenots)

  print(chisq.test(obs_pred, t1_freqv))
}
