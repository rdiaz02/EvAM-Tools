## Issues to document

## a) different score with the same theta: makes sense if this is the
## penalized likelihood

## b) CV choosing ever greater values of lambda: but I think when using CV
## you want the largest likelihood, not largest penalized likelihood


luad <- read.csv("~/Downloads/LUAD_n12.csv")

luad100 <- luad[1:100, ]

opts0_MHN_python <- list(
  omp_threads = 1, ##
  Type = "cMHN", ## "oMHN" won't be used for now
  Penalty = "L1", ## SYM_SPARSE, L1 or L2. We'll use SYM and L1
  seed = 1,
  steps = 10,
  nfolds = 5,
  maxit = 5000,
  reltol = 1e-07,
  round_result = "True", ## Rounding was used in the 2020 paper
  ## and it is the default in the Python examples
  show_progressbar = "True",
  return_lambda_scores = "True",
  pick_1se = "False" ## For real, probably wants True
)

setup_python_MHN()

r10 <- run_MHN_python(luad100,
                      opts = c(opts0_MHN_python,
                                lambda_min = 10,
                                lambda_max = 10))
r10$out$lambda_used
r10$out$score
r10$out$theta


r1 <- run_MHN_python(luad100,
                     opts = c(opts0_MHN_python,
                               lambda_min = 1,
                              lambda_max = 1))
r1$out$lambda_used
r1$out$score
r1$out$thetapy_env$lambda_scores


r08 <- run_MHN_python(luad100,
                      opts = c(opts0_MHN_python,
                               lambda_min = .08,
                               lambda_max = .08))
r08$out$lambda_used
r08$out$score
r08$out$theta



## Look at CV scores


rx <- run_MHN_python(luad100,
                     opts = c(opts0_MHN_python,
                              lambda_min = 0.01,
                              lambda_max = 20))
rx$out$lambda_used
rx$out$score
rx$out$theta








luad <- read.csv("~/Downloads/LUAD_n12.csv")

luad100 <- luad[1:100, ]

luad200 <- luad[1:200, ]


setup_python_MHN()

t1 <- do_MHN_python(luad100,
                    opts = c(opts0_MHN_python,
                             lambda_min = 0.02,
                             lambda_max = 27))



t2 <- do_MHN_python(luad100,
                    opts = c(opts0_MHN_python,
                             lambda_min = 0.01,
                             lambda_max = 5))


t3 <- do_MHN_python(luad100,
                    opts = c(opts0_MHN_python,
                             lambda_min = 0.01,
                             lambda_max = 8))
t3$result$meta
t3$lambda_scores


t4 <- do_MHN_python(luad200,
                    opts = c(opts0_MHN_python,
                             lambda_min = 0.02,
                             lambda_max = 1))
t4$result$meta
t4$lambda_scores

## This is a good example
t4 <- do_MHN_python(luad200,
                    opts = c(opts0_MHN_python,
                             lambda_min = 0.02,
                             lambda_max = 1))
t4$result$meta
t4$lambda_scores

t4b <- do_MHN_python(luad200,
                     opts = c(opts0_MHN_python,
                              lambda_min = 0.01,
                              lambda_max = 1))
t4b$result$meta
t4b$lambda_scores


t5 <- do_MHN_python(luad200,
                    opts = c(opts0_MHN_python,
                             lambda_min = 0.02,
                             lambda_max = 0.02))
t5$result$meta
t5$lambda_scores
