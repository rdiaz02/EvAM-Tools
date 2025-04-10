## Example data generation
rc <- random_evam(model = "CBN", ngenes = 5)
rd <- sample_evam(rc, N = 1000, obs_noise = 0.05)
ddd <- rd[[2]]


## Session-wide
setup_python_MHN()



opts0_MHN_python <- list(
  omp_threads = 1, ##
  Type = "cMHN", ## "oMHN" won't be used for now
  Penalty = "SYM_SPARSE", ## or L1 or L2. We'll use SYM and L1
  seed = NA,
  steps = 4, # 10,
  nfolds = 3, # 5,
  ## Default is 5000
  maxit = 100, ## 5000
  show_progressbar = "False",
  verbose = TRUE
)

## ## Each run
## do_MHN_python(ddd, opts)




r1 <- run_MHN_python(ddd, opts = c(opts0_MHN_python,
                                   lambda_min = 0.1/nrow(ddd),
                                   lambda_max = 100/nrow(ddd)))

## FIXME
## return the lambda used
## return the maxit, message, nit, etc
## compare the predicted genotypes
