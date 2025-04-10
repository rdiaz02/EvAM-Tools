## Example data generation
rc <- random_evam(model = "CBN", ngenes = 5)
rd <- sample_evam(rc, N = 1000, obs_noise = 0.05)
ddd <- rd[[2]]


## Session-wide
setup_python_MHN()



opts <- list()
opts$MHN_python <- list(
  Type = "cMHN", ## "oMHN" won't be used for now
  Penalty = "SYM_SPARSE", ## or L1 or L2. We'll use SYM and L1
  seed = NA,
  ## FIXME: use the usual apparatus
  lambda_min = 0.001, ## 0.1/nrow(x),
  lambda_max = 0.1, ## 100/nrow(x),
  steps = 4, # 10,
  nfolds = 3, # 5,
  ## Default is 5000
  maxit = 100, ## 5000
  show_progressbar = "True"
)

## Each run
do_MHN_python(ddd, opts)
