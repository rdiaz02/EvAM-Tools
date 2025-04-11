### What is this?

## Examples of uses of the MHN_python code, including minimal tests
## that things are correct. These tests include

## - that results agree with those from the evamtools code for fixed
##   lambda

## - that the theta matrix is ordered the way we think
##   (this ain't explicit in their document, as far as I can tell)

## - that we can run the code parallelized: my initial versions of this
##   code had weird behavior because some objects remamind in the
##   environment.



### Set the python environment.

setup_python_MHN()


### Example calls

## Example data generation
rc <- random_evam(model = "CBN", ngenes = 8)
rd <- sample_evam(rc, N = 20000, obs_noise = 0.05)
ddd <- rd[[2]]

## This is a sensible set of defaults for default runs
opts_MHN_python <- list(
  omp_threads = 1,
  Type = "cMHN", ## "oMHN" won't be used for now
  Penalty = "L1", ## SYM_SPARSE, L1 or L2. We'll use SYM and L1
  seed = NA,
  steps = 10,
  nfolds = 5,
  maxit = 5000,
  reltol = 1e-07,
  round_result = "True", ## Rounding was used in the 2020 paper
  ## and it is the default in the Python examples
  show_progressbar = "False",
  return_lambda_scores = "False",
  pick_1se = "True" ## For real, probably wants True
)


## Can get a warning here because best lambda at the boundary
r1 <- run_MHN_python(ddd[1:100, ],
                     opts = c(opts_MHN_python,
                              lambda_min = 0.1/100,
                              lambda_max = 0.5))
r1$out$lambda_used





### Minimal checks

#### evamtools and Python calls

## Compare the Python and the evamtools implementations
## with a fixed lambda. Results are not exactly identical
## but close enough.

rc <- random_evam(model = "CBN", ngenes = 8)
rd <- sample_evam(rc, N = 10000, obs_noise = 0.05)
x <- rd[[2]]



out_mhn_evam <- evam(x, methods = "MHN",
                     mhn_opts = list(lambda = 1/nrow(x),
                                     omp_threads = 1))


## Much slower for large number of samples
out_mhn_python <- run_MHN_python(x,
                                 opts = c(opts_MHN_python,
                                          lambda_min = 1/nrow(x),
                                          lambda_max = 1/nrow(x)))

## Check lambda used
stopifnot(identical(out_mhn_python$out$lambda_used,
                    1/nrow(x)))

## Are we maxing the iterations?
out_mhn_python$out$nit
out_mhn_python$out$message

cbind(out_mhn_python$out$predicted_genotype_freqs,
      out_mhn_evam$MHN_predicted_genotype_freqs)

out_mhn_python$out$theta
out_mhn_evam$MHN_theta




#### Using the correct theta matrix

## A simple check that we are using the correct theta matrix
## Compare genotype frequencies of a simulated sample using the
## Python with the predicted frequencies we get after
## reading the theta matrix, obtaining the  transition rate matrix
## and from it the genotype frequencies

## Actually, a moot point as we know we are doing things
## correctly based on the comparison of evam and Python code above.


f_compare_pred_sampled <- function(N = 5e6, k = 6) {
  rc <- random_evam(model = "CBN", ngenes = k)
  rd <- sample_evam(rc, N = 100, obs_noise = 0.05)
  x <- rd[[2]]

  t1 <- do_MHN_python(x, opts = c(opts0_MHN_python,
                                  lambda_min = 0.1/nrow(x),
                                  lambda_max = 100/nrow(x)))
  t1_data <- t1$result$sample_artificial_data(N, as_dataframe = "True")
  t1_freq <- OncoSimulR::sampledGenotypes(t1_data)
  t1_freqv <- t1_freq$Freq
  names(t1_freqv) <- t1_freq$Genotype
  t1_freqv <- evamtools:::reorder_to_standard_order(t1_freqv)
  t1_freqvp <- t1_freqv/sum(t1_freqv)

  ## Exactly as in run_MHN_python
  t1_theta <- t1$result$log_theta
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

f_compare_pred_sampled()


### Examples parallelized call


library(parallel)

## Common options used. Some are too verbose for routine use
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


## We will use the LUAD data, from their repo.
## https://github.com/spang-lab/LearnMHN/blob/main/demo/demo.ipynb
## We do this in other places too, skip if you have the data.
luad <- read.csv("~/Downloads/LUAD_n12.csv")
luad100 <- luad[1:100, ]

#### Verify parallel runs can use different data
dd <- list(luad100[, 1:4],
           luad100[, 6:12],
           luad100[1:50, c(1, 3, 6, 7)],
           luad100[1:50, ],
           luad100[1:50, c(2, 3, 4, 9, 8)]
           )

pa1 <- mclapply(dd, function(x)
  run_MHN_python(x, opts = c(opts0_MHN_python,
                             lambda_min = 0.01,
                             lambda_max = 0.08)),
  mc.cores = 2 ## Force core reuse
  )

## Note the dimensions
pa1[[1]]$out$theta
pa1[[2]]$out$theta
pa1[[3]]$out$theta

## Minimal test
lapply(seq_along(dd),
       function(i)
         stopifnot(identical(colnames(pa1[[i]]$out$theta),
                             colnames(dd[[i]]))))




#### Verify parallel runs can use different options

opt1 <- c(opts0_MHN_python,  lambda_min = 0.01,  lambda_max = 0.02)
opt2 <- c(opts0_MHN_python,  lambda_min = 0.015,  lambda_max = 0.015)
opt3 <- c(opts0_MHN_python,  lambda_min = 0.02,  lambda_max = 0.18)
opt4 <- c(opts0_MHN_python,  lambda_min = 0.31,  lambda_max = 0.31)
opt5 <- c(opts0_MHN_python,  lambda_min = 0.06,  lambda_max = 0.08)
opt6 <- c(opts0_MHN_python,  lambda_min = 0.04,  lambda_max = 0.68)

opts <- list(opt1, opt2, opt3, opt4, opt5, opt6)

## We also check different data
dd2 <- c(dd, list(luad100))

pa2 <- mclapply(1:6, function(i)
  run_MHN_python(dd2[[i]], opts = opts[[i]]),
  mc.cores = 2 ## Force core reuse
  )


## Check data used correct
lapply(seq_along(dd2),
       function(i)
         stopifnot(identical(colnames(pa2[[i]]$out$theta),
                             colnames(dd2[[i]]))))

## This we know should be true
## For runs 2 and 4, no CV, and thus the used lambda
stopifnot(pa2[[2]]$out$lambda_used == 0.015)
stopifnot(pa2[[4]]$out$lambda_used == 0.31)

## For the rest, we know the range of the lambda scores
## Yes, we need the all equal.
stopifnot(all.equal(pa2[[1]]$out$lambda_scores[1, 1], 0.01))
stopifnot(all.equal(pa2[[3]]$out$lambda_scores[1, 1], 0.02))
stopifnot(all.equal(pa2[[5]]$out$lambda_scores[1, 1], 0.06))
stopifnot(all.equal(pa2[[6]]$out$lambda_scores[1, 1], 0.04))

stopifnot(all.equal(pa2[[1]]$out$lambda_scores[10, 1], 0.02))
stopifnot(all.equal(pa2[[3]]$out$lambda_scores[10, 1], 0.18))
stopifnot(all.equal(pa2[[5]]$out$lambda_scores[10, 1], 0.08))
stopifnot(all.equal(pa2[[6]]$out$lambda_scores[10, 1], 0.68))


### Miscell other example calls: possible non-intuitive behavior of CV choice

## You can play with this, changing the settings of
## lambdas and pick_1se to see that it behaves as it should
## but that with pick_1se = True one can get undesired
## behavior with huge lambda_max

## I use do_MHN, not run_MHN to get rid off the overhead
## of the post-run transformations.

t4b <- do_MHN_python(luad100,
                     opts = c(opts0_MHN_python,
                              lambda_min = 0.01,
                              lambda_max = 1))
t4b$result$meta
t4b$lambda_scores


opts0_MHN_python2 <- list(
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
  pick_1se = "True" ## For real, probably wants True
)


t5b <- do_MHN_python(luad100,
                     opts = c(opts0_MHN_python2,
                              lambda_min = 0.02,
                              lambda_max = 39))
t5b$result$meta
t5b$lambda_scores



### Running the LUAD data: this is slow!
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
