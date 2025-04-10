## Use the Python implementation of MHN

## I try to use Javier Lopez de Lema's code structure
## But since this might not yet be fully integrated into evam
## I may add a few things here that will later be moved out


library(evamtools)
library("reticulate")

## This is somewhat time consuming, so do it once per session
## Pass the venv you use, most likely 3.12
setup_python_MHN <- function(python_version = "~/.local/python3.12-venv/bin/python3.12") {
  use_python(python = python_version)
  mhn <- import("mhn")
  mhn_version <- py_get_attr(mhn, "__version__")
  message("Using MHN version ", mhn_version,
          " with Python ", python_version)
  py_run_string("from mhn.optimizers import Optimizer")
  ## For data handling
  py_run_string("
import numpy as np
import pandas as pd
")
}

## Might change during session
set_MHN_python_options <- function(opts) {
  ## FIXME: we need error checking for options
  ## Type: one of cMHN or oMHN
  py_run_string(paste0("opt = Optimizer(Optimizer.MHNType.", opts$Type, ")"))
  ## Penalty: one of L1, L2, SYM_SPARSE
  py_run_string(paste0("opt.set_penalty(opt.Penalty.", opts$Penalty, ")"))
  ## Leave init_theta at default. So do not touch it
}

## After reading the data
set_run_cv_MHN_python <- function(opts) {
  if(!is.null(opts$seed) && !is.na(opts$seed)) {
    py_run_string(paste0("mhn.set_seed(", opts$seed, ")"))
  }
  py_run_string(paste0("cv_lambda = opt.lambda_from_cv(",
                       "lambda_min=", opts$lambda_min, ",",
                       "lambda_max=", opts$lambda_max, ",",
                       "steps=", opts$steps, ",",
                       "nfolds=", opts$nfolds, ",",
                       "show_progressbar=", opts$show_progressbar,
                       ")"))
}

## matrix, list -> output
##  opts is the opts$MHN_python list
run_MHN_python <- function(x, opts) {
  check_MHN_python_options(opts)
  RhpcBLASctl::omp_set_num_threads(opts$omp_threads)
  time_out <- system.time({
    outp <- do_MHN_python(x, opts)
    out <- list()
    theta <- outp$log_theta
    colnames(theta) <- rownames(theta) <- colnames(x)
    out$theta <- theta
    rm(theta)
    out$transitionRateMatrix <- theta_to_trans_rate_3(out$theta)
    out <- c(out,
             lambda_used = outp$meta$lambda,
             maxit       = outp$meta$maxit,
             reltol       = outp$meta$reltol,
             score = outp$meta$score,
             message = outp$meta$message,
             status = outp$meta$status,
             nit = outp$meta$nit,
             predicted_genotype_freqs =
               list(probs_from_trm(out$transitionRateMatrix))
             )
  })["elapsed"]

  return(list(time_out = time_out, out = out))
}


train_MHN_python <- function(opts) {
  py_run_string(paste0("opt.train(",
                       "lam=cv_lambda,",
                       "maxit=", opts$maxit, ",",
                       "round_result=False,",
                       ")"
                       ))

}


data_2_MHN_load_data_matrix <- function(x) {
  py <- reticulate::py
  py$x <- x
  py$colnames_x <- colnames(x)
  py_run_string("
data_matrix = pd.DataFrame(np.array(x, dtype=np.int32), columns=colnames_x)
opt.load_data_matrix(data_matrix)
  ")
}


## matrix, list (this is the MHN_python list component)
do_MHN_python <- function(x,
                          opts) {
  if (opts$verbose) message("Starting.")
  if (opts$verbose) message("     setting options")
  set_MHN_python_options(opts)
  if (opts$verbose) message("     loading data")
  data_2_MHN_load_data_matrix(x)
  if (opts$verbose) message("     starting run_cv")
  set_run_cv_MHN_python(opts)
  if (opts$verbose) message("     starting training")
  train_MHN_python(opts)
  if (opts$verbose) message("Finished run.")
  ## MHN_result <- py$opt$result
  return(py$opt$result)
}


check_MHN_python_options <- function(opts) {
  stopifnot(is.numeric(opts$omp_threads))
  stopifnot(opts$Type %in% c("oMHN", "cMHN"))
  stopifnot(opts$Penalty %in% c("L1", "L2", "SYM_SPARSE"))
  stopifnot(is.numeric(opts$lambda_min))
  stopifnot(is.numeric(opts$lambda_max))
  stopifnot(is.numeric(opts$steps))
  stopifnot(is.numeric(opts$nfolds))
  stopifnot(is.numeric(opts$maxit))
  stopifnot(is.logical(opts$verbose))
}


## The next two make code possibly harder to understand
## to set lambda_min and lambda_max from the object itself
eval_MHN_python_opts <- function(opts, x){
  return(c(opts, lambda_min = 0.1/nrow(x), lambda_max = 100/nrow(x)))
}

## wrapper of a single argument
wrap_run_MHN_python <- function(x, opts) {
  run_MHN_python(ddd, opts = eval_MHN_python_opts(opts, x))
}


library(codetools)
checkUsageEnv(env = .GlobalEnv)
