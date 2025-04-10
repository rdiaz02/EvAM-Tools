## Use the Python implementation of MHN

## I try to use Javier Lopez de Lema's code structure
## But since this might not yet be fully integrated into evam
## I may add a few things here that will later be moved out

## FIXME: pass the path to python
## e.g., ~/.local/python3.12-venv/bin/


library(evamtools)
library("reticulate")



## This is somewhat time consuming, so do it once per session?
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


run_MHN_python <- function(x, opts) {
  RhpcBLASctl::omp_set_num_threads(opts$omp_threads)
  time_out <- system.time({
    out <- do_MHN_python(x, FIXME_lambda = opts$lambda)
    out <- c(out,
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



do_MHN_python <- function(x,
                          opts
                          ## dirname = NULL,
                          ## filename = "mhn_input.csv",
                          ## silent = TRUE,
                          ## rmfile = TRUE
                          ) {
  set_MHN_python_options(opts$MHN_python)
  data_2_MHN_load_data_matrix(x)
  set_run_cv_MHN_python(opts$MHN_python)
  train_MHN_python(opts$MHN_python)
  MHN_result <- py$opt$result

}
