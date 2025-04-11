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
}

## matrix, list -> output
##  opts is the opts$MHN_python list

## A wrapper to do_MHN_python with the standard additional
## manipulations (transition rate matrix and
## predicted genotype frequencies, mainly)
run_MHN_python <- function(x, opts) {
  check_MHN_python_options(opts)
  RhpcBLASctl::omp_set_num_threads(opts$omp_threads)
  time_out <- system.time({
    time_py_call <- system.time({
      outp <- do_MHN_python(x, opts)}, gcFirst = FALSE)["elapsed"]
    out <- list()
    theta <- outp$result$log_theta
    colnames(theta) <- rownames(theta) <- colnames(x)
    out$theta <- theta
    rm(theta)
    out$transitionRateMatrix <- theta_to_trans_rate_3(out$theta)

    out <- c(out,
             lambda_used = outp$result$meta$lambda,
             maxit       = outp$result$meta$maxit,
             reltol      = outp$result$meta$reltol,
             score       = outp$result$meta$score,
             message     = outp$result$meta$message,
             status      = outp$result$meta$status,
             nit         = outp$result$meta$nit,
             predicted_genotype_freqs =
               list(probs_from_trm(out$transitionRateMatrix)),
             opts = list(opts),
             lambda_scores = list(outp$lambda_scores),
             time_py_call = unname(time_py_call)
             )
  })["elapsed"]

  return(list(time_out = time_out, out = out))
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
  stopifnot(is.numeric(opts$reltol))
  stopifnot(is.character(opts$round_result))
  stopifnot(is.character(opts$return_lambda_scores))
  stopifnot(is.character(opts$pick_1se))
}



## matrix, list -> output
##  opts is the opts$MHN_python list
do_MHN_python <- function(x, opts) {
  ## Define a Python function that encapsulates all the logic
  ## This ensures everything runs in a properly scoped environment
  reticulate::py_run_string('
def run_mhn_analysis(x, colnames_x, opts):
    """Run MHN analysis in an isolated scope"""
    import mhn
    import numpy as np
    import pandas as pd
    from mhn.optimizers import Optimizer

    # Create optimizer
    opt = Optimizer(Optimizer.MHNType[opts["Type"]])
    opt.set_penalty(opt.Penalty[opts["Penalty"]])

    # Convert data to Python
    data_matrix = pd.DataFrame(np.array(x, dtype=np.int32), columns=colnames_x)
    opt.load_data_matrix(data_matrix)

    # Run CV if needed
    if opts["lambda_min"] != opts["lambda_max"]:
        # Set seed before CV if provided
        if opts["seed"] is None or pd.isna(opts["seed"]):
            np.random.seed()  # No seed provided, use system time
        else:
            np.random.seed(int(opts["seed"]))  # Set the provided seed

        # Convert string to boolean
        pick_1se_bool = opts["pick_1se"] == "True"

        if opts["return_lambda_scores"] == "False":
            cv_lambda = opt.lambda_from_cv(
                lambda_min=opts["lambda_min"],
                lambda_max=opts["lambda_max"],
                steps=opts["steps"],
                nfolds=opts["nfolds"],
                show_progressbar=opts["show_progressbar"] == "True",
                pick_1se=pick_1se_bool
            )
            lambda_scores = None
        else:
            cv_lambda, lambda_scores = opt.lambda_from_cv(
                lambda_min=opts["lambda_min"],
                lambda_max=opts["lambda_max"],
                steps=opts["steps"],
                nfolds=opts["nfolds"],
                show_progressbar=opts["show_progressbar"] == "True",
                return_lambda_scores=True,
                pick_1se=pick_1se_bool
            )
    else:
        cv_lambda = opts["lambda_min"]
        lambda_scores = None

    # Train
    opt.train(
        lam=cv_lambda,
        maxit=opts["maxit"],
        reltol=opts["reltol"],
        round_result=opts["round_result"] == "True"
    )

    # Return results
    return {
        "result": opt.result,
        "lambda_scores": lambda_scores
    }
  ')

  ## Convert options to a Python dictionary
  py_opts <- list(
    Type = opts$Type,
    Penalty = opts$Penalty,
    lambda_min = as.numeric(opts$lambda_min),
    lambda_max = as.numeric(opts$lambda_max),
    steps = as.integer(opts$steps),
    nfolds = as.integer(opts$nfolds),
    show_progressbar = opts$show_progressbar,
    maxit = as.integer(opts$maxit),
    reltol = as.numeric(opts$reltol),
    round_result = opts$round_result,
    return_lambda_scores = opts$return_lambda_scores,
    seed = opts$seed,  # Pass seed to Python
    pick_1se = opts$pick_1se
  )

  ## Call the Python function
  py_result <- reticulate::py$run_mhn_analysis(x, colnames(x), py_opts)

  if ((opts$return_lambda_scores == "False") ||
        (opts$lambda_min == opts$lambda_max)) {
    lambda_scores <- NA
  } else {
    lambda_scores <- py_result$lambda_scores
  }

  return(list(result = py_result$result,
              lambda_scores = lambda_scores))
}

library(codetools)
checkUsageEnv(env = .GlobalEnv)
