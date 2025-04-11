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
  ## mhn <- import("mhn")
  ## mhn_version <- py_get_attr(mhn, "__version__")
  ## message("Using MHN version ", mhn_version,
  ##         " with Python ", python_version)
  ##   py_run_string("from mhn.optimizers import Optimizer")
  ##   ## For data handling
  ##   py_run_string("
  ## import numpy as np
  ## import pandas as pd
  ## ")
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
  stopifnot(is.logical(opts$verbose))
  stopifnot(is.character(opts$round_result))
  stopifnot(is.character(opts$return_lambda_scores))
}


## Sensible opts
## opts0_MHN_python <- list(
##   omp_threads = 1, ##
##   Type = "cMHN", ## "oMHN" won't be used for now
##   Penalty = "SYM_SPARSE", ## SYM_SPARSE, L1 or L2. We'll use SYM and L1
##   seed = NA,
##   steps = 10,
##   nfolds = 5,
##   maxit = 5000,
##   reltol = 1e-07,
##   round_result = "True", ## Rounding was used in the 2020 paper
##                          ## and it is the default in the Python examples
##   show_progressbar = "False",
##   verbose = TRUE,
##   return_lambda_scores = "False"
## )


reset_python_env <- function() {
  py_run_string('
# Get all user-defined variables (skip built-ins and modules)
user_vars = [var for var in globals().keys()
             if not var.startswith("__")
                and not var in ("pd", "np", "mhn", "Optimizer")
                and not callable(globals()[var])]

# Delete all user variables
for var in user_vars:
    del globals()[var]
')
}




do_MHN_python <- function(x, opts) {
  if (opts$verbose) message("Starting.")

  # Define a Python function that encapsulates all the logic
  # This ensures everything runs in a properly scoped environment
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
        if opts["return_lambda_scores"] == "False":
            cv_lambda = opt.lambda_from_cv(
                lambda_min=opts["lambda_min"],
                lambda_max=opts["lambda_max"],
                steps=opts["steps"],
                nfolds=opts["nfolds"],
                show_progressbar=opts["show_progressbar"] == "True"
            )
            lambda_scores = None
        else:
            cv_lambda, lambda_scores = opt.lambda_from_cv(
                lambda_min=opts["lambda_min"],
                lambda_max=opts["lambda_max"],
                steps=opts["steps"],
                nfolds=opts["nfolds"],
                show_progressbar=opts["show_progressbar"] == "True",
                return_lambda_scores=True
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

  # Convert options to a Python dictionary
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
    return_lambda_scores = opts$return_lambda_scores
  )

  # Call the Python function
  py_result <- reticulate::py$run_mhn_analysis(x, colnames(x), py_opts)

  # Return the result in the format expected
  list(
    result = py_result$result,
    lambda_scores = if (opts$return_lambda_scores == "False") NA else py_result$lambda_scores
  )
}

library(codetools)
checkUsageEnv(env = .GlobalEnv)
