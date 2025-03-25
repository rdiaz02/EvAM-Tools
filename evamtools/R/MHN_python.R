## Use the Python implementation of MHN

## I try to use Javier Lopez de Lema's code structure
## But since this might not yet be fully integrated into evam
## I may add a few things here that will later be moved out

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





do_MHN_python <- function(x,
                          dirname = NULL,
                          addname = NULL,
                          silent = TRUE,
                          rmfile = TRUE) {
  ## Setup the tmp directory structure, like for cbn
  ## Check the python executable is there



  if(is.null(dirname)) {
    dirname <- tempfile()
    dirname0 <- NULL
    if(!is.null(addname)) {
      dirname0 <- dirname
      dirname <- paste0(dirname, "/",
                        "_cbn_", init.poset, "_",
                        addname)
    }
    if(!silent)
      message(paste("\n Using dir", dirname))
    if(dir.exists(dirname)) {
      stop("dirname ", dirname, "exists")
    }
    dir.create(dirname, recursive = TRUE)
  }





  ## FIXME: adapt all this to MHN
  if(rmfile) {
    files.created <- paste(dirname, c(".pat", ".prf", ".log", ".poset",
                                      ".lambda"), sep = "")
    file.remove(files.created)
    try(file.remove(paste0(dirname, ".lambda_i"))) ## for the hack that produces lambda_i
    try(file.remove(paste0(dirname, ".logliknew"))) ## for the 2nd hack that produces logliknew
    ## to get the lambda, iterative, file
    file.remove(paste(dirname, "/00000.poset", sep = ""))
    file.remove(dirname)
    if(!is.null(dirname0))
      file.remove(dirname0)
  }
}
