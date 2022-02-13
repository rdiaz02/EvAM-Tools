\name{evam}
\alias{evam}
\title{Runs the CPMs}
\usage{
evam(x,
     methods = c("CBN", "OT", "HESBCN", "MHN", "OncoBN", "MCCBN"),
     max_cols = 15,
     cores = detectCores(),
     mhn_opts = list(lambda = 1/nrow(x),
                     omp_threads = ifelse(cores > 1, 1, detectCores())),
     ot_opts = list(with_errors_dist_ot = TRUE),
     cbn_opts = list(
                     omp_threads = ifelse(cores > 1, 1, detectCores()),
                     init_poset = "OT"
                 ),
     hesbcn_opts = list(
                     steps = 100000,
                     seed = NULL,
                     reg = c("bic", "aic", "loglik")
                 ),
     oncobn_opts = list(
                     model = "DBN",
                     algorithm = "DP",
                     epsilon = min(colMeans(x)/2),
                     silent = TRUE
                 ),
     mccbn_opts = list(
                     model = "OT-CBN",
                     tmp_dir = NULL,
                     addname = NULL,
                     silent = TRUE,
                     L = 100,
                     sampling = c("forward", "add-remove",
                               "backward", "bernoulli", "pool"),
                     max.iter = 100L,
                     update.step.size = 20L,
                     tol = 0.001,
                     max.lambda.val = 1e+06,
                     T0 = 50,
                     adap.rate = 0.3,
                     acceptance.rate = NULL,
                     step.size = NULL,
                     max.iter.asa = 10000L,
                     neighborhood.dist = 1L,
                     adaptive = TRUE,
                     thrds = 1L,
                     verbose = FALSE,
                     seed = NULL)
)
}
\arguments{

\item{x}{cross sectional data}

\item{methods}{Methods to use. A vector with one or more of the
  following strings, \dQuote{OT}, \dQuote{OncoBN}, \dQuote{CBN}, \dQuote{MCCBN}, \dQuote{MHN}, \dQuote{HESBCN}.}


\item{max_cols}{Maximum number of columns to use in the analysis. If x has > max_cols, selected columns are those with the largest number of events.} 

\item{cores}{If larger than 1, use mclapply to run all methods. This is
  the default. If you use mclapply, MHN and MCCBN should not use OMP
  (i.e., the number of threads for OMP for MHN and MCCBN should be 1).}

\item{mhn_opts}{A list with the named argument \code{lambda}: the
  penalty for MHN. Defaults to 1/nrow(data). (These are not the lambdas
  as the estimated parameters for the rates of the continuous-time
  Markov chains for MHN or CBN or HESBCN.). \code{omp_threads}: OMP threads
  for MHN.} 


\item{ot_opts}{A list with the single named argument \code{with_errors_dist_ot}: value for option with \code{with.errors} in the call to \code{\link[Oncotree]{distribution.oncotree}}. A value of TRUE means to incorporate the false positive and negative errors when returning the probabilities of genotypes under OT. Note that for large models using a value of TRUE can result in very long computing times. Default is TRUE.} 


\item{cbn_opts}{A named list with arguments passed to CBN. 
\itemize{
\item{omp_threads:}{OMP threads to be used by CBN (set via the
  environment variable OMP\_NUM\_THREADS). Defaults to 1. In contrast to
  MCCBN and MHN, you can set
  this to a number larger than one even if you set \code{cores} to a
  number larger than one (i.e., if we use
  \code{\link[parallel]{mclapply}}). It is unclear, though, that it will
speed things much or what is the best number of threads to use. }
\item{cbn_init_poset:}{Initial poset for CBN; one of "linear" or "OT" (default).}
}
}

 \item{hesbcn_opts}{Named list of arguments used in the fit of H-ESBCN. See \code{\link{do_HESBCN}}. } 

 \item{oncobn_opts}{Named list of arguments used in the fit of OncoBN. See \code{\link{do_OncoBN}}. } 



\item{mccbn_opts}{Named list of arguments used in the fit of
    MC-CBN. These are \code{model} (one of \code{OT-CBN} or
    \code{H-CBN2}); the
    rest are arguments passed to \code{mccbn_hcbn_proc} (i.e., used only
    if fitting a \code{H-CBN2} model). See \code{\link{do_MCCBN_HCBN2}}.

    Note: do not pass \code{thrds} > 1 with \code{cores} > 1: as with
  MHN, do not use OpenMP threads from forked process from 
  \code{\link[parallel]{mclapply}}.
  }
}

\description{
Executes all CPMS given a cross sectional data set.
}


\examples{
    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:6]
    out <- suppressMessages(evam(Dat1,
                                 methods = c("CBN", "OT", "OncoBN",
                                             "MHN", "HESBCN", "MCCBN")))
}
