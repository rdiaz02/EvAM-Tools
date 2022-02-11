\name{evam}
\alias{evam}
\title{Runs the CPMs}
\usage{
evam(x,
     methods = c("CBN", "OT", "HESBCN", "MHN", "OncoBN", "MCCBN"),
     max_cols = 15,
     mhn_opts = list(lambda = 1/nrow(x)),
                 ot_opts = list(with_errors_dist_ot = TRUE),
                 cbn_opts = list(
                     cores = 1,
                     init_poset = "OT"
                 ),
     hesbcn_opts = list(
                     steps = 100000,
                     seed = NULL
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

\item{methods}{Methods to use. For now, the only thing that matters is whether the vector includes the strings "MCCBN" and "OncoBN".}


\item{max_cols}{Maximum number of columns to use in the analysis. If x has > max_cols, selected columns are those with the largest number of events.} 


\item{mhn_opts}{A list with the single named argument \code{lambda}: the
  penalty for MHN. Defaults to 1/nrow(data). (These are not the lambdas
  as the estimated parameters for the rates of the continuous-time
  Markov chains for MHN or CBN or HESBCN.)} 


\item{ot_opts}{A list with the single named argument \code{with_errors_dist_ot}: value for option with \code{with.errors} in the call to \code{\link[Oncotree]{distribution.oncotree}}. A value of TRUE means to incorporate the false positive and negative errors when returning the probabilities of genotypes under OT. Note that for large models using a value of TRUE can result in very long computing times. Default is TRUE.} 


\item{cbn_opts}{A named list with arguments passed to CBN. 
\itemize{
\item{cores:}{Cores to be used by CBN (set via the environment variable OMP\_NUM\_THREADS). Defaults to 1.}
\item{cbn_init_poset:}{Initial poset for CBN; one of "linear" or "OT" (default).}
}
}

 \item{oncobn_opts}{Named list of arguments used in the fit of OncoBN. See \code{\link{do_OncoBN}}. } 


 \item{hesbcn_opts}{Named list of arguments used in the fit of H-ESBCN. See \code{\link{do_HESBCN}}. } 

\item{mccbn_opts}{Named list of arguments used in the fit of
    MC-CBN. These are \code{model} (one of \code{OT-CBN} or
    \code{H-CBN2}); the
    rest are arguments passed to \code{mccbn_hcbn_proc} (i.e., used only
    if fitting a \code{H-CBN2} model). See \code{\link{do_MCCBN_HCBN2}}.
  }
}

\description{
Executes all CPMS given a cross sectional data set
}


\examples{
    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:6]
    out <- suppressMessages(evam(Dat1,
                                 methods = c("CBN", "OT", "OncoBN",
                                             "MHN", "HESBCN", "MCCBN")))
}
