\name{evam}
\alias{evam}
\title{Runs the CPMs}
\usage{
evam(x, methods = c("CBN", "OT", "HESBCN", "MHN", "OncoBN"),
     max_cols = 15,
     cbn_cores = 1,
     cbn_init_poset = "OT",
     hesbcn_steps = 100000,
     hesbcn_seed = NULL,
     mhn_lambda = 1/nrow(x),
     oncobn_model = "DBN",
     oncobn_algorithm = "DP",
     oncobn_epsilon = min(colMeans(x)/2),
     oncobn_silent = TRUE,
     ot_with_errors_dist_ot = TRUE
)
}
\arguments{
\item{x}{cross sectional data}

\item{methods}{Methods to use. For now, the only thing that matters is
  whether the vector includes the strings "MCCBN" and "OncoBN". The default
  is not to use MCCBN, so c("CBN", "OT", "HESBCN", "MHN", "OncoBN")}

\item{max_cols}{Maximum number of columns to use in the analysis. If x
  has > max_cols, selected columns are those with the largest number of events.}

\item{cbn_cores}{Cores to be used by CBN (set via the environment
  variable OMP\_NUM\_THREADS. Defaults to 1.}
\item{cbn_init_poset}{Initial poset for CBN; one of "linear" or "OT"
  (default).}
\item{hesbcn_steps}{Number of steps for H-ESBCN; default is 100000}.
\item{hesbcn_seed}{Seed used for H-ESBCN; default is NULL. If you need
  reproducible results, fix it to some number, and note this differs
  from the R seed.}
\item{mhn_lambda}{Lambda (penalty) for MHN. Defaults to 1/nrow(data).}

\item{oncobn_model}{One of "DBN" (default) or "CBN".}

\item{oncobn_algorithm}{"DP" (default) or "GA". See \code{\link[OncoBN]{fitCPN}}.}

\item{oncobn_epsilon}{Penalty term. See \code{\link[OncoBN]{fitCPN}}.}

\item{oncobn_silent}{Whether to show messages produced by OncoBN.}

\item{ot_with_errors_dist_ot}{Value for option with \code{with.errors} in
the call to \code{\link[Oncotree]{distribution.oncotree}}. A value of
TRUE means to incorporate the false positive and negative errors when
returning the probabilities of genotypes under OT. Note
that for large models using a value of TRUE can result in very long
computing times. Default is TRUE.}
}

\description{
Executes all CPMS given a cross sectional data set
}
