\name{do_MCCBN_HCBN2}
\alias{do_MCCBN_HCBN2}
\title{Run MCCBN with H-CBN2.
Using H-CBN2, as in example from \url{https://github.com/cbg-ethz/MC-CBN}}
\usage{
do_MCCBN_HCBN2(x,
               mccbn_hcbn2_opts = list(
                                       tmp_dir = NULL,
                                       addname = NULL,
                                       silent = TRUE,
                                       L = 100,
                                       sampling = c("forward",
                                                   "add-remove",
                                                   "backward",
                                                   "bernoulli",
                                                   "pool"),
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
  \item{mccbn_hcbn2_opts}{Options passed to
    \code{\link[mccbn]{adaptive.simulated.annealing}}. See the help of
    \code{\link[mccbn]{adaptive.simulated.annealing}} for details.
    In addition, the following options:
    
    \itemize{
      \item\code{tmp_dir}. Directory name where the oput is      located. This is  passed to      \code{\link[mccbn]{adaptive.simulated.annealing}}, as  argument \code{outdir}, with \code{addname} added, if provided. 
 
      \item\code{addname}. String to append to the temporary directory name. Default is NULL.

      \item\code{silent}.Whether to run show message showing the directory name  where MCCBN is run. This is different from  \code{mccbn_hcbn2_opts$verbose}.
    }
  }
  
  \item{x}{Cross sectional data. Matrix of genes (columns) and individuals (rows)}
}

\value{
A list with the adjacency matrix, the lambdas, the parent set
and a data.frame with From-To edges and associated lambdas.
}

\description{
Run MCCBN with H-CBN2

Using H-CBN2, as in example from \url{https://github.com/cbg-ethz/MC-CBN}.

L is not given a default value in the call; using here the one shown in
the web page example.
}
