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

\value{A list with named components (that should be self-explanatory). The pattern is method_component.
  
  \itemize{
        \item{OT_model}{Data frame with parent and descendant edges,
  edge weight, and observed and predicted frequencies of genes.}
        \item{OT_f_graph}{The fitness graph, as a sparse matrix, with weights  obtained from the edge weights (this is not a transition rate  matrix). See full documentation for details.}
        \item{OT_trans_mat}{Transition matrix between genotypes. This is
  really an abuse of what an untimed OT provides. See full documentation for details.}
        \item{OT_genots_predicted}{Probabilities of genotypes from the
  OT model, as a data frame.}

        \item{CBN_model}{Similar to the ouput from OT, but with
        lambdas. The lambda to be used is "rerun_lambda".}
        \item{CBN_trans_rate_mat}{Transition rate matrix as a sparse matrix.}
        \item{CBN_trans_mat}{Transition matrix between genotypes,
        obtained from the transition rate matrix using competing exponentials.}
        \item{CBN_td_trans_mat}{Time-discretized transition matrix,
        using the uniformization method; see full documentation for details.}
      \item{CBN_genots_predicted}{Named vector of probabilities of genotypes predicted
        by the CBN model (under a model where sampling times are
        distributed as an exponential of rate 1).}
      
        \item{MCCBN_model}{As for CBN, only with one column of $\lambda$s.}
        \item{MCCBN_trans_rate_mat}{As for CBN.}
        \item{MCCBN_trans_mat}{As for CBN.}
        \item{MCCBN_td_trans_mat}{As for CBN.}
       \item{MCCBN_genots_predicted}{As for CBN.}
      
        \item{MHN_theta}{Matrix of estimated thetas2.}
        \item{MHN_trans_rate_mat}{As for CBN.}
        \item{MHN_trans_mat}{As for CBN.}
        \item{MHN_td_trans_mat}{As for CBN.}
        \item{MHN_exp_theta}{Matrix of the exponential of thetas; the
        matrix $\Theta$ in Schill et al. (just
        each theta, exponentiated; not the matrix exponential of the matrix of thetas).}
        \item{MHN_genots_predicted}{As for CBN.}
      
        \item{OncoBN_model}{Similar to the ones above (but with a column
        named Thetas, instead of lambdas or edge weights), with the
        additional column dQuote{Relation}, that can take values OR (if
        fitting model DBN) or AND (if fitting model CBN); Single
        indicates nodes with a single ancestor (where OR or AND make no difference).}
        \item{OncoBN_likelihood}{Likelihood of the OncoBN model.}
        \item{OncoBN_f_graph}{As for OT.}
	\item{OncoBN_trans_mat}{As for OT.}
        \item{OncoBN_genots_predicted}{As for OT.}
        \item{OncoBN_fitted_model}{DBN or CBN, depending on what you chose.}
        \item{OncoBN_epsilon}{Epsilon (this is an argument of the call
        to evam, but it is evaluated after possibly having modified the
        input data; see below).}
        \item{OncoBN_parent_set}{.}

        \item{HESBCN_model}{As for CBN.}
        \item{HESBCN_parent_set}{As for CBN.}
        \item{HESBCN_trans_rate_mat}{As for CBN.}
        \item{HESBCN_trans_mat}{As for CBN.}
        \item{HESBCN_td_trans_mat}{As for CBN.}
	\item{HESBCN_genots_predicted}{As for CBN.}
	
        \item{original_data}{The original data.}
        \item{analyzed_data}{The data that were actually analyzed. Can
        differ from the original data because of the data pre-processing
        steps.}
  }
}


\note{For some methods, such as MHN and OncoBN, some parameters tipically depend on the data (lambda and epsilon for MHN and OncoBN, respectively). Since we first examine and possibly modify the input data, the values might not be the ones you thought you entered, as the options should be evaluated after the data are processed.

  The \strong{data pre-processing} involves, in sequence, these steps:


  \itemize{ 

  \item{Adding pseudosamples}{If any gene (column) is always observed
    mutated (i.e, has a value of 1 for all observations), we add one
    observation that has no gene mutated.} 
  
  \item{Removing genes with no mutations}{Any column that has value of 0 for all observations is removed.}

  \item{Merging identical columns}{Any identical columns are replaced by a single one (with a new identifier, the result of pasting the names of the fused columns separated by a _). 

  \item{No more than max_cols}{If the \quote{max_cols} argument is not NULL, and if the data set has more columns that max_cols, we keep only max_cols columns of data, those with the largest number of mutations.}} 
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
