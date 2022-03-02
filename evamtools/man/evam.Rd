\name{evam}
\alias{evam}
\title{Runs the CPMs}
\usage{
evam(x,
     methods = c("CBN", "OT", "HESBCN", "MHN", "OncoBN", "MCCBN"),
     max_cols = 15,
     cores = detectCores(),
     paths_max = FALSE,
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
                     reg = c("bic", "aic", "loglik"),
                     silent = TRUE
                 ),
     oncobn_opts = list(
                     model = "DBN",
                     algorithm = "DP",
                     k = 3,
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

\item{paths_max}{If TRUE, return all paths to the maxim/maxima, with
  their log probabilities.}

\item{mhn_opts}{A list with two named arguments.
  \itemize{

 \item{lambda: }{The penalty for MHN. Defaults to 1/nrow(data). (These are not the lambdas as the estimated parameters for the rates of the continuous-time Markov chains for MHN or CBN or HESBCN.) }

 \item{omp_threads: }{Number of OMP threads for MHN. Do not pass
   \code{thrds} > 1 with \code{cores} > 1: as with MCCBN, do not use OpenMP threads from forked process from \code{\link[parallel]{mclapply}}. }
  }

}

\item{ot_opts}{A list with the single named argument \code{with_errors_dist_ot}: value for option with \code{with.errors} in the call to \code{\link[Oncotree]{distribution.oncotree}}. A value of TRUE means to incorporate the false positive and negative errors when returning the probabilities of genotypes under OT. Note that for large models using a value of TRUE can result in very long computing times. Default is TRUE.} 
\item{cbn_opts}{A named list with arguments passed to CBN. 
\itemize{
\item{omp_threads: }{OMP threads to be used by CBN (set via the
  environment variable OMP\_NUM\_THREADS). Defaults to 1. In contrast to
  MCCBN and MHN, you can set
  this to a number larger than one even if you set \code{cores} to a
  number larger than one (i.e., if we use
  \code{\link[parallel]{mclapply}}). It is unclear, though, that it will
speed things much or what is the best number of threads to use. }
\item{cbn_init_poset: }{Initial poset for CBN; one of "linear" or "OT" (default).}
}
}

\item{hesbcn_opts}{Named list of arguments used in the fit of
  H-ESBCN (details in \url{https://github.com/danro9685/HESBCN}). 
  \itemize{
\item{n_steps: }{Number of steps to run. Default: 100000}

\item{reg: }{Regularization: one of \code{bic} (default), \code{aic},
  \code{loglik}.} 

%% \item{tmp_dir: }{Directory name where the oput is located.}

\item{seed: }{Seed to run the experiment}

%\item{addname: }{String to append to the temporary folder name. Default NULL}

\item{silent: }{Whether to run show message showing the folder name where HESBCN is run}
  }
  }
  

 \item{oncobn_opts}{Named list of arguments used in the fit of OncoBN. See \code{\link[OncoBN]{fitCPN}}. } 


\item{mccbn_opts}{Named list of arguments used in the fit of
    MC-CBN. These are \code{model} (one of \code{OT-CBN} or
    \code{H-CBN2}). The rest are options passed to
    \code{\link[mccbn]{adaptive.simulated.annealing}}; see the help of
    \code{\link[mccbn]{adaptive.simulated.annealing}} for details.
    In addition, the following options:
    \itemize{
       \item{tmp_dir: }{Directory name where the oput is      located. This is  passed to      \code{\link[mccbn]{adaptive.simulated.annealing}}, as  argument \code{outdir}, with \code{addname} added, if provided. }
 
      \item{addname: }{String to append to the temporary directory name. Default is NULL.}

      \item{silent: }{Whether to show a message with the name of the
    directory where MCCBN is run. This \code{silen} is different from  \code{mccbn_hcbn2_opts$verbose}.}
    }
    %% \code{mccbn_hcbn_proc} (i.e., used only
    %% if fitting a \code{H-CBN2} model). See \code{\link{do_MCCBN_HCBN2}}.
    
    Note: do not pass \code{thrds} > 1 with \code{cores} > 1: as with
  MHN, do not use OpenMP threads from forked process from 
  \code{\link[parallel]{mclapply}}.
  }
}

\value{A list with named components (that should be self-explanatory). The pattern is method_component.
  
  \itemize{
        \item{OT_model: }{Data frame with parent and descendant edges,
  edge weight, and observed and predicted frequencies of genes.}
        \item{OT_f_graph: }{The fitness graph, as a sparse matrix, with weights  obtained from the edge weights (this is not a transition rate  matrix). See full documentation for details.}
        \item{OT_trans_mat: }{Transition matrix between genotypes. This is
  really an abuse of what an untimed OT provides. See full documentation for details.}
        \item{OT_predicted_genotype_freqs: }{Probabilities of genotypes from the
  OT model, as a data frame.}
	\item{OT_paths_max: }{If \code{paths_max} is TRUE, a list of two components, \code{paths} and \code{weights}. The \code{paths} list is a list of \code{igraph.vs} (igraph vertex sequences) objects, one for each path; the \code{weights} is vector of log-probabilities of each path. If \code{paths_max} is FALSE, the default, NA.}




        \item{CBN_model: }{Similar to the ouput from OT, but with
        lambdas. The lambda to be used is "rerun_lambda".}
        \item{CBN_trans_rate_mat: }{Transition rate matrix as a sparse matrix.}
        \item{CBN_trans_mat: }{Transition matrix between genotypes,
        obtained from the transition rate matrix using competing exponentials.}
        \item{CBN_td_trans_mat: }{Time-discretized transition matrix,
        using the uniformization method; see full documentation for details.}
      \item{CBN_predicted_genotype_freqs: }{Named vector of probabilities of genotypes predicted
        by the CBN model (under a model where sampling times are
        distributed as an exponential of rate 1).}
      \item{CBN_paths_max: }{As for OT.}

      
        \item{MCCBN_model: }{As for CBN, only with one column of \eqn{\lambda}s.}
        \item{MCCBN_trans_rate_mat: }{As for CBN.}
        \item{MCCBN_trans_mat: }{As for CBN.}
        \item{MCCBN_td_trans_mat: }{As for CBN.}
       \item{MCCBN_predicted_genotype_freqs: }{As for CBN.}
       \item{MCCBN_paths_max: }{As for OT.}

       
        \item{MHN_theta: }{Matrix of estimated thetas2.}
        \item{MHN_trans_rate_mat: }{As for CBN.}
        \item{MHN_trans_mat: }{As for CBN.}
        \item{MHN_td_trans_mat: }{As for CBN.}
        \item{MHN_exp_theta: }{Matrix of the exponential of thetas; the
        matrix \eqn{\Theta} in Schill et al. (just
        each theta, exponentiated; not the matrix exponential of the matrix of thetas).}
        \item{MHN_predicted_genotype_freqs: }{As for CBN.}
	\item{MHN_paths_max: }{As for OT.}

	
        \item{OncoBN_model: }{Similar to the ones above (but with a column
        named Thetas, instead of lambdas or edge weights), with the
        additional column dQuote{Relation}, that can take values OR (if
        fitting model DBN) or AND (if fitting model CBN); Single
        indicates nodes with a single ancestor (where OR or AND make no difference).}
        \item{OncoBN_likelihood: }{Likelihood of the OncoBN model.}
        \item{OncoBN_f_graph: }{As for OT.}
	\item{OncoBN_trans_mat: }{As for OT.}
        \item{OncoBN_predicted_genotype_freqs: }{As for OT.}
        \item{OncoBN_fitted_model: }{DBN or CBN, depending on what you chose.}
        \item{OncoBN_epsilon: }{Epsilon (this is an argument of the call
        to evam, but it is evaluated after possibly having modified the
        input data; see below).}
        \item{OncoBN_parent_set: }{.}
	\item{OncoBN_paths_max: }{As for OT.}

	
        \item{HESBCN_model: }{As for CBN.}
        \item{HESBCN_parent_set: }{As for CBN.}
        \item{HESBCN_trans_rate_mat: }{As for CBN.}
        \item{HESBCN_trans_mat: }{As for CBN.}
        \item{HESBCN_td_trans_mat: }{As for CBN.}
	\item{HESBCN_predicted_genotype_freqs: }{As for CBN.}
	\item{HESBCN_paths_max: }{As for OT.}
	
	
        \item{original_data: }{The original data.}
        \item{analyzed_data: }{The data that were actually analyzed. Can
        differ from the original data because of the data pre-processing
        steps.}
      \item{genotype_id_ordered: }{A named vector, from 1:number of genotypes, with names the genotypes. This can be useful for sorting; WT has value 1, and genotypes are ordered by increasing number of mutations and, withing number of mutations, alphanumerically.}
      \item{all_options: }{All of the options used, as a list of lists.}
  }
}


\note{For some methods, such as MHN and OncoBN, some parameters tipically depend on the data (lambda and epsilon for MHN and OncoBN, respectively). Since we first examine and possibly modify the input data, the values might not be the ones you thought you entered, as the \strong{options should be evaluated after the data are pre-processed}.

  The \strong{data pre-processing} involves, in sequence, these steps:


  \itemize{ 

  \item{Adding pseudosamples: }{If any gene (column) is always observed
    mutated (i.e, has a value of 1 for all observations), we add one
    observation that has no gene mutated.} 
  
  \item{Removing genes with no mutations: }{Any column that has value of 0 for all observations is removed.}

  \item{Merging identical columns: }{Any identical columns are replaced by a single one (with a new identifier, the result of pasting the names of the fused columns separated by a _). 

  \item{No more than max_cols: }{If the \dQuote{max_cols} argument is not NULL, and if the data set has more columns that max_cols, we keep only max_cols columns of data, those with the largest number of mutations.}} 
 }

\strong{Changing only some options:} Often, you will want to change only
some of the options. You can enter, in the list, only the options you
want changed (not the remaining ones). There is one example below.


During the execution, and as messages, the elapsed time of each
procedure is reported. This includes executing the model itself and
possible additional operations, such as obtaining the transition rate
matrix, etc. (So, for example, the time for estimating the
matrix of thetas for MHN is much smaller than the reported time here,
which also includes building the transition rate matrix).


You can \strong{convert paths and their probabilities} to a data frame easily. For example, using the non-exported simple function \code{pats_and_probs_to_df}. See example below. By default, we do not return paths, as their number can grow very quickly with number of genotypes and you only need if, well, you care about them.

}



\description{
Executes all CPMS given a cross sectional data set.
}

\references{

\bold{OT}

- Szabo, A., & Boucher, K. M. (2008). Oncogenetic Trees. In W. Tan, & L. Hanin
  (Eds.), Handbook of Cancer Models with Applications (pp. 1–24). : World
  Scientific.


- Desper, R., Jiang, F., Kallioniemi, O. P., Moch, H., Papadimitriou, C. H., &
  Sch\"affer, A A (1999). Inferring tree models for oncogenesis from comparative
  genome hybridization data. J Comput Biol, 6(1), 37–51.


\bold{CBN and MCCBN}

- Beerenwinkel, N., & Sullivant, S. (2009). Markov models for accumulating
  mutations. Biometrika, 96(3), 645.

- Gerstung, M., Baudis, M., Moch, H., & Beerenwinkel, N. (2009). Quantifying
  cancer progression with conjunctive Bayesian networks. Bioinformatics, 25(21),
  2809–2815. \url{http://dx.doi.org/10.1093/bioinformatics/btp505}


- Gerstung, M., Eriksson, N., Lin, J., Vogelstein, B., & Beerenwinkel,
  N. (2011). The Temporal Order of Genetic and Pathway Alterations in
  Tumorigenesis. PLoS ONE, 6(11),
  27136. \url{http://dx.doi.org/10.1371/journal.pone.0027136}


- Montazeri, H., Kuipers, J., Kouyos, R., B\"oni, J\"urg, Yerly, S., Klimkait, T., Aubert, V., … (2016). Large-scale inference of conjunctive Bayesian networks. Bioinformatics, 32(17), 727–735. \url{http://dx.doi.org/10.1093/bioinformatics/btw459}


\bold{MHN}


- Schill, R., Solbrig, S., Wettig, T., & Spang, R. (2020). Modelling cancer
  progression using Mutual Hazard Networks. Bioinformatics, 36(1),
  241–249. \url{http://dx.doi.org/10.1093/bioinformatics/btz513}

\bold{HESBCN (PMCE)}


- Angaroni, F., Chen, K., Damiani, C., Caravagna, G., Graudenzi, A., &
  Ramazzotti, D. (2021). PMCE: efficient inference of expressive models of cancer
  evolution with high prognostic power. Bioinformatics, 38(3),
  754–762. \url{http://dx.doi.org/10.1093/bioinformatics/btab717}

   (About terminology: we will often refer to HESBCN, as that is the program we use, as shown here: \url{https://github.com/danro9685/HESBCN}. H-ESBCN is part of the PMCE procedure).

\bold{OncoBN (DBN)}

- Nicol, P. B., Coombes, K. R., Deaver, C., Chkrebtii, O., Paul, S., Toland,
  A. E., & Asiaee, A. (2021). Oncogenetic network estimation with disjunctive
  Bayesian networks. Computational and Systems Oncology, 1(2),
  1027. \url{http://dx.doi.org/10.1002/cso2.1027}

\bold{Conditional prediction of genotypes and probabilities of paths from CPMs}

- Hosseini, S., Diaz-Uriarte, Ramon, Markowetz, F., & Beerenwinkel,
  N. (2019). Estimating the predictability of cancer evolution. Bioinformatics,
  35(14), 389–397. \url{http://dx.doi.org/10.1093/bioinformatics/btz332}

- Diaz-Uriarte, R., & Vasallo, C. (2019). Every which way? On predicting tumor
  evolution using cancer progression models. PLOS Computational Biology, 15(8),
  1007246. \url{http://dx.doi.org/10.1371/journal.pcbi.1007246}

- Diaz-Colunga, J., & Diaz-Uriarte, R. (2021). Conditional prediction of
  consecutive tumor evolution using cancer progression models: What genotype
  comes next? PLOS Computational Biology, 17(12),
  1009055. \url{http://dx.doi.org/10.1371/journal.pcbi.1009055}

 
}


\examples{
data(every_which_way_data)
Dat1 <- every_which_way_data[[16]][1:40, 2:6]
out1 <- evam(Dat1,
            methods = c("CBN", "OT", "OncoBN",
                        "MHN", "HESBCN"))


## Running only some methods and changing some options
## (this example is not necessarily sensible!)
## Of course, we must use the name of the data in an option that is
## data-dependent

out2 <- evam(Dat1,
            methods = c("CBN", "OT", "OncoBN",
                        "MHN", "HESBCN", "MCCBN"),
            mhn_opts = list(lambda = 5/nrow(Dat1)),
            oncobn_opts = list(model = "CBN"))


## Getting paths to maximum/maxima

out3 <- evam(Dat1,
            methods = c("CBN", "OT", "OncoBN",
                        "HESBCN"),
            paths_max = TRUE)


## Convert probabilities of paths into a data frame,
## with strings and the probs. of the paths.
## Same code as in evamtools:::paths_and_probs_to_df

f_prob_paths <- function(path_prob) {
    return(data.frame(Path = vapply(path_prob$paths,
                             function(x) paste(igraph::as_ids(x),
                                                collapse = " -> "), ""),
                      Prob = exp(path_prob$weights)))
}

f_prob_paths(out3$HESBCN_paths_max)


}
