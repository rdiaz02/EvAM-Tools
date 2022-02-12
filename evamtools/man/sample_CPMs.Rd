\name{sample_CPMs}
\alias{sample_CPMs}
\title{Obtain genotype frequencies from the CPM models and, optionally,
  counts of genotype transitions.}
\usage{
sample_CPMs(cpm_output,  N,
                methods = c("OT", "OncoBN",
                            "CBN", "MCCBN",
                            "MHN", "HESBCN"),
                obs_genotype_transitions = TRUE
)
}
\arguments{
\item{cpm_output}{Output from calling all_methods2trans_mat}

\item{N}{Number of samples to generate}

\item{methods}{List of methods that we want to sample}

\item{obs_genotype_transitions}{If TRUE, return the number of observed transitions between pairs of genotypes. Not available for OT and OncoBN.}
}


\value{
  A list, with a many entries as methods * number of components
  requested. If all, then for each method among CBN, HESBCN, and MHN:
   \itemize{
    \item{genotype_freqs:}{frequencies of genotypes (final genotypes from      the transitions)}.
    \item{obs_genotype_transitions:}{Number of observed transitions between
      genotypes (as a sparse matrix).}
    \item{state_counts}{Number of times each genotype was visited during
      the transitions. Column sums of observed transitions are equal to state counts.}
  }
  
  If \code{obs_genotype_transitions = FALSE},
  \code{obs_genotype_transitions} will not be returned.

  Observed genotype transitions are not to best way to obtain estimates
  of  transition probabilities. The transition probabilities given by
  each   method are already available from the output of evam itself. These  genotype transitions are the observed transitions during the
  simulation of the sampling process and, thus, have additional noise.

  For OT and OncoBN, only the \code{genotype_freqs} component (the other
  two are not available).
  
}


\description{
 Obtain samples (of size \code{N}) for each method. For methods that
 return a continuous-time Markov Chain, we sample from the transition
 rate matrix, assuming that sampling times come from an exponential
 distribution of rate 1. For OT and OncoBN, we use multinomial sampling,
 of size \code{N}, with probabilities of genotypes as given by the
 method. 
 
  Observed genotype transitions, if requested, are obtained by counting
  the transitions between pairs of genotypes when simulating from the
  continuous-time Markov Chain. State counts are also obtained by
  counting from this process how many times a genotype was visited.
}


\examples{
    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:6]
    out <- suppressMessages(evam(Dat1,
                                 methods = c("CBN", "OT", "OncoBN",
                                             "MHN", "HESBCN", "MCCBN")))
    outs <- sample_CPMs(out, N = 1000)
    outs_no_transit <- sample_CPMs(out, N = 1000,
                                       obs_genotype_transitions = FALSE)
}
