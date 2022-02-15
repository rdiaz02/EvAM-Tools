\name{sample_CPMs}
\alias{sample_CPMs}
\title{Obtain samples of genotypes from the CPM models and, optionally,
  counts of genotype transitions.}
\usage{
sample_CPMs(cpm_output,  N,
                methods = c("OT", "OncoBN",
                            "CBN", "MCCBN",
                            "MHN", "HESBCN"),
                output = c("genotype_freqs")
)
}
\arguments{
\item{cpm_output}{Output from calling all_methods2trans_mat}

\item{N}{Number of samples to generate}

\item{methods}{List of methods that we want to sample}

\item{output}{A vector with one or more of the following possible outputs: genotype_freqs, obs_genotype_transitions, state_counts. Even if requested, obs_genotype_transitions and state_counts are not available for OT and OncoBN. } 

}


\value{
  A list, with a many entries as methods * number of components
  requested. If all, then for each method among CBN, MCCBN, HESBCN, and MHN:
   \itemize{
    \item{genotype_freqs:}{frequencies of genotypes (final genotypes from  the transitions)}.
    \item{obs_genotype_transitions:}{Number of observed transitions between   genotypes (as a sparse matrix).}
    \item{state_counts}{Number of times each genotype was visited during
      the transitions. Column sums of observed transitions are equal to state counts.}
  }
  
  Observed genotype transitions are not the way to obtain estimates
  of  transition probabilities. The transition probabilities given by
  each   method are already available from the output of evam itself. These  genotype transitions are the observed transitions during the
  simulation of the sampling process and, thus, have additional noise.

  For OT and OncoBN, only the \code{genotype_freqs} component is
  available (the other  two are not available).
  
}


\description{
  Obtain samples of genotypes from the CPM models and, optionally,
  counts of genotype transitions.

  If you request only genotype_freqs, obtain samples (of size \code{N}) for each method using the predicted genotype frequencies. If you request also obs_genotype_transitions and state_counts, for OT and OncoBN you still only get samples (of size \code{N}) for each method using the predicted genotype frequencies. For the remaining methods, as they return the transition rate of a continuous-time Markov Chain, we sample from the transition rate matrix, assuming that sampling times come from an exponential  distribution of rate 1. 
 
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

## Sample from the predicted genotype frequencies
## for all methods in the output out

outS1 <- sample_CPMs(out, N = 1000)

## Only CBN and will simulate sampling from the transition
## rate matrix. 

outS2 <- sample_CPMs(out, N = 1000, methods = "CBN",
                     output = "obs_genotype_transitions")


## No output available for OT
## For CBN and MHN simulate from the transition rate matrix

outS3 <- sample_CPMs(out, N = 1000, methods = c("CBN", "OT", "MHN"), 
                     output = c("obs_genotype_transitions",
                                "state_counts"))


## OT sampled from the predicted genotype frequencies
## No obs_genotype_transitions available for OT
## CBN and OT simulate from the transition rate matrix, for consistency

outS4 <- sample_CPMs(out, N = 1000, methods = c("CBN", "OT", "MHN"), 
                     output = c("obs_genotype_transitions",
                                "genotype_freqs"))
}
