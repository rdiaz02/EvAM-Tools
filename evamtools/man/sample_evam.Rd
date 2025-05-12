\name{sample_evam}
\alias{sample_evam}
\alias{sample_CPMs}
\title{Obtain samples of genotypes from the EvAM (CPM) models and, optionally,
  counts of genotype transitions.}
\usage{
sample_evam(cpm_output,  N,
                methods = NULL,
                output = c("sampled_genotype_counts"),
                obs_noise = 0,
                genotype_freqs_as_data = TRUE
            )

sample_CPMs(cpm_output,  N,
                methods = NULL,
                output = c("sampled_genotype_counts"),
                obs_noise = 0,
                genotype_freqs_as_data = TRUE
            )
}
\arguments{
\item{cpm_output}{Output from calling all_methods2trans_mat}

\item{N}{Number of samples to generate}

\item{methods}{Vector of strings with the names of methods that we want to sample. If NULL, all
  methods with output in cpm_output. The list of available methods is
  OT, OncoBN, CBN, MCCBN, MHN, HESBCN.}

\item{output}{A vector with one or more of the following possible
  outputs: sampled_genotype_counts, obs_genotype_transitions,
  state_counts. Even if requested, obs_genotype_transitions and
  state_counts are not available for OT and OncoBN. }

\item{obs_noise}{When obtaining a sample, should we add observation
  noise (for example, genotyping error) to the data?   If larger than 0,
  this \code{obs_noise} proportion of entries in the sampled matrix
  will be flipped (i.e., 0s turned to 1s and 1s turned to 0s).}


\item{genotype_freqs_as_data}{If TRUE, return a matrix where each row is a "sampled genotype", where 0 denotes no alteration and 1 alteration in the gene of the corresponding column.}

}


\value{A list, with a many entries as methods times number of components requested. For each method among CBN, MCCBN, HESBCN, and MHN:

    \item{sampled_genotype_counts}{Counts, or absolute genotype frequencies,
  obtained by sampling from the predicted frequencies. See also
  Description, below.}
    \item{obs_genotype_transitions}{Number of observed transitions between   genotypes (as a sparse matrix).}
    \item{state_counts}{Number of times each genotype was visited
  during the transitions. Column sums of observed genotype transitions are equal to state counts.}
    \item{sampled_genotype_counts_as_data}{The genotypes in a matrix of 0/1. This can directly be passed as an argument for \code{\link{evam}}, as the input data.}



  Observed genotype transitions are not the way to obtain estimates of transition probabilities. The transition probabilities given by each method are already available from the output of evam itself. These genotype transitions are the observed transitions during the simulation of the sampling process and, thus, have additional noise.

  For OT and OncoBN, only the \code{sampled_genotype_counts} and
  \code{sampled_genotype_counts_as_data} components are available (the other two are not available).

}

\note{
\strong{sample_CPMs has been deprecated}. Use sample_evam.
}

\description{
  Obtain samples of genotypes from the CPM models and, optionally, counts of genotype transitions.

  For OT and OncoBN we always obtain the absolute genotype frequencies
  by drawing samples of size \code{N}, with replacement, using as
  probabilities the predicted genotype frequencies.

  For the remaining methods, that is also what we do, unless you request
  also obs_genotype_transitions and state_counts. In this case, since we
  need to simulate sampling from the continuous-time Markov Chain (with
  transition rates given by the transition rate matrix) to obtain
  state counts and observed genotype transitions, we use this same
  sampling to obtain the absolute genotype frequencies. (The results
  are, of course, equivalent, but sampling directly from the predicted
  frequencies is much faster). Note that the option to request
  obs_genotype_transitions was removed from the web app, as it was
  rarely used, but lead to confusion and could increase without good
  reason running times. So, from the web app, we sample
  without using obs_genotype_transitions.

  Observed genotype transitions, if requested, are obtained by counting the transitions between pairs of genotypes when simulating from the continuous-time Markov Chain. State counts are also obtained by counting from this process how many times a genotype was visited.
}

\seealso{
  \code{\link{random_evam}}

}

\examples{

data(every_which_way_data)
## Note that the last column of data is all 0's, so it will be dropped
Dat1 <- every_which_way_data[[16]][1:40, 2:6]
## For faster execution, use only some methods
out <- suppressMessages(evam(Dat1,
                             methods = c("CBN", "OT", "OncoBN",
                                         "MHN")))

## Sample from the predicted genotype frequencies
## only for OT
outS1_ot <- sample_evam(out, N = 1000, methods = "OT")

## Sample from the predicted genotype frequencies
## for OT and HESBCN. But the later was not in the output
## so we get a warning-
outS1_ot_2 <- sample_evam(out, N = 1000, methods = c("OT", "HESBCN"))


## Sample from the predicted genotype frequencies
## for all methods in the output out

outS1 <- sample_evam(out, N = 1000)

## Same, but adding observation error
outS1e <- sample_evam(out, N = 1000, obs_noise = 0.1)

## Only CBN and will simulate sampling from the transition
## rate matrix.

outS2 <- sample_evam(out, N = 1000, methods = "CBN",
                     output = "obs_genotype_transitions")


## No output available for OT
## For CBN and MHN simulate from the transition rate matrix

outS3 <- sample_evam(out, N = 1000, methods = c("CBN", "OT", "MHN"),
                     output = c("obs_genotype_transitions",
                                "state_counts"))


## OT sampled from the predicted genotype frequencies
## No obs_genotype_transitions available for OT
## CBN and OT simulate from the transition rate matrix, for consistency

outS4 <- sample_evam(out, N = 1000, methods = c("CBN", "OT", "MHN"),
                     output = c("obs_genotype_transitions",
                                "sampled_genotype_counts"))

## Only CBN,  will simulate sampling from the transition
## rate matrix and add observation error to the genotype frequencies.

outS5 <- sample_evam(out, N = 1000, methods = "CBN",
                     output = c("obs_genotype_transitions", "sampled_genotype_counts"), obs_noise = 0.1)


## A few HyperTraPS examples.
## This is an example; length is probably too short
outht <- evam(Dat1,
              methods = c("HyperTraPS"),
              hyper_traps_opts = list(length = 2,
                                      model = -1,
                                      seed = -1))

cbind(cond = outht$HyperTraPS_conditional_genotype_freqs,
      pred = outht$HyperTraPS_predicted_genotype_freqs)


## Recall last column of Dat1 is all 0s, so it is dropped.
## Thus, prob.set must be of length 5. Obtain
## predicted probability only at time 0
outht_t1 <- evam(Dat1,
                 methods = c("HyperTraPS"),
                 hyper_traps_opts = list(length = 2,
                                         model = -1,
                                         seed = -1,
                                         prob.set = c(1, rep(0, 4))))

## conditional should be similar to above, but predicted
## should only show WT
cbind(cond = outht_t1$HyperTraPS_conditional_genotype_freqs,
      pred = outht_t1$HyperTraPS_predicted_genotype_freqs)

stopifnot(all.equal(outht_t1$HyperTraPS_predicted_genotype_freqs,
                    c(1, rep(0, 15)), check.attributes = FALSE))

## Of course, we can obtain a sample at a different time without rerunning
## All sampled at 0 mutations
hypertrapsct::state.probs(outht$HyperTraPS_primary_output, c(1, rep(0, 4)))
## Uniform sampling
hypertrapsct::state.probs(outht$HyperTraPS_primary_output, NA)


outHT_SG <- sample_evam(outht, N = 1000)
outHT_SG$HyperTraPS_sampled_genotype_counts
}
