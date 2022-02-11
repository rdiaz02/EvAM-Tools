\name{sample_CPMs}
\alias{sample_CPMs}
\title{Obtain genotype frequencies from the CPM models and, optionally,
  counts of genotype transitions.}
\usage{
sample_CPMs(cpm_output,  N,
                methods = c("OT", "OncoBN",
                            "CBN", "MCCBN",
                            "MHN", "HESBCN"),
                genotype_transitions = TRUE
)
}
\arguments{
\item{cpm_output}{Output from calling all_methods2trans_mat}

\item{N}{Number of samples to generate}

\item{methods}{List of methods that we want to sample}

\item{genotype_transitions}{If TRUE, return the number of transitions
  between pairs of genotypes. Not available for OT and OncoBN.}
}


\value{
  List with observed frequencies of genotypes for the methods used. If
  \code{genotype_transitions = TRUE}, vector of counts of genotype transitions.
}


\description{
 Obtain samples (of size \code{N}) for each method. For methods that
 return a continuous-time Markov Chain, we sample from the transition
 rate matrix, assuming that sampling times come from an exponential
 distribution of rate 1. For OT and OncoBN, we use multinomial sampling,
 of size \code{N}, with probabilities of genotypes as given by the
 method. 
 
  Genotype transitions, if requested, are obtained by counting the transitions between pairs of genotypes when simulating from the continuous-time Markov Chain.
}


\examples{
    data(every_which_way_data)
    Dat1 <- every_which_way_data[[16]][1:40, 2:6]
    out <- suppressMessages(evam(Dat1,
                                 methods = c("CBN", "OT", "OncoBN",
                                             "MHN", "HESBCN", "MCCBN")))
    outs <- sample_CPMs(out, N = 1000)
    outs_no_transit <- sample_CPMs(out, N = 1000,
                                       genotype_transitions = FALSE)
}
