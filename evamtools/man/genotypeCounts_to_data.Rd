\name{genotypeCounts_to_data}
\alias{genotypeCounts_to_data}
\title{Obtain a matrix of 0/1 with subjects as rows and genes as columns
  from a named vector with genotype frequencies. Optionally, add
  observation error.}
\usage{
genotypeCounts_to_data(x, e = 0.01)

}
\arguments{
\item{x}{Named vector that contains the counts of each genotype. This vector is obtained, for example, from \code{\link{sample_CPMs}}.}
\item{e}{The proportion of observations with error (for
    instance, genotyping error). If larger than 0, this proportion of
    entries in the sampled matrix will be flipped (i.e., 0s turned
    to 1s and 1s turned to 0s).}
}


\value{A matrix. Each row is a "sampled genotype", where 0 denotes no
  alteration   and 1 alteration. If \code{e > 0}, each observed genotype
  might have been altered by random noise, where the probability of
  turning a 1 into a 0 is the same as turning a 0 into a 1.}

\seealso{
  \code{\link{sample_CPMs}}, \code{\link{generate_random_evam}}
 
}

\description{Obtain a matrix of 0/1 with subjects as rows and genes as columns
  from a data frame with genotype frequencies. Optionally, add
  observation error.}

\examples{
rmhn <- generate_random_evam(model = "MHN", ngenes = 5)
sample_mhn <- sample_CPMs(rmhn, N = 12)

## No noise.
d1 <- genotypeCounts_to_data(sample_mhn$MHN_sampled_genotype_freqs, 0)

## Compare
d1
sample_mhn$MHN_sampled_genotype_freqs

## A lot of noise added
(d2 <- genotypeCounts_to_data(sample_mhn$MHN_sampled_genotype_freqs, 0.3))

}
