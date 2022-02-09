\name{sample_to_pD_order}
\alias{sample_to_pD_order}
\title{Count genotypes 

Take a sample (a vector), with genotypes as "A, B", etc
and return a vector of frequencies (counts) in the exact same
order as used by MHN.}
\usage{
sample_to_pD_order(x, ngenes, gene_names = NULL)
}
\arguments{
\item{x}{vector of genotypes}

\item{ngenes}{total number of genes}

\item{gene_names}{List of gene names. If NULL, genes will be named alphabetically}
}
\value{
counts of all genotypes in same order as used by MHN
}
\description{
Count genotypes 

Take a sample (a vector), with genotypes as "A, B", etc
and return a vector of frequencies (counts) in the exact same
order as used by MHN.

gene_names is always sorted inside the function to ensure
results do not differ by gene_names order.

}
