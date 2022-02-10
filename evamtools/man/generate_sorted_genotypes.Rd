\name{generate_sorted_genotypes}
\alias{generate_sorted_genotypes}
\title{List of sorted genotypes}
\usage{
generate_sorted_genotypes(n_genes, gene_names = NULL)
}
\arguments{
\item{n_genes}{Integer > 0}

\item{gene_names}{Vector of strings
If NULL, defaults gene_names will the firt n_genes letters of the alphabet}
}
\value{
Vector with the sorted genotypes
}
\description{
  Returns all sorted genotypes for a given number of genes
   with sorting the same as given by sample_to_pD_order.
    gene_names is always sorted inside the function to
    ensure results do not differ by gene_names order.
}
