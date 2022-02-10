\name{do_OncoBN}
\alias{do_OncoBN}
\title{Run OncoBN}
\usage{
do_OncoBN(
  data,
  model = "DBN",
  algorithm = "DP",
  epsilon = min(colMeans(data)/2),
  silent = TRUE)
}
\arguments{
\item{data}{Cross secitonal data. Matrix of genes (columns)
and individuals (rows)}

\item{model}{One of "DBN" (default) or "CBN".}

\item{algorithm}{"DP" (default) or "GA". See \code{\link[OncoBN]{fitCPN}}.}

\item{epsilon}{Penalty term. See \code{\link[OncoBN]{fitCPN}}.}

\item{silent}{Whether to show messages produced by OncoBN.}
}
\value{
A list with  a data.frame with From-To edges and associated thetas, the
thetas, the likelihood and epsilon of the model, and the predicted
probabilities of genotypes.
}
\description{
Run OncoBN
}
