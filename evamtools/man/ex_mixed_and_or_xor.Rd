\name{ex_mixed_and_or_xor}
\alias{ex_mixed_and_or_xor}
\docType{data}
\title{
  Small example data set that shows, for HESBCN, both AND and OR, or AND
  and XOR, in some runs.
}

\description{
  Synthetic data set used for testing and plotting. HESBCN, with some
  seeds, will infer  both AND and OR, or AND
  and XOR.
  

}
\usage{
data("ex_mixed_and_or_xor")
}

\format{
A data frame with \dQuote{genes} in columns and \dQuote{patients} in rows. 
}



\examples{
data("ex_mixed_and_or_xor")

out_AND_OR <- evam(ex_mixed_and_or_xor,
               methods = c("OT", "HESBCN", "MHN", "OncoBN"),
               hesbcn_opts = list(seed = 26))

out_AND_XOR <- evam(ex_mixed_and_or_xor,
               methods = c("OT", "HESBCN", "MHN", "OncoBN"),
               hesbcn_opts = list(seed = 31))

plot_CPMs(out_AND_OR, plot_type = "trans_mat")

plot_CPMs(out_AND_XOR, plot_type = "trans_mat")

}
\keyword{datasets}
