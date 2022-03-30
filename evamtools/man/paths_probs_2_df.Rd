\name{paths_probs_2_df}
\alias{paths_probs_2_df}
\title{Utility to convert list with paths and their log-probabilities to
  data frame with paths as strings and probabilities.}
\usage{
paths_probs_2_df(x, order = c("prob", "path"))
}
\arguments{
  \item{x}{A list with two components, \code{paths} and
    \code{weights}. \code{paths} is a list of \code{igraph.vs} objects,
    and \code{weights} is a vector of log-probabilities. This is the
    type of return object from the internal function
    \code{trans_mat_2_paths_probs} and the type of object in all of the
    \code{*_paths_max} returned by \code{\link{evam}}.}

\item{order}{How should the output be ordered. If \code{prob}, order by
  decreasing probability (solving ties by increasing alphanumeric order
  of the path string). If \code{path}, order increasingly according to path.} 

}


\value{
  A data frame with two columns, \code{Path} and \code{Prob}.
}


\description{
  Utility to convert list with paths and their log-probabilities to
  data frame with paths as strings and probabilities and return output
  ordered according to two possible criteria.
}

\seealso{
  \code{\link{evam}}
}

\examples{

data(every_which_way_data)
d11 <- every_which_way_data[[16]][1:40, 2:5]

## For more interesting output, add HESBCN to the methods
out11 <- evam(d11,
              methods = c("MHN", "OncoBN"),
              paths_max = TRUE)

paths_probs_2_df(out11$MHN_paths_max, order = "prob")
paths_probs_2_df(out11$MHN_paths_max, order = "path")

}
