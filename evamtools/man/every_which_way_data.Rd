\name{every_which_way_data}
\alias{every_which_way_data}
\docType{data}
\title{
  Cancer data sets from Diaz-Uriarte and Vasallo, 2019; also used in
  Diaz-Colunga and Diaz-Uriarte, 2021.
}

\description{

    

}
\usage{data("every_which_way_data")}

\format{
  A list of length 22. Each element of the list is a data set, with
  subjects in rows and genes/probes in columns.
}

\references {


}


\examples{
\dontrun{
data(every_which_way_data)
lapply(every_which_way_data, colnames)
lapply(every_which_way_data, dim)
}
}
\keyword{datasets}
