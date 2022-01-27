\name{evam}
\alias{evam}
\title{Runs the CPMs}
\usage{
evam(x, cores_cbn = 1, methods = c("CBN", "OT", "HESBCN", "MHN", "DBN"),
     max.cols = 15)
}
\arguments{
\item{x}{cross sectional data}

\item{cores_cbn}{How many cores to use for CBN}

\item{methods}{Methods to use. For now, the only thing that matters is
  whether the vector includes the strings "MCCBN" and "DBN". The default
  is not to use MCCBN, so c("CBN", "OT", "HESBCN", "MHN", "DBN")}

\item{max.cols}{Maximum number of columns to use in the analysis. If x has >
max.cols, selected columns are those with the largest number of events.}
}
\description{
Executes all CPMS given a cross sectional data set
}
