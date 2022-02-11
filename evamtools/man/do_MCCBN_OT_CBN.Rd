\name{do_MCCBN_OT_CBN}
\alias{do_MCCBN_OT_CBN}
\title{Run MCCBN's OT-CBN.}
\usage{
do_MCCBN_OT_CBN(x)
}
\arguments{

\item{x}{Cross secitonal data. Matrix of genes (columns)
and individuals (rows)}
}
\value{
A list with the adjacency matrix, the lambdas, the parent set
and a data.frame with From-To edges and associated lambdas.
}
\description{
Run MCCBN's OT-CBN.

Using OT-CBN, as in example from \url{https://github.com/cbg-ethz/MC-CBN}, calling MC-CBN's \code{learn_network}, without the additional \code{sampling_times} argument (they are unknown in all of our cases). Function \code{learn_network} accepts more parameters, but the function is undocumented and the additional parameters are not used in the examples.
}
