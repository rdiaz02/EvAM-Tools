\name{do_MCCBN_OT_CBN}
\alias{do_MCCBN_OT_CBN}
\title{Run MCCBN's OT-CBN.}
\usage{
do_MCCBN_OT_CBN(
  x
)
}
\arguments{
%% \item{max.iter.asa}{Int. Number of steps to run. Default: 100000}

%% \item{ncores}{Int. Number of threads to use}

%% \item{L}{Int. Number of samples to be drawn from the proposal in the E-step}

%% \item{tmp_dir}{Directory name where the oput is located.}

%% \item{addname}{String to append to the temporary folder name. Default NULL}

%% \item{silent}{Whether to run show message showing the folder name where MCCBN is run}

\item{x}{Cross secitonal data. Matrix of genes (columns)
and individuals (rows)}
}
\value{
A list with the adjacency matrix, the lambdas, the parent set
and a data.frame with From-To edges and associated lambdas.
}
\description{
Run MCCBN's OT-CBN.

Using OT-CBN, as in example from https://github.com/cbg-ethz/MC-CBN , 
calling MC-CBN's \code{learn_network}, without the additional
\code{sampling_times} argument (they are unknown in all of our
cases). Function  \code{learn_network} accepts more parameters,
but the function is undocumented and the additional parameters are not used
in the examples.
}
