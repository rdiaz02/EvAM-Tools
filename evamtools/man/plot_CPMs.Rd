\name{plot_CPMs}
\alias{plot_CPMs}
\title{Plot results from CPMs}
\usage{
plot_CPMs(
  cpm_output,
  samples = NULL,
  orientation = "horizontal",
  models = c("OT", "CBN", "OncoBN", "MCCBN", "MHN", "HESBCN"),
  plot_type = "trans_mat",
  fixed_vertex_size = FALSE,
  top_paths = NULL
)
}
\arguments{
\item{cpm_output}{Output from the cpm}

\item{samples}{Output from a call to \link{sample_CPMs}. Necessary if
  you request plot type \code{transitions}.}

\item{orientation}{String. If it is not "vertical" it will be displayed with an horizontal layout. Optional.}

\item{models}{Models for which you want output.}

\item{plot_type}{One of:
  \itemize{
    \item{trans_mat:}{Transition matrix between genotypes (see
  supporting information for OT and OncoBN). This plots the object of
  name \code{trans_mat} from the output of \code{\link{evam}}.}
    \item{trans_rate_mat:}{Transition rate matrix between genotypes;
  unavailable for OT and OncoBN. This plots the object of
  name \code{trans_rate_mat} from the output of \code{\link{evam}}.}
    \item{obs_genotype_transitions:}{Observed transitions during the
  simulation of the sampling process. This plots the object called  \code{obs_genotype_transitions} from the output of \code{\link{sample_CPMs}}.}
  }
}

\item{fixed_vertex_size}{Boolean. If TRUE, all nodes with have the same
  size; otherwise, scale them proportional to frequencies or observed data.}

\item{top_paths}{Number of most relevant paths to plot. Default NULL 
will plot all paths}
}

\description{
By default it creates a top row with the DAG of the CPM 
or de transtionRateMatrix for MHN
The bottom row has a custom plot for trans_mat, trans_rate_mat, or obs_genotype_transitions.
}

\note{The color and design of figures in the bottom row, depicting
  transition matrices, transition rate matrices, and observed genotype
  transitions are heavily inspired by (a blatant copy of) some of the
  representations in Greenbury et al., 2020.}

\references{
  Greenbury, S. F., Barahona, M., & Johnston, I. G. (2020). HyperTraPS:
  Inferring Probabilistic Patterns of Trait Acquisition in Evolutionary
  and Disease Progression Pathways. Cell Systems, 10(1),
  39–51–10. http://dx.doi.org/10.1016/j.cels.2019.10.009
}
\examples{
dB_c1 <- matrix(
 c(
     rep(c(1, 0, 0, 0, 0), 30) #A
   , rep(c(0, 0, 1, 0, 0), 30) #C
   , rep(c(1, 1, 0, 0, 0), 20) #AB
   , rep(c(0, 0, 1, 1, 0), 20) #CD
   , rep(c(1, 1, 1, 0, 0), 10) #ABC
   , rep(c(1, 0, 1, 1, 0), 10) #ACD
   , rep(c(1, 1, 0, 0, 1), 10) #ABE
   , rep(c(0, 0, 1, 1, 1), 10) #CDE
   , rep(c(1, 1, 1, 0, 1), 10) #ABCE
   , rep(c(1, 0, 1, 1, 1), 10) #ACDE
   , rep(c(1, 1, 1, 1, 0), 5) # ABCD
   , rep(c(0, 0, 0, 0, 0), 1) # WT
 ), ncol = 5, byrow = TRUE
)
colnames(dB_c1) <- LETTERS[1:5]

out <- evam(dB_c1,
            methods = c("CBN", "OT", "HESBCN", "MHN", "OncoBN", "MCCBN"))

plot_CPMs(out, plot_type = "trans_mat")

plot_CPMs(out, plot_type = "trans_rate_mat")

out_samp <- sample_CPMs(out, 1000)

plot_CPMs(out, out_samp, plot_type = "obs_genotype_transitions")

}
