\name{plot_CPMs}
\alias{plot_CPMs}
\title{Plot results from CPMs}
\usage{
plot_CPMs(
  cpm_output,
  samples = NULL,
  orientation = "horizontal",
  methods = c("OT", "OncoBN", "CBN", "MCCBN", "HESBCN", "MHN"),
  plot_type = "trans_mat",
  label_type = "genotype",
  fixed_vertex_size = FALSE,
  top_paths = NULL
)
}
\arguments{
\item{cpm_output}{Output from the cpm}

\item{samples}{Output from a call to \link{sample_CPMs}. Necessary if
  you request plot type \code{transitions}.}

\item{orientation}{String. If it is not "vertical" it will be displayed with an horizontal layout. Optional.}

\item{methods}{Methods for which you want output.}

\item{plot_type}{One of:
  \itemize{
    \item{trans_mat: }{Transition matrix between genotypes (see
  supporting information for OT and OncoBN). This plots the object of
  name \code{trans_mat} from the output of \code{\link{evam}}.}
    \item{trans_rate_mat: }{Transition rate matrix between genotypes;
  unavailable for OT and OncoBN. This plots the object of
  name \code{trans_rate_mat} from the output of \code{\link{evam}}.}
    \item{obs_genotype_transitions: }{Observed transitions during the
  simulation of the sampling process. This plots the object called  \code{obs_genotype_transitions} from the output of \code{\link{sample_CPMs}}.}
  }
}

\item{label_type}{Type of label to show. One of:
  \itemize{
    \item{genotype: }{Displays all genes mutated}
    \item{acquisition: }{Only displays the last gente mutated}
  }
}

\item{fixed_vertex_size}{Boolean. If TRUE, all nodes with have the same
  size; otherwise, scale them proportional to frequencies of observed data.}

\item{top_paths}{Number of most relevant paths to plot. Default NULL 
will plot all paths}
}

\description{
By default it creates a top row with the DAG of the CPM 
or the transtionRateMatrix for MHN
The bottom row has a custom plot for trans_mat, trans_rate_mat, or
obs_genotype_transitions.

 In the bottom row plots, unless \code{fixed_vertex_size = TRUE}, the size of genotype nodes is proportional to
the observed frequency of genotypes for all plots except for
\code{plot_type = 'obs_genotype_transitions'}, where it is
proportional to the genotype frequency as obtained by the sampling from
the predictions of each method; if a node (genotype) has no observations,
its size is of a fixed size. Thus, for all plots except \code{plot_type = 'obs_genotype_transitions'} the size of the genotype nodes is
the same among methods, but the size of the genotype nodes can differ
between methods for \code{plot_type = 'obs_genotype_transitions'}. 


In the top plots, in the DAGs, when a node has two or more incoming
edges, color depends on the type of relationship. (With a single
incoming edge, there is no difference in model behavior with type of
edge and, for consistency with OT and CBN, nodes with a Single parent
have edges colored the same way as nodes with two or more parents and
AND relationship).

In the bottom row plots, non-observed genotypes are shown in light
green, to differentiate them from the observed genotypes (shown in orange).

The criteria to decide which are the most relevant paths with the top_paths option
is the following:
1) get all the leaves from the graph.
2) calculate all the paths leading from "WT" to all leaves.
3) select n paths with highest cumulate weight. The nature of the weights
depends on the plot_type: for trans_mat it will be probabilities, 
for trans_rate_mat it will be rates, and for obs_genotype_transitions it will be
raw counts.  

Plots can be more readable with a combination of top_paths and
label_type. If \code{label_type = 'acquisition'} node labels will dissapear
and edge label will be shown instead. They will display the information of
the last gene mutated.

}



\note{The color and design of figures in the bottom row, depicting
  transition matrices, transition rate matrices, and observed genotype
  transitions are heavily inspired by (a blatant copy of) some of the
  representations in Greenbury et al., 2020.

  It is easy to get the plots to display poorly (overlapping names in
  nodes, overlapping labels between plots, etc) if you use long gene
  names. For best results, try to use short gene names.

}

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

plot_CPMs(out, plot_type = "trans_rate_mat", top_paths=2)
plot_CPMs(out, plot_type = "trans_rate_mat", top_paths=2
  , label_type ="acquisition")

out_samp <- sample_CPMs(out, 1000, output = c("sampled_genotype_freqs", "obs_genotype_transitions"))

plot_CPMs(out, out_samp, plot_type = "obs_genotype_transitions")

## Only showing new gene mutated respect with its parent
plot_CPMs(out, out_samp, plot_type = "obs_genotype_transitions", 
  label_type = "acquisition")

plot_CPMs(out, out_samp, plot_type = "obs_genotype_transitions", 
  label_type = "acquisition", top_paths = 3)


## Examples with mixed AND and OR and AND and XOR for HESBCN
data("ex_mixed_and_or_xor")

out_AND_OR_XOR <- evam(ex_mixed_and_or_xor,
                       methods = c("OT", "HESBCN", "MHN", "OncoBN"),
                       hesbcn_opts = list(seed = 26,
                                          silent = FALSE))

plot_CPMs(out_AND_OR_XOR,plot_type = "trans_mat")

         

## Only two methods, but one not fitted
plot_CPMs(out_AND_OR_XOR, methods = c("CBN", "HESBCN"),
           plot_type = "trans_mat")

## Only one method
plot_CPMs(out_AND_OR_XOR, methods = c("MHN"),
          plot_type = "trans_mat")

plot_CPMs(out_AND_OR_XOR, plot_type = "trans_mat", top_paths = 3)



plot_CPMs(out_AND_OR_XOR, methods = c("MHN", "HESBCN"),
          plot_type = "trans_mat")

plot_CPMs(out_AND_OR_XOR, methods = c("MHN", "HESBCN"),
          plot_type = "trans_mat", label_type="acquisition"
          , top_paths=3)
plot_CPMs(out_AND_OR_XOR, methods = c("MHN", "HESBCN"),
          plot_type = "trans_mat", label_type="genotype"
          , top_paths=3)
}
