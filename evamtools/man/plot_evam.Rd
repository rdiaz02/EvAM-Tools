\name{plot_evam}
\alias{plot_evam}
\alias{plot_CPMs}
\alias{plot_BML_all}
\alias{plot_HT_as_evam_trans_mat}
\title{Plot results from EvAMs (CPMs)}
\description{Plots fitted EvAMs (CPMs), both fitted model and custom
  plots for transition rates and transition probabilities.}

\usage{
plot_evam(
  cpm_output,
  samples = NULL,
  orientation = "horizontal",
  methods = NULL,
  plot_type = "trans_mat",
  label_type = "genotype",
  fixed_vertex_size = FALSE,
  top_paths = NULL
)

plot_CPMs(
  cpm_output,
  samples = NULL,
  orientation = "horizontal",
  methods = NULL,
  plot_type = "trans_mat",
  label_type = "genotype",
  fixed_vertex_size = FALSE,
  top_paths = NULL
)

## This is a temporary kludge for BML
plot_BML_all(x)

## What about HyperTraPS? Use the plotting functions provided
## by the hypertrapsct package. See Note and examples.
## But if you want to see the transition matrix as for the rest of evam you can use

plot_HT_as_evam_trans_mat(
  cpm_output,
  label_type = "genotype",
  fixed_vertex_size = FALSE,
  top_paths = NULL
)


}
\arguments{
  \item{cpm_output}{Output from the cpm}

  \item{samples}{Output from a call to \link{sample_evam}. Necessary if
  you request plot type \code{transitions}.}

  \item{orientation}{String. If it is not "vertical" it will be displayed with an horizontal layout. Optional.}

  \item{methods}{Vector of strings with the names of methods that we want to plot. If NULL, all
  methods with output in cpm_output. The list of available methods is
  OT, OncoBN, CBN, MCCBN, MHN, HESBCN.}

  \item{plot_type}{One of:
  \describe{
    \item{trans_mat: }{Transition matrix between genotypes (see
  supporting information for OT and OncoBN). This plots the object of
  name \code{trans_mat} from the output of \code{\link{evam}}.}
    \item{trans_rate_mat: }{Transition rate matrix between genotypes;
  unavailable for OT and OncoBN. This plots the object of
  name \code{trans_rate_mat} from the output of \code{\link{evam}}.}
    \item{obs_genotype_transitions: }{Observed transitions during the
  simulation of the sampling process. This plots the object called  \code{obs_genotype_transitions} from the output of \code{\link{sample_evam}}.}
  }
}

\item{label_type}{Type of label to show. One of:
  \describe{
    \item{genotype: }{Displays all genes mutated}
    \item{acquisition: }{Only displays the last gente mutated}
  }
}

\item{fixed_vertex_size}{Boolean. If TRUE, all nodes with have the same
  size; otherwise, scale them proportional to frequencies of observed data.}

\item{top_paths}{Number of most relevant paths to plot. Default NULL
will plot all paths. See below, Description, for details about
relevance. With many genes, and particularly for MHN (or other methods,
when there are no restrictions in the order of accumulation of
mutations), using NULL can lead to a very long time to plot.}

\item{x}{An object that contains the output from runs that used BML.}
}

\value{
By default this function creates a top row with the DAG of the CPM
or the log-Theta matrix for MHN.
The bottom row has a custom plot for the transition matrix, or the transition
rate matrix, or the observed genotype transitions.

 In the bottom row plots, unless \code{fixed_vertex_size = TRUE}, the size of genotype nodes is proportional to
the observed frequency of genotypes for all plots except for
\code{plot_type = 'obs_genotype_transitions'}, where it is
proportional to the genotype frequency as obtained by the sampling from
the predictions of each method; if a node (genotype) has no observations,
its size is of fixed size. Thus, for all plots except \code{plot_type = 'obs_genotype_transitions'} the size of the genotype nodes is
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


For BML, we are using, for the moment, a temporary kludge. It returns
figs as in Fig. 2 and Fig. 3a and 3b in Misra et al., 2014.



}



\note{The color and design of figures in the bottom row, depicting
  transition matrices, transition rate matrices, and observed genotype
  transitions are heavily inspired by (a blatant copy of) some of the
  representations in Greenbury et al., 2020.

  It is easy to get the plots to display poorly (overlapping names in
  nodes, overlapping labels between plots, etc) if you use long gene
  names. For best results, try to use short gene names.


The criteria to decide which are the most relevant paths with the
top_paths option is the following:
1) Get all the leaves from the graph.
2) Calculate all the paths leading from "WT" to all leaves.
3) Select n paths with highest cumulative weighted sum. The weights used
depend on the type of plot (\code{plot_type}): for trans_mat it will be
log probabilities (so the most relevant paths are the most likely
paths), for trans_rate_mat it will be rates, and for obs_genotype_transitions it will be raw counts.


Plots can be more readable with a combination of top_paths and
label_type. If \code{label_type = 'acquisition'} node labels will dissapear
and edge label will be shown instead. They will display the information of
the last gene mutated.

What about HyperTraPS? Use the plots provided by the
\code{\link[hypertrapsct:hypertrapsct]{hypertrapsct}}
package (see also \url{https://github.com/StochasticBiology/hypertraps-ct/tree/bioconductor}) passing it the component \code{$HyperTraPS_all_output} of the evam object. For example,
\code{\link[hypertrapsct:plotHypercube.influencegraph]{plotHypercube.influencegraph}}
or
\code{\link[hypertrapsct:plotHypercube.summary]{plotHypercube.summary}},
etc. See an example below.


\strong{plot_CPMs has been deprecated}. Use plot_evam.

}

\references{
  Greenbury, S. F., Barahona, M., & Johnston, I. G. (2020). HyperTraPS:
  Inferring Probabilistic Patterns of Trait Acquisition in Evolutionary
  and Disease Progression Pathways. Cell Systems, 10(1),
  39–51–10. \url{http://dx.doi.org/10.1016/j.cels.2019.10.009}

  Misra, N., Szczurek, E., & Vingron, M. (2014). Inferring the paths of somatic evolution in cancer. Bioinformatics (Oxford, England), 30(17), 2456–2463. \url{https://doi.org/10.1093/bioinformatics/btu319}

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

## Use MCCBN only if installed
MCCBN_INSTALLED <- requireNamespace("mccbn", quietly = TRUE)
methods <- c("CBN", "OT", "OncoBN", "MHN", "HESBCN")
if (MCCBN_INSTALLED) {
    methods <- c(methods, "MCCBN")
}

out <- evam(dB_c1,
            methods = methods)

plot_evam(out, plot_type = "trans_mat")

plot_evam(out, plot_type = "trans_rate_mat")

plot_evam(out, plot_type = "trans_rate_mat", top_paths=2)
plot_evam(out, plot_type = "trans_rate_mat", top_paths=2
  , label_type ="acquisition")

out_samp <- sample_evam(out, 1000, output = c("sampled_genotype_counts", "obs_genotype_transitions"))

plot_evam(out, out_samp, plot_type = "obs_genotype_transitions")

## Only showing new gene mutated respect with its parent
plot_evam(out, out_samp, plot_type = "obs_genotype_transitions",
  label_type = "acquisition")

plot_evam(out, out_samp, plot_type = "obs_genotype_transitions",
  label_type = "acquisition", top_paths = 3)


## Examples with mixed AND and OR and AND and XOR for HESBCN
data("ex_mixed_and_or_xor")

out_AND_OR_XOR <- evam(ex_mixed_and_or_xor,
                       methods = c("OT", "HESBCN", "MHN", "OncoBN"),
                       hesbcn_opts = list(seed = 26))

plot_evam(out_AND_OR_XOR,plot_type = "trans_mat",
          top_paths = 3)


## Asking for a method not in the output will give a warning
plot_evam(out_AND_OR_XOR, plot_type = "trans_mat",
          methods = c("OT", "OncoBN"),
          top_paths = 4)


## Only two methods, but one not fitted
plot_evam(out_AND_OR_XOR, methods = c("CBN", "HESBCN"),
           plot_type = "trans_mat")

## Only one method
plot_evam(out_AND_OR_XOR, methods = c("MHN"),
          plot_type = "trans_mat",
          top_paths = 5)

plot_evam(out_AND_OR_XOR, plot_type = "trans_mat", top_paths = 3)


plot_evam(out_AND_OR_XOR, methods = c("MHN", "HESBCN"),
          plot_type = "trans_mat", label_type="acquisition"
          , top_paths=3)

plot_evam(out_AND_OR_XOR, methods = c("MHN", "HESBCN"),
          plot_type = "trans_mat", label_type="genotype"
          , top_paths=3)


## HyperTraPS: using hypertrapsct package function
## and showing the transition matris like the rest of evam
out_h_1 <- evam(dB_c1,
                methods = c("HyperTraPS"),
                hyper_traps_opts = list(length = 2,
                                        model = -1,
                                        seed = -1))
hypertrapsct::plotHypercube.summary(out_h_1$HyperTraPS_all_output)
plot_HT_as_evam_trans_mat(out_h_1, top_paths = 4)

## Three BML examples
out_B1 <- evam(dB_c1,
               methods = c("BML"),
             bml_opts = list(ntree = 3, rep = 7))

out_B_nb <- evam(dB_c1,
                 methods = c("BML"),
                bml_opts = list(ntree = 3, rep = 0))

out_B_OT <- evam(dB_c1,
                 methods = c("BML", "OT"),
              bml_opts = list(ntree = 3, rep = 7))

plot_BML_all(out_B1)
plot_BML_all(out_B_nb)
plot_BML_all(out_B_OT)


}
