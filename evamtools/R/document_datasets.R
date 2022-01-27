#' Cross sectional data sets
#'
#' A list of cross sectional data set to be used as example inputs.
#'
#' @format A list with 3 items:
#' \describe{
#'   \item{csd}{csd, with 17 different scenarios. They represents different relationships: AND, OR, and XOR, but there are also more complex set ups, like sign epitasis and reciprocal sign epistasis. }
#'   \item{dag}{dag, It has a simple example of an AND relationships}
#'   \item{matrix}{matrix, Simple random example}
#' }
"examples_csd"


#' Output of CPMS
#'
#' A list containing the output of several CPMs for 10 cross sectional data sets.
#'
#' @format Each element cotains 27 elements. This includes different outputs for 
#' each of the following CPMs: OT, CBN, MHN, HESBCN. DBN and MCCBN are being considered.
#' Some of the specific ouputs are:
#' \describe{
#'   \item{csd_data}{csd_data, Dataframe with the cross sectional data that generated the results}
#'   \item{*_trans_rate_mat}{Transition rate matrix.}
#'   \item{*_f_graph}{For OT (and DBN?) what we call weighted fitness graph. These are not really probabilities, though they are obtained with an algorithm identical to the one used for the transition rate matrix.}
#'   \item{*_model}{Contains the edges of the gene relationship DAG}
#'   \item{*_trans_mat}{Transition matrix between genotypes.}
#'   \item{*_td_trans_mat}{Time discretized transition matrix.}
#'   \item{OT_genots_predicted}{Genotype frequency for OT model}
#' }
"cpm_output"

#' Output from the sampling process
#'
#' A list containing the output of several CPMs for 10 cross sectional data sets.
#'
#' @format Each element cotains 27 elements. This includes different outputs for 
#' each of the following CPMs: OT, CBN, MHN, HESBCN. DBN and MCCBN are being considered.
#' Some of the specific ouputs are:
#' \describe{
#'   \item{trajectory}{trajectory, Every nested list shows genotypes jumps seen in a single sample}
#'   \item{obs_events}{obs_events, List with the laste genotypes of the trajectory: the observed events}
#' }
"sampling_output"
