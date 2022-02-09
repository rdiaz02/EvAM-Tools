## Copyright 2021, 2022 Pablo Herrera Nieto

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Set up: add PATH_TO_HYPERTRAPS_REPO/src/python to your PATH
## Set up: add PATH_TO_HYPERTRAPS_REPO/bin/ to your PATH
## Set up: operate in a conda environment: Follow installation in README.md
#library(imager)


if(FALSE) {

#' Runs HyperTraPS model
#' This function runs command line tools to call to HyperTraps.
#' Then processes the output
#' 
#' HOW SIMULATIONS ARE PERFORMED ()
#' Greenbury, Sam F. Barahona, Mauricio Johnston, Iain G. 2020
#' Simulated Walks to Illustrate Order of Acquisition (pages e3-e4)
#' The inference process above yields inferred posterior distributions on the hypercubic edge weightsW.
#' We can query these posteriors in a number of ways to gain descriptive and predictive information about the mechanisms generating observed states. 
#' First, we pro- duce a parsimonious and intuitive representation of the dynamic pathways supported by the inferred posteriors. 
#' Here, we simulate an ensemble of random walkers generating complete trajectories on hypercubes with sets of transition probabilities sampled from the inferred posterior. 
#' This ensemble reflects the likely dynamic pathways supported by the dynamic transition model after parameter- ization. 
#' We simulate an ensemble of random walks in two ways:
#'  Walk Simulation 1 (WS1), with walkers that run from f0gL to f1gL where a feature is acquired at every time step, and 
#' Walk Simulation 2 (WS2) which only simulates trajectories corresponding to tran- sitions observed in the dataset. 
#' In each case, we record every transition between states allowing the construction of a weighted directed graph of all states and transitions encountered. From this graph, the frequency fij with which feature i is gained at step j.
#' Code for sampling is found in the HypeTraPS repository
#' HyperTraPS/src/{run_pw.cpp,updater_fast.cpp}
#' @param data Str data file to read
#' @param data Dataframe with all the cross sectional data
#' @param tmp_folder Str. Folder to store and read all the results 
#' @param runs Int. Number of posterior samples
#' @param bi Int. Burn in from MCMC sampler
#' @param r Int. HyperTraPS trajectories per transition
#' @param seed Int. Ramdom seed for the code to use.
#' @param prob_type "joint" or "conditional" probaility type for transitions
#' @param dry_run Boolean. Whether to run the plot command.
#' @return List$edges
#' @return List$feature_transitions
#' @examples 
#'\dontrun{
#'dB_OR <- matrix(
#'  c(
#'     rep(c(1, 0, 0, 0), 200) #A
#'     , rep(c(1, 0, 1, 0), 100) #AC
#'     , rep(c(1, 1, 0, 0), 100) #AB
#'     , rep(c(1, 1, 1, 0), 50) #ABC
#'     , rep(c(1, 1, 0, 1), 50) #ABD
#'     , rep(c(1, 0, 1, 1), 50) #ACD
#'     , rep(c(1, 1, 1, 1), 10) #ABCD
#'     , rep(c(0, 0, 0, 0), 10) #WT
#'   ), ncol = 4, byrow = TRUE
#' )
#' colnames(dB_OR) <- LETTERS[1:4]
#' 
#' do_HyperTraPS(tmp_data, 
#'       sprintf("./HP_%s", dataset), 
#'       runs = runs, bi = bi, dry_run = dry_run, show_plot = FALSE )
#' }

do_HyperTraPS <- function(data = NULL, tmp_folder = NULL, 
  runs = 1000, bi = 50000, r = 100, 
  seed = 0, prob_type = "joint", 
  # plot = TRUE, show_plot = TRUE, 
  dry_run = FALSE){
  
  warning("You should run the script within a conda enviroment with a working installation of HyperTraPS")
  orig_folder <- getwd()

  if (!(prob_type %in% c("joint", "conditional"))) {
    stop(sprintf("Invalid probability type %s.
      Supported types are: 1.- joint and 2.- conditional", prob_type))
  }
  
  if (!dry_run & !is.null(data)) {

    ## Running HyperTraps
    ### Handling data 
    if(typeof(data) == "character") {
      tryCatch(expr = {
        data <- read.table(data, sep = ",", header = T)
      },
      error = {
        function(e){
          stop(sprintf("Could not locate the file", data))
        }
      })
    }
    
    ## Create tmp folder with random name
    if(is.null(tmp_folder)) {
      dateTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
      random_letters <- paste(c("_", 
        LETTERS[floor(runif(4, min = 1, max = 26))]), collapse = "")
      tmp_folder_name <- paste(c(dateTime, random_letters), collapse = "")
      tmp_folder <- file.path("/tmp", tmp_folder_name)
    }
    dir.create(tmp_folder, recursive = TRUE) 

    setwd(tmp_folder)
    print(getwd())

    ## Writing data to tmp folder

    output_name <- "data.csv"
    write.csv(data, output_name)
    
    ### Transforming data
    ## TODO modify input_type
    print("Converting Data")
    system("convert_data_to_transitions.py -data data.csv -input_type cross-sectional")

    ### MCM sampler
    ## TODO check the meaning of all parameters and allow to set them
    print("Sampling posterior")
    time_posterior <- system.time(
      system(sprintf("RUN_MCMC_SAMPLER -f transitions.txt -M second-order -N %1.f -b %1.f -r %1.f -S %1.f", runs, bi, r, seed))
    )["elapsed"]
    cat("Elapsed time for sampling posterior ", time_posterior, "\n")

    ### Run walkers
    print("Generating Paths")
    system(sprintf("RUN_PW -w match-data -b %1.f", bi))
    # system(sprintf("RUN_PW -w zero-one -b %1.f -R 100", bi))

    ##Here I force to use conditional probabilities only
    system(sprintf("hypertraps_get_output.py -f 'forwards.txt' -prob_type %s", "conditional")) 
  } else { setwd(tmp_folder) }

  ### Plots
  # if (plot){
  #   system(sprintf("plot_mc_stats.py -b %1.f", bi))
  #   system(sprintf("custom_plot_hypercube_graph.py -f 'forwards_list-pord-match-data.csv' -outfile_graph 'forwards-hypercube-graph-mach-data-g0' -transition_data 'transitions.txt' -labels 'labels.csv' -label_type 'None' -labels_fontsize 4 -aspect 0.9 -width 3.4 -out_type 'png'"))
      
  #   system(sprintf("plot_feature_graph.py -f 'forwards.txt' -layout_type 'circular' -prob_type %s -data_type 'match-data' -width 4 -fontsize 10 -any_time 0 -node_size 100 -connection_style 'arc3,rad=-0.3' -outfile_type 'png'", prob_type))

  #   system(sprintf("plot_ordering_histogram.py -f 'forwards_list-pord-zero-one.csv' -f2 'forwards_list-pord-match-data.csv' -transition_data 'transitions.txt' -labels 'labels.csv' -fontsize 6 -xevery 1 -aspect 0.9 -verbose 'no' -outfile 'forwards-ws1-ws2' -out_type 'png'"))
    
  # }

  # if(show_plot){
  #   library(imager)
  #   stats_names <- unlist(strsplit(list.files(pattern = "^stats.*png"), " "))
  #   files <- list(
  #     stats_names[1], 
  #     "forwards-hypercube-graph-mach-data-g0.png"
  #     , "forwards-ws1-ws2.png")
  #   for(f in files){
  #     im <- load.image(f)
  #     plot(im)
  #   }
  # }

  ## Reading features 
  features <- read.csv("feature_transitions.csv", header = FALSE)
  genes <- c("Root", LETTERS[1: (ncol(features) - 1)])
  colnames(features) <- genes
  rownames(features) <- genes
  model <- features2model(features)

  ## Compute Transition Rate Matrix

  ## Compute Transition Count Matrix
  genotypes_transition_counts <- read.table("forwards_list-pord-edge-list-long-match-data.txt", header = TRUE)

  genotypes_transition_counts$from <- vapply(genotypes_transition_counts$from, int2str, character(1))
  genotypes_transition_counts$to <- vapply(genotypes_transition_counts$to, int2str, character(1))
  g <- graph_from_data_frame(genotypes_transition_counts)
  transition_count_matrix <- igraph::as_adjacency_matrix(g, attr="weight")
  states <- setdiff(colnames(transition_count_matrix), "WT")
  sorted_states <- c("WT", states[order(vapply(states, nchar, numeric(1)))])
  transition_count_matrix <- transition_count_matrix[sorted_states, sorted_states]

  ## Compute Transition Probability Matrix
  transition_probability_matrix <- transition_count_matrix / rowSums(transition_count_matrix)

  # Cleaning 
  setwd(orig_folder)
  return(list(
    feature_transitions = features,
    transition_count_matrix = transition_count_matrix,
    transition_probability_matrix = transition_probability_matrix,
    edges = model))
}


features2model <- function(data){  
  genes <- colnames(data)
  positions <- data.frame(which(data > 0, arr.ind = TRUE))
  parents <- genes[positions$row]
  childrens <- genes[positions$col]
  edges <- paste(parents, childrens, sep = "->")
  probabilities <- mapply(function(a,b) data[a,b], positions$row, positions$col)

  model <- data.frame(From = parents,
                      To = childrens,
                      Edge = edges,
                      Probabilities = probabilities,
                      stringsAsFactors = FALSE)

  return(model)
}


create_transition_probabi <- function() {}
}
