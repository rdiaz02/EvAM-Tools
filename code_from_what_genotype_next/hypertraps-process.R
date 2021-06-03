## Set up: add PATH_TO_HYPERTRAPS_REPO/src/python to your PATH
## Set up: add PATH_TO_HYPERTRAPS_REPO/bin/ to your PATH
## Set up: operate in a conda environment: Follow installation in README.md
library(reticulate) ## To activate the conda environment

do_HyperTraPS <- function(data, tmp_folder="", 
  runs=1000, bi=50000, r=100, 
  seed=0, conda_env_name="HyperTraPS", plot = T, dry_run = F){
  
  orig_folder <- getwd()
  
  if(!dry_run){
    ## Activating conda env
    warning("You should run the script within a conda enviroment with a working installation of HyperTraPS")
    use_condaenv(conda_env_name) ## This does not seem to work

    ## Running HyperTraps
    ### Handling data 
    # browser()getw
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
    if(tmp_folder == "") {
      dateTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
      random_letters <- paste(c("_", 
        LETTERS[floor(runif(4, min = 1, max = 26))]), collapse = "")
      tmp_folder_name <- paste(c(dateTime, random_letters), collapse = "")
      tmp_folder <- file.path("/tmp", conda_env_name, tmp_folder_name)
    }

    dir.create(tmp_folder, recursive = TRUE)
    setwd(tmp_folder)
    print(tmp_folder)

    ## Writin data to tmp folder

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
    system(sprintf("RUN_PW -w match-data -b %1.f -R 100", bi))
    system(sprintf("RUN_PW -w zero-one -b %1.f -R 100", bi))
  } else { setwd(tmp_folder) }

  ### Plots
  system(sprintf("plot_mc_stats.py -b %1.f", bi))
  system(sprintf("custom_plot_hypercube_graph.py -f 'forwards_list-pord-match-data.csv' -outfile_graph 'forwards-hypercube-graph-mach-data-g0' -transition_data 'transitions.txt' -labels 'labels.csv' -label_type 'None' -labels_fontsize 4 -layout_type 'spring' -aspect 0.9 -width 3.4 -out_type 'png'"))
    
  system(sprintf("plot_feature_graph.py -f 'forwards.txt' -layout_type 'circular' -prob_type 'joint' -data_type 'match-data' -width 4 -fontsize 10 -any_time 0 -node_size 100 -connection_style 'arc3,rad=-0.3' -outfile_type 'png'"))

  system(sprintf("plot_ordering_histogram.py -f 'forwards_list-pord-zero-one.csv' -f2 'forwards_list-pord-match-data.csv' -transition_data 'transitions.txt' -labels 'labels.csv' -fontsize 6 -xevery 1 -aspect 0.9 -verbose 'no' -outfile 'forwards-ws1-ws2' -out_type 'png'"))
  
  system(sprintf("hypertraps_get_output.py -f 'forwards.txt' -prob_type 'conditional'"))

  if(plot){
    library(imager)
    stats_names <- unlist(strsplit(list.files(pattern = "^stats.*png"), " "))
    files <- list(
      stats_names[1], 
      "forwards-hypercube-graph-mach-data-g0.png"
      , "forwards-ws1-ws2.png")
    for(f in files){
      im <- load.image(f)
      plot(im)
    }
  }
  ### Extract data from the model:
  ## Posets and transitions numbers
  ## This shoulb be easy because all this data is plotted already
  ## So I should generate a python script that mimics the functions
  ## that extract the data in order to save it

  ## plot_feature_graph.py includes MakeFeatureAdjacencyJoint
  ## 

  # Cleaning 
  setwd(orig_folder)
  return(NULL)
}


