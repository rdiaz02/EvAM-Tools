## Set up: add PATH_TO_HYPERTRAPS_REPO/src/python to your PATH
## Set up: add PATH_TO_HYPERTRAPS_REPO/bin/ to your PATH
## Set up: operate in a conda environment: Follow installation in README.md

do_HyperTraPS <- function(data, tmp_folder="", runs=1000, bi=50000, r=100, seed=42){
  ## Create tmp folder with random nameS
  dateTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  random_letters <- paste(c("_", tmp_folder, "_", LETTERS[floor(runif(4, min=1, max=26))]), sep="", collapse="")
  tmp_folder <- paste(c("/tmp/", dateTime, random_letters), collapse="")

  dir.create(tmp_folder)
  orig_folder <- getwd()
  setwd(tmp_folder)
  print(tmp_folder)

  ## Running HyperTraps
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
    # system("RUN_MCMC_SAMPLER -f transitions.txt -M second-order -N 1000 -r 20 -n zero -p 20 -k mcmc-apm -s 0.025 -b 50000 -i 100 -t 4 -q 0")
    system(sprintf("RUN_MCMC_SAMPLER -f transitions.txt -M second-order -N %s -b %s -r %s -S %s", runs, bi, r, seed))
  )["elapsed"]
  cat("Elapsed time for sampling posterior ", time_posterior, "\n")

  ### Run walkers
  print("Generating Paths")
  system(sprintf("RUN_PW -w match-data -b %s -R 100", bi))

  system(sprintf("plot_mc_stats.py -b %s", bi))
  setwd(orig_folder)
  # Cleaning 
  return(time_posterior)
}


