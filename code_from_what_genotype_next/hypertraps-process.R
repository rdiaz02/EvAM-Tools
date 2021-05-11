## Set up: add PATH_TO_HYPERTRAPS_REPO/src/python to your PATH
## Set up: add PATH_TO_HYPERTRAPS_REPO/bin/ to your PATH
## Set up: operate in a conda environment: Follow installation in README.md

do_HyperTraPS <- function(data, tmp_folder="", runs=1000, bi=50000, r=100, seed=42, conda_env_name="HyperTraPS"){
  ## Create tmp folder with random nameS
  dateTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  random_letters <- paste(c("_", tmp_folder, "_", LETTERS[floor(runif(4, min=1, max=26))]), collapse="")

  tmp_folder_name <- paste(c(dateTime, random_letters), collapse="")
  tmp_folder <- file.path("/tmp", conda_env_name, tmp_folder_name), collapse="")

  dir.create(tmp_folder, recursive=TRUE)
  orig_folder <- getwd()
  setwd(tmp_folder)
  print(tmp_folder)

  ## Activating conda env
  system(sprintf("conda activate %s", conda_env_name))

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
    system(sprintf("RUN_MCMC_SAMPLER -f transitions.txt -M second-order -N %1.f -b %1.f -r %1.f -S %1.f", runs, bi, r, seed))
  )["elapsed"]
  cat("Elapsed time for sampling posterior ", time_posterior, "\n")

  ### Run walkers
  print("Generating Paths")
  system(sprintf("RUN_PW -w match-data -b %1.f -R 100", bi))

  system(sprintf("plot_mc_stats.py -b %1.f", bi))
  setwd(orig_folder)
  # Cleaning 
  return(time_posterior)
}


