## Set up: add PATH_TO_HYPERTRAPS_REPO/src/python to your PATH
## Set up: add PATH_TO_HYPERTRAPS_REPO/bin/ to your PATH
## Set up: operate in a conda environment: Follow installation in README.md

do_HyperTraPS <- function(data){
  ## Create tmp folder with random nameS
  dateTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  random_letters <- paste(c("_", LETTERS[floor(runif(4, min=1, max=26))]), sep="", collapse="")
  tmp_folder <- paste(c("/tmp/", dateTime, random_letters), collapse="")

  dir.create(tmp_folder)
  orig_folder <- getwd()
  setwd(tmp_folder)

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
    system("RUN_MCMC_SAMPLER -f transitions.txt -M second-order")
  )["elapsed"]
  print("Elapsed time for sampling posterior ", time_posterior)

  ### Run walkers
  print("Generating Paths")
  system("RUN_PW -w match-data -b 50000 -R 100")

  setwd(orig_folder)
  # Cleaning 
  return(tmp_folder)
}

N <- 100
na <- N
nc <- N + round( 10 * runif(1))
nab <- 1.6 * N + round( 10 * runif(1))
ncd <- 1.5 * N + round( 10 * runif(1))
n00 <- round( 10 * runif(1))
dB <- matrix(
  c(
    rep(c(1, 0, 0, 0), na) 
    , rep(c(0, 0, 1, 0), nc)
    , rep(c(1, 1, 0, 0), nab)
    , rep(c(0, 0, 1, 1), ncd)        
    , rep(c(0, 0, 0, 0), n00)
  ), ncol = 4, byrow = TRUE
)
colnames(dB) <- LETTERS[1:4]

do_HyperTraPS(dB)
