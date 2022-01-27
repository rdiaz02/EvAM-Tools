#!/usr/local/bin/Rscript
library(devtools)
load_all()


library("optparse")
 
option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = NULL, 
              help = "dataset file name", metavar = "character")
)
 
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

how2run <- function(opt){
    if(opt$file == "NULL"){ 
        runShiny()
        print("I will run Shiny")
    } else { ## Running command mode
        print("Running in command mode")
        tryCatch({
            csd <- read.csv(paste("/app/outside/", opt$file, sep = ""))
        }, error = function(e){
            stop("File could not be read. Place your file in the shared volume.")
        })

        if(!all(unique(unlist(csd)) %in% c(0, 1))){
            stop("Binary data should only contain 0 and 1")
        }

        out <- evam(csd, do_MCCBN = TRUE)
        out_with_simulations <- sample_all_CPMs(out, 100, 5)
        saveRDS(out_with_simulations, "/app/outside/cpm_out_with_simulations.rds")
    }
}

how2run(opt)

