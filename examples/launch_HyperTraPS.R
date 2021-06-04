source("HyperTraPS_data.R")
pwd0 <- getwd()
setwd("../code_from_what_genotype_next/")
source("hypertraps-process.R")
setwd(pwd0)
rm(pwd0)

sample_freqs <- function(data, save_file = NULL){
  if(!is.null(save_file)){
    jpeg(save_file, height = 350, units = "px")
  }
  par(mar = c(4, 2.5, 1, 0), las = 2)
  genes <- LETTERS[1 : ncol(data)]
  genes_freq <- table(apply(data, 1, function(x) {
      tmp_name <- paste0(genes[which(x == 1)], collapse = "")
      if(tmp_name == "") tmp_name <- "WT"
      return(tmp_name)
    }))
  barplot(genes_freq)
  dev.off()
}

bi = 20000
runs = 1000

for (dataset in names(all_examples)){
  print(sprintf("HyperTraPS_examples/HP_%s", dataset))
  tmp <- all_examples[dataset]
  do_HyperTraPS(tmp_data, 
    sprintf("HyperTraPS_examples/HP_%s", dataset), 
    runs = runs, bi = bi, dry_run = FALSE, plot = FALSE )
  sample_freqs(tmp_data, 
    sprintf("HyperTraPS_examples/HP_%s/freqs.jpg", dataset))
  ## generate markdown here
}



