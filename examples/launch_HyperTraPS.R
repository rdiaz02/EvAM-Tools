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
  genes_names <- names(genes_freq)
  try(genes_names <- genes_names[-which(genes_names == "WT")])
  ordered_names <- c("WT", genes_names[order(sapply(genes_names, nchar))])
  ordered_genotypes <- genes_freq[ordered_names]
  barplot(ordered_genotypes)
  dev.off()
}

bi = 20000
runs = 1000
dry_run = FALSE
for (dataset in names(all_examples)){
  print(sprintf("HyperTraPS_examples/HP_%s", dataset))
  tmp_data <- all_examples[dataset]
  do_HyperTraPS(tmp_data, 
    sprintf("HyperTraPS_examples/HP_%s", dataset), 
    runs = runs, bi = bi, dry_run = dry_run, plot = FALSE )
  sample_freqs(tmp_data, 
    sprintf("HyperTraPS_examples/HP_%s/freqs.jpg", dataset))
  ## generate markdown here
}



