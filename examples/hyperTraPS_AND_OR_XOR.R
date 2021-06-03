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
      return(tmp_name))
    })
  browser()
  barplot(genes_freq)
  dev.off()
}


dB_AND <- matrix(
  c(
    rep(c(1, 0, 0, 0), 100)
    , rep(c(1, 1, 0, 0), 100)
    , rep(c(1, 1, 1, 0), 100)
    , rep(c(1, 1, 1, 1), 100)
    , rep(c(0, 0, 0, 0), 10)
  ), ncol = 4, byrow = TRUE
)
colnames(dB_AND) <- LETTERS[1:4]

# do_HyperTraPS(dB_AND, "HyperTraPS_examples/HP_AND", runs = 500, bi=200, dry_run = TRUE)
sample_freqs(dB_AND, "HyperTraPS_examples/HP_AND/freqs.jpg")

dB_OR <- matrix(
  c(
    rep(c(1, 0, 0, 0), 200)
    , rep(c(1, 0, 1, 0), 100)
    , rep(c(1, 1, 0, 0), 100)
    , rep(c(1, 1, 1, 0), 50)
    , rep(c(1, 1, 0, 1), 50)
    , rep(c(1, 0, 1, 1), 50)
    , rep(c(1, 1, 1, 1), 10)
    , rep(c(0, 0, 0, 0), 10)
  ), ncol = 4, byrow = TRUE
)
colnames(dB_OR) <- LETTERS[1:4]

# do_HyperTraPS(dB_OR, "HyperTraPS_examples/HP_OR", runs = 500, bi=200,  dry_run = TRUE)
sample_freqs(dB_OR, "HyperTraPS_examples/HP_OR/freqs.jpg")

dB_XOR <- matrix(
  c(
    rep(c(1, 0, 0, 0), 200)
    , rep(c(1, 1, 0, 0), 100)
    , rep(c(1, 0, 1, 0), 100)
    , rep(c(1, 1, 0, 1), 50)
    , rep(c(1, 0, 1, 1), 50)
    , rep(c(0, 0, 0, 0), 10)
  ), ncol = 4, byrow = TRUE
)
colnames(dB_XOR) <- LETTERS[1:4]

# do_HyperTraPS(dB_XOR, "HyperTraPS_examples/HP_XOR", runs = 500, bi=200, dry_run = TRUE)
sample_freqs(dB_XOR, "HyperTraPS_examples/HP_XOR/freqs.jpg")

# A --> B ; C --> D; B XOR D for E 

dB_c1 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), 300) #A
    , rep(c(0, 0, 1, 0, 0), 300) #C
    , rep(c(1, 1, 0, 0, 0), 200) #AB
    , rep(c(0, 0, 1, 1, 0), 200) #CD
    , rep(c(1, 1, 1, 0, 0), 100) #ABC
    , rep(c(1, 0, 1, 1, 0), 100) #ACD
    , rep(c(1, 1, 0, 0, 1), 100) #ABE
    , rep(c(0, 0, 1, 1, 1), 100) #CDE
    , rep(c(1, 1, 1, 0, 1), 100) #ABCE
    , rep(c(1, 0, 1, 1, 1), 100) #ACDE
    , rep(c(1, 1, 1, 1, 0), 50) # ABCD
    , rep(c(0, 0, 0, 0, 0), 10) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c1) <- LETTERS[1:5]

# do_HyperTraPS(dB_c1, "HyperTraPS_examples/HP_c1", runs = 500, bi=200, dry_run = TRUE)
sample_freqs(dB_c1, "HyperTraPS_examples/HP_c1/freqs.jpg")
 
##Considerer here another examples to run

# ((A AND B) or C) to reach D

dB_c2 <- matrix(
  c(
      rep(c(1, 0, 0, 0), 300) #A
    , rep(c(0, 1, 0, 0), 300) #B
    , rep(c(0, 0, 1, 0), 300) #C
    , rep(c(1, 1, 0, 0), 200) #AB
    , rep(c(1, 0, 1, 0), 200) #AC
    , rep(c(0, 1, 1, 0), 200) #BC
    , rep(c(1, 1, 1, 0), 100) #ABC
    , rep(c(1, 1, 0, 1), 200) #ABD
    , rep(c(0, 0, 1, 1), 250) #CD
    , rep(c(1, 1, 1, 1), 200) # ABCD
    , rep(c(0, 0, 0, 0), 10) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c2) <- LETTERS[1:4]

# do_HyperTraPS(dB_c2, "HyperTraPS_examples/HP_c2", runs = 500, bi=200, dry_run = TRUE)
sample_freqs(dB_c2, "HyperTraPS_examples/HP_c2/freqs.jpg")
 
## WT --> A AND B --> E
## XOR 
## WT --> C AND D --> E
## ((A AND B) XOR (C AND D) to reach E
dB_c3 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), 300) #A
    , rep(c(0, 1, 0, 0, 0), 300) #B
    , rep(c(0, 0, 1, 0, 0), 300) #C
    , rep(c(0, 0, 0, 1, 0), 300) #D
    , rep(c(1, 1, 0, 0, 0), 200) #AB
    , rep(c(0, 0, 1, 1, 0), 200) #CD
    , rep(c(1, 1, 0, 0, 1), 100) #ABE
    , rep(c(0, 0, 1, 1, 1), 100) #CDE
    , rep(c(0, 0, 0, 0, 0), 10) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c3) <- LETTERS[1:5]

# do_HyperTraPS(dB_c3, "HyperTraPS_examples/HP_c3", runs = 500, bi=200, dry_run = TRUE)
sample_freqs(dB_c3, "HyperTraPS_examples/HP_c3/freqs.jpg")

## WT --> A --> B
## WT --> C --> D
dB_c4 <- matrix(
  c(
      rep(c(1, 0, 0, 0), 300) #A
    , rep(c(0, 0, 1, 0), 300) #C
    , rep(c(1, 1, 0, 0), 200) #AB
    , rep(c(0, 0, 1, 1), 200) #CD
    , rep(c(1, 1, 1, 0), 100) #ABC
    , rep(c(1, 0, 1, 1), 100) #ACD
    , rep(c(1, 1, 1, 1), 50) #ABCD
    , rep(c(0, 0, 0, 0), 10) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c4) <- LETTERS[1:4]

# do_HyperTraPS(dB_c4, "HyperTraPS_examples/HP_c4", runs = 500, bi=200, dry_run = TRUE)
sample_freqs(dB_c4, "HyperTraPS_examples/HP_c4/freqs.jpg")

## WT --> (A AND B AND C) --> D

dB_c5 <- matrix(
  c(
      rep(c(1, 0, 0, 0), 300) #A
    , rep(c(0, 1, 0, 0), 300) #B
    , rep(c(0, 0, 1, 0), 300) #C
    , rep(c(1, 1, 0, 0), 200) #AB
    , rep(c(1, 0, 1, 0), 200) #AC
    , rep(c(0, 1, 1, 0), 200) #BC
    , rep(c(1, 1, 1, 0), 100) #ABC
    , rep(c(1, 1, 1, 1), 50)  #ABCD
    , rep(c(0, 0, 0, 0), 10) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c5) <- LETTERS[1:4]

# do_HyperTraPS(dB_c5, "HyperTraPS_examples/HP_c5", runs = 500, bi=200, dry_run = TRUE)
sample_freqs(dB_c5, "HyperTraPS_examples/HP_c5/freqs.jpg")

## Sign epistasis
dB_se <- matrix(
  c(
      rep(c(1, 0), 50) #A
    , rep(c(0, 1), 300) #B
    , rep(c(1, 1), 400) #AB
    , rep(c(0, 0), 200) # WT
  ), ncol = 2, byrow = TRUE
)
colnames(dB_se) <- LETTERS[1:2]
do_HyperTraPS(dB_se, "HyperTraPS_examples/HP_se", runs = 500, bi=200, dry_run = TRUE)
sample_freqs(dB_se, "HyperTraPS_examples/HP_se/freqs.jpg")

## Reciprocal sign epistasis
dB_rse <- matrix(
  c(
      rep(c(1, 0), 50) #A
    , rep(c(0, 1), 50) #B
    , rep(c(1, 1), 400) #AB
    , rep(c(0, 0), 200) # WT
  ), ncol = 2, byrow = TRUE
)
colnames(dB_rse) <- LETTERS[1:2]
do_HyperTraPS(dB_rse, "HyperTraPS_examples/HP_rse", runs = 500, bi=200, dry_run = FALSE)
sample_freqs(dB_rse, "HyperTraPS_examples/HP_rse/freqs.jpg")


##Considerer here another examples to run

## ((A AND B) OR (C AND D) to reach E

dB_c7 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), 300) #A
    , rep(c(0, 1, 0, 0, 0), 300) #B
    , rep(c(0, 0, 1, 0, 0), 300) #C
    , rep(c(0, 0, 0, 1, 0), 300) #D
    , rep(c(1, 1, 0, 0, 0), 200) #AB
    , rep(c(1, 0, 1, 0, 0), 200) #AC
    , rep(c(1, 0, 0, 1, 0), 200) #AD
    , rep(c(0, 1, 1, 0, 0), 200) #BC
    , rep(c(0, 1, 0, 1, 0), 200) #BD
    , rep(c(0, 0, 1, 1, 0), 200) #CD
    , rep(c(1, 1, 1, 0, 0), 100) #ABC
    , rep(c(1, 1, 0, 1, 0), 100) #ABD
    , rep(c(0, 1, 1, 1, 0), 100) #BCD
    , rep(c(1, 1, 0, 0, 1), 200) #ABE
    , rep(c(0, 0, 1, 1, 1), 200) #CDE
    , rep(c(1, 1, 1, 0, 1), 200) #ABCE
    , rep(c(1, 1, 0, 1, 1), 200) # ABDE
    , rep(c(1, 0, 1, 1, 1), 200) # ACDE
    , rep(c(0, 1, 1, 1, 1), 200) # BCDE
    , rep(c(1, 1, 1, 1, 1), 200) # ABCDE
    , rep(c(0, 0, 0, 0, 0), 10) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c7) <- LETTERS[1:5]

# do_HyperTraPS(dB_c7, "HyperTraPS_examples/HP_c7", runs = 500, bi=200, dry_run = TRUE)
sample_freqs(dB_c7, "HyperTraPS_examples/HP_c7/freqs.jpg")

## Make a test with the real world data of the pmce paper

# do_HyperTraPS("Bladder_Urothelial_Carcinoma.csv", "HyperTraPS_examples/HP_pmce", runs = 1000, bi=500)
# sample_freqs(dB_OR, "HyperTraPS_examples/HP_OR/freqs.jpg")
