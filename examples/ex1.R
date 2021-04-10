library(dplyr)
library(Oncotree)
library(readr)
library(data.table)
library(igraph)
## library(TRONCO)
## library(parallel)
## library(foreach)
## library(doParallel)
library(Rgraphviz)
library(stringr)
library(pryr)
library(OncoSimulR)
library(testthat)
library(Matrix)

library(plot.matrix) ## for the plot of the theta matrix


pwd0 <- getwd()
setwd("../code_from_what_genotype_next/")
source("code-all-methods-minimal.R")
source("schill-trans-mat.R") ## slow: we run many tests
source("pre-process.R")
source("ot-process.R")
source("cbn-process.R")
setwd(pwd0)
rm(pwd0)

igraph_options(vertex.size = 14)
igraph_options(vertex.label.cex = 1.5)
igraph_options(edge.label.cex = 1.5)




## 
N <- 200
na <- N
nb <- N + round( 10 * runif(1))
nab <- N + round( 10 * runif(1))
n00 <- N/10 + round( 10 * runif(1))
dB <- matrix(
    c(
        rep(c(1, 0), na) 
      , rep(c(0, 1), nb)
      , rep(c(1, 1), nab)
      , rep(c(0, 0), n00)
    ), ncol = 2, byrow = TRUE
)
colnames(dB) <- LETTERS[1:2]
storage.mode(dB) <- "integer"
sampledGenotypes(dB)

out <- all_methods_2_trans_mat(dB)
plot_DAG_fg(out, dB)




N <- 200
na <- N
nb <- N + round( 10 * runif(1))
nab <- 1.6 * N + round( 10 * runif(1))
n00 <- round( 10 * runif(1))
dB <- matrix(
    c(
        rep(c(1, 0), na) 
      , rep(c(0, 1), nb)
      , rep(c(1, 1), nab)
      , rep(c(0, 0), n00)
    ), ncol = 2, byrow = TRUE
)
colnames(dB) <- LETTERS[1:2]
out <- all_methods_2_trans_mat(dB)
plot_DAG_fg(out, dB)



N <- 100
na <- N
nc <- N + round( 10 * runif(1))
nab <- 1.6 * N + round( 10 * runif(1))
nac <- 1.6 * N + round( 10 * runif(1))
nabc <- 1.4 * N + round( 10 * runif(1))
n00 <- round( 10 * runif(1))
dB <- matrix(
    c(
        rep(c(1, 0, 0), na) 
      , rep(c(0, 0, 1), nc)
      , rep(c(1, 1, 0), nab)
      , rep(c(1, 0, 1), nac)        
      , rep(c(1, 1, 1), nabc)
      , rep(c(0, 0, 0), n00)
    ), ncol = 3, byrow = TRUE
)
colnames(dB) <- LETTERS[1:3]

out <- all_methods_2_trans_mat(dB)
plot_DAG_fg(out, dB)
