pwd0 <- getwd()
setwd("../code_from_what_genotype_next/")
source("code-all-methods-minimal.R")
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




N <- 50
na <- N
nb <- N + round( 10 * runif(1))
nab <- .6 * N + round( 10 * runif(1))
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



## Very weird behaviour with CBN or with OT
## very sensitive to the n00
N <- 100
na <- N + round( 10 * runif(1))
nc <- .73 * N  + round( 10 * runif(1))
nab <- .8 * N + round( 10 * runif(1))
nac <- 1.6 * N + round( 10 * runif(1))
nabc <- 1.2 * N + round( 10 * runif(1))
n00 <- 19 ## 1.95 * N  ## 1000 + round( 10 * runif(1))
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

db3 <- remove_WT(dB, 1)
db4 <- add_WT(db3, 10 * nrow(db3))

out3 <- all_methods_2_trans_mat(db3)
out4 <- all_methods_2_trans_mat(db4)

plot_DAG_fg(out3, db3)
plot_DAG_fg(out4, db4)

## Very weird behaviour with CBN
N <- 100
na <- N
nc <- 2*N  + round( 10 * runif(1))
nb <- 0
nab <- 1.6 * N + round( 10 * runif(1))
nac <- 2.3 * N + round( 10 * runif(1))
nabc <- N + round( 10 * runif(1))
n00 <- 10  ## 1000 + round( 10 * runif(1))
dB <- matrix(
    c(
        rep(c(1, 0, 0), na)
      , rep(c(0, 1, 0), nb)         
      , rep(c(0, 0, 1), nc)
      , rep(c(1, 1, 0), nab)
      , rep(c(1, 0, 1), nac)        
      , rep(c(1, 1, 1), nabc)
      , rep(c(0, 0, 0), n00)
    ), ncol = 3, byrow = TRUE
)
colnames(dB) <- LETTERS[1:3]
## dB <- dB[, sample(ncol(dB))]
## dB <- dB[sample(nrow(dB)), ]
out <- all_methods_2_trans_mat(dB)
plot_DAG_fg(out, dB)



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
out <- all_methods_2_trans_mat(dB)
plot_DAG_fg(out, dB)



## Very weird behaviour with CBN
## very sensitive to the n00
N <- 100
na <- N + round( 10 * runif(1))
nc <- N  + round( 10 * runif(1))
nabc <- 1.2 * N + round( 10 * runif(1))
n00 <- 19 ## 1.95 * N  ## 1000 + round( 10 * runif(1))
dB <- matrix(
    c(
        rep(c(1, 0, 0), na) 
      , rep(c(0, 0, 1), nc)
      , rep(c(1, 1, 1), nabc)
      , rep(c(0, 0, 0), n00)
    ), ncol = 3, byrow = TRUE
)
colnames(dB) <- LETTERS[1:3]
db3 <- remove_WT(dB, 1)
db4 <- add_WT(db3, 10 * nrow(db3))
db5 <- add_WT(db3, 5 * nrow(db3))



## Very weird behaviour with CBN
## very sensitive to the n00
N <- 100
na <- N + round( 10 * runif(1))
nc <- N  + round( 10 * runif(1))
nac <- .5 * N  + round( 10 * runif(1))
nabc <- .2 * N + round( 10 * runif(1))
n00 <- 19 ## 1.95 * N  ## 1000 + round( 10 * runif(1))
dB <- matrix(
    c(
        rep(c(1, 0, 0), na) 
      , rep(c(0, 0, 1), nc)
      , rep(c(1, 0, 1), nac)
      , rep(c(1, 1, 1), nabc)
      , rep(c(0, 0, 0), n00)
    ), ncol = 3, byrow = TRUE
)
colnames(dB) <- LETTERS[1:3]
db3 <- remove_WT(dB, 1)
db4 <- add_WT(db3, 10 * nrow(db3))
db5 <- add_WT(db3, 5 * nrow(db3))










