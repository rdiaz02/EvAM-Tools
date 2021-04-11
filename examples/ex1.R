pwd0 <- getwd()
setwd("../code_from_what_genotype_next/")
source("code-all-methods-minimal.R")
setwd(pwd0)
rm(pwd0)

igraph_options(vertex.size = 14)
igraph_options(vertex.label.cex = 1.5)
igraph_options(edge.label.cex = 1.5)


## Remove a fraction of the WT
rm_wt <- function(x, rm = 0.9) {
    which_wt <- which(rowSums(x) == 0)
    if(length(which_wt) > 0) {
        rmwt <- which_wt[1:(round(rm * length(which_wt)))]
    }
    return(x[-rmwt, ])
}


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


## Very weird behaviour with CBN
N <- 100
na <- N
nc <- 2*N  + round( 10 * runif(1))
nb <- 0
nab <- 1.6 * N + round( 10 * runif(1))
nac <- 2.3 * N + round( 10 * runif(1))
nabc <- N + round( 10 * runif(1))
n00 <- 5  ## 1000 + round( 10 * runif(1))
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
N <- 100
na <- N
ne <- N + round(10 * runif(1))
nab <- 2*N  + round( 10 * runif(1))
nabc <- N
nabcd <- 1.6 * N + round( 10 * runif(1))
nef <- 2.3 * N + round( 10 * runif(1))
naef <- 2.3 * N + round( 10 * runif(1))
nae <- 1.9 * N + round( 10 * runif(1))

nabc <- N + round( 10 * runif(1))
n00 <- 5  ## 1000 + round( 10 * runif(1))
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




###################################
## Simulating data under posets

## OT fails in most simple cases
## But note how the samples look in terms of genotype freqs.
## This emphasizes role of the WT, specially for CBN


true_p1 <- matrix(0, nrow = 3, ncol = 3)
true_p1[1, 2] <- 1
lambda_s <- 1
lambdas <- runif(3, 1/3*lambda_s, 3*lambda_s)
simGenotypes <- mccbn::sample_genotypes(700, true_p1,
                                        sampling_param = lambda_s,
                                        lambdas = lambdas)
db2 <- simGenotypes$obs_events
colnames(db2) <- LETTERS[1:3]
sampledGenotypes(db2)
out <- all_methods_2_trans_mat(db2)
plot_DAG_fg(out, db2)





## Simulating data under posets
true_p1 <- matrix(0, nrow = 4, ncol = 4)
true_p1[1, 2] <- 1
true_p1[3, 4] <- 1
lambda_s <- 1
lambdas <- runif(4, 1/3*lambda_s, 3*lambda_s)
simGenotypes <- mccbn::sample_genotypes(20000, true_p1,
                                        sampling_param = lambda_s,
                                        lambdas = lambdas)
db2 <- simGenotypes$obs_events
colnames(db2) <- LETTERS[1:4]
sampledGenotypes(db2)
out <- all_methods_2_trans_mat(db2)
plot_DAG_fg(out, db2)





## Simulating data under posets
true_p1 <- matrix(0, nrow = 6, ncol = 6)
true_p1[1, 2] <- 1
true_p1[2, 3] <- 1
true_p1[3, 4] <- 1
true_p1[5, 6] <- 1
lambda_s <- 1
lambdas <- runif(6, 1/6*lambda_s, 6*lambda_s)
simGenotypes <- mccbn::sample_genotypes(20000, true_p1,
                                        sampling_param = lambda_s,
                                        lambdas = lambdas)
db2 <- simGenotypes$obs_events
colnames(db2) <- LETTERS[1:6]
sampledGenotypes(db2)
out <- all_methods_2_trans_mat(db2)
plot_DAG_fg(out, db2)

db3 <- rm_wt(db2, 1)
out3 <- all_methods_2_trans_mat(db3)
plot_DAG_fg(out3, db3)

## Simulating data under posets
true_p1 <- matrix(0, nrow = 5, ncol = 5)
true_p1[1, 2] <- 1
true_p1[2, 3] <- 1
true_p1[4, 5] <- 1
lambda_s <- 1
lambdas <- runif(5, 1/5*lambda_s, 5*lambda_s)
simGenotypes <- mccbn::sample_genotypes(20000, true_p1,
                                        sampling_param = lambda_s,
                                        lambdas = lambdas)
db2 <- simGenotypes$obs_events
colnames(db2) <- LETTERS[1:5]
sampledGenotypes(db2)
out <- all_methods_2_trans_mat(db2)
plot_DAG_fg(out, db2)

db3 <- rm_wt(db2, 1  )
out3 <- all_methods_2_trans_mat(db3)
plot_DAG_fg(out3, db3)
