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
lambdas <- runif(6, 1/2*lambda_s, lambda_s)
simGenotypes <- mccbn::sample_genotypes(20000, true_p1,
                                        sampling_param = lambda_s,
                                        lambdas = lambdas)
db2 <- simGenotypes$obs_events
colnames(db2) <- LETTERS[1:6]
sampledGenotypes(db2)
out <- all_methods_2_trans_mat(db2)
plot_DAG_fg(out, db2)

db3 <- remove_WT(db2, 1)
out3 <- all_methods_2_trans_mat(db3)
plot_DAG_fg(out3, db3)

db4 <- add_WT(db2, 100000)
out4 <- all_methods_2_trans_mat(db4)
plot_DAG_fg(out4, db4)



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

db3 <- remove_WT(db2, 1  )
out3 <- all_methods_2_trans_mat(db3)
plot_DAG_fg(out3, db3)





## We cannot simulate a case where there both A and B strong
## effect, but joint increases fitness very little
## MHN can capture that
## Simulating data under posets
true_p1 <- matrix(0, nrow = 4, ncol = 4)
true_p1[3, 4] <- 1
lambda_s <- 1
simGenotypes <- mccbn::sample_genotypes(20000, true_p1,
                                        sampling_param = 1,
                                        lambdas = c(5, 5, .5, .4))
db2 <- simGenotypes$obs_events
colnames(db2) <- LETTERS[1:4]
sampledGenotypes(db2)
out <- all_methods_2_trans_mat(db2)
plot_DAG_fg(out, db2)


## Two great examples!!

## Why MHN lacks any mechanistic/causal interpretation.

## We get the negative from MHN, but transitions
## are screwed up

## This is a case where A and B each has huge fitness effects but further
## additions do not add much.  D depends on C.

## Interesting that we get a negative coeff.
## for C on A and B. But mild for A or B on C.
## Should be the other way around.

## Also: Very similar estimates of transitions, with some interesting
## diffrences.


N <- 1000
na <- 4 * N
nb <- 4* N + round( 100 * runif(1))
nc <- N + round( 100 * runif(1))
nab <- N + round( 100 * runif(1))
nac <- N + round( 100 * runif(1))
nbc <- N + round( 100 * runif(1))
ncd <- N + round( 100 * runif(1))
nabc <- N + round( 100 * runif(1))
nacd <- N + round( 100 * runif(1))
nbcd <- N + round( 100 * runif(1))
nabcd <- N + round( 100 * runif(1))
n00 <- 20 # und( 10 * runif(1))
dB <- matrix(
    c(
        rep(c(1, 0, 0, 0), na) 
      , rep(c(0, 1, 0, 0), nb)
      , rep(c(0, 0, 1, 0), nc)
      , rep(c(1, 1, 0, 0), nab)
      , rep(c(1, 0, 1, 0), nac)
      , rep(c(0, 1, 1, 0), nbc)
      , rep(c(0, 0, 1, 1), ncd)
      , rep(c(1, 1, 1, 0), nabc)
      , rep(c(1, 0, 1, 1), nacd)
      , rep(c(0, 1, 1, 1), nbcd)
      , rep(c(1, 1, 1, 1), nabcd)      
      , rep(c(0, 0, 0, 0), n00)
    ), ncol = 4, byrow = TRUE
)
colnames(dB) <- LETTERS[1:4]
out <- all_methods_2_trans_mat(dB)
plot_DAG_fg(out, dB)






## Note that highly positive D on C effect, but it is C that depends on D:
##   very high \theta_C_D

N <- 1000
na <- 12 * N
nb <- 13 * N + round( 100 * runif(1))
nc <- 9 * N + round( 100 * runif(1))
nab <- .0 * N + round( 100 * runif(1))
nac <- .5 * N + round( 100 * runif(1))
nbc <- .5 * N + round( 100 * runif(1))
ncd <- N + round( 100 * runif(1))
nabc <- .1 * N + round( 100 * runif(1))
nacd <- N + round( 100 * runif(1))
nbcd <- N + round( 100 * runif(1))
nabcd <- N + round( 100 * runif(1))
n00 <- 20 # und( 10 * runif(1))
dB <- matrix(
    c(
        rep(c(1, 0, 0, 0), na) 
      , rep(c(0, 1, 0, 0), nb)
      , rep(c(0, 0, 1, 0), nc)
      , rep(c(1, 1, 0, 0), nab)
      , rep(c(1, 0, 1, 0), nac)
      , rep(c(0, 1, 1, 0), nbc)
      , rep(c(0, 0, 1, 1), ncd)
      , rep(c(1, 1, 1, 0), nabc)
      , rep(c(1, 0, 1, 1), nacd)
      , rep(c(0, 1, 1, 1), nbcd)
      , rep(c(1, 1, 1, 1), nabcd)      
      , rep(c(0, 0, 0, 0), n00)
    ), ncol = 4, byrow = TRUE
)
colnames(dB) <- LETTERS[1:4]
out <- all_methods_2_trans_mat(dB)
plot_DAG_fg(out, dB)


