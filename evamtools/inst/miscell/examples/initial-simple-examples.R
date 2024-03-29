## One of a set of initial examples that can help
## get an idea of some interesting patterns.
## Function plot_DAG_fg no longer works, but can be replaced
## by the current working function.
## FIXME: do that: replace plot_DAG_fg by current working function.


igraph_options(vertex.size = 14)
igraph_options(vertex.label.cex = 1.5)
igraph_options(edge.label.cex = 1.5)


## I would say CBN is not working properly 
## it also gives different results depending of the run
## we have 200 A and 200 B and 200 AB 
## I would expect A and B are independent 
## but CBN sometimes says that B comes before A or vice versa
N <- 200
na <- N
nb <- N + round( 10 * runif(1))
nab <- N + round( 10 * runif(1))
n00 <- N/10 + round( 10 * runif(1))
dB <- matrix(
    c(
        rep(c(0, 1), nb)
      , rep(c(1, 0), na) 
      , rep(c(1, 1), nab*3)
      , rep(c(0, 0), n00)
    ), ncol = 2, byrow = TRUE
)
colnames(dB) <- LETTERS[1:2]
storage.mode(dB) <- "integer"
data_to_counts(dB, out = "data.frame", omit_0 = TRUE)
out <- evam(dB)
## plot_DAG_fg(out, dB)



## OT and MHN have identical transition matrices
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
out <- evam(dB)
## plot_DAG_fg(out, dB)



## Very weird behaviour with CBN or with OT
## very sensitive to the n00
## --> but OT gives the right result, CBN is the one failing.
## --> changing n00 we get different trees with ot and cbn.
## Why do we see theses changes?
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
out <- evam(dB)

db3 <- remove_WT(dB, 1)
db4 <- add_WT(db3, 10 * nrow(db3))
db5 <- add_WT(db3, floor(na/2))

par(mfrow=c(1,4))
plot_sampled_genots(dB)
plot_sampled_genots(db3)
plot_sampled_genots(db4)
plot_sampled_genots(db5)
par(mfrow=c(1,1))
out3 <- evam(db3)
out4 <- evam(db4)
out5 <- evam(db5)

## CBN and OT have very similar transition matrices
## regardless of the WT frequency
## plot_DAG_fg(out, dB, matrix=TRUE)
## plot_DAG_fg(out3, db3)
## plot_DAG_fg(out4, db4) 
## plot_DAG_fg(out5, db5)
## db4 is quite rare:
## it goes from 
## WT --> a -->b
##    --> c
## to
## WT --> a -->b
##          --> c

## Very weird behaviour with CBN
## the DAG is ok but the transition matrix is not 
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
out <- evam(dB)
## plot_DAG_fg(out, dB)


## What level of theta_ij would be indicative of 
## mutual exclusivity
## This is a toy model for mutual exclusivity: 
## we either go  with one branch or the otherversion
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
out <- evam(dB)
## plot_DAG_fg(out, dB)



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








