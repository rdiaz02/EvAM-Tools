## See file OT_transition_matrices.org for the place
## where these numbers are used.
library(Oncotree)
library(OncoSimulR)

pa <- 0.3
pb <- 0.6
N <- 1000
na <- N * pa * (1 - pb)
nb <- N * pb * (1 - pa)
nab <- N * pa * pb
n00 <- N * (1 - pa) * (1 - pb)
dB <- matrix(
    c(
        rep(c(0, 1), nb)
      , rep(c(1, 0), na) 
      , rep(c(1, 1), nab)
      , rep(c(0, 0), n00)
    ), ncol = 2, byrow = TRUE
)
colnames(dB) <- LETTERS[1:2]
storage.mode(dB) <- "integer"
sampledGenotypes(dB)

ot.fit <- oncotree.fit(dB)

## Weights. This is a very sample example
## so the next two are identical
## and identical to the pa and pb given.
ot.fit$parent$est.weight
ot.fit$parent$obs.weight

## Predicted genotypes
distribution.oncotree(ot.fit, with.probs = TRUE)



