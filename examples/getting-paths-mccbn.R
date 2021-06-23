library(mccnb)

set.seed(1)
true_p1 <- mccbn::random_poset(6)
true_p1
lambda_s <- 1
lambdas <- runif(6, 1/6*lambda_s, 6*lambda_s)
set.seed(1)
simGenotypes <- mccbn::sample_genotypes(10, true_p1,
                                        sampling_param = lambda_s,
                                        lambdas = lambdas)

simGenotypes$T_sampling

simGenotypes$T_events

## The first subject has been sampled at 0.7552
## the mutations take place at less than that time could have taken place
## (we also need to check T_sum_events?)
simGenotypes$T_events[1, ] < simGenotypes$T_sampling[1]
## all except the 2nd are OK
## And, from $T_sum_events
## we see all, except the second, are < simGenotypes$T_sampling[1]

## In what order? As given by the mutation time
order(simGenotypes$T_events[1, ])
## So the path is 6, 5, 4, 1, 3.   2 will not happen.

## Let's repeat this with the 3rd individual

simGenotypes$T_events[3, ] < simGenotypes$T_sampling[3]
## only genes 3 and 6

order(simGenotypes$T_events[3, ])
## First 3 and then 6


## OJO! We would need to check that our output of paths is identical to the
## genotype given.
