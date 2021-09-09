library(mccbn)

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
simGenotypes$T_sum_events[1, ] < simGenotypes$T_sampling[1]
## all except the 2nd are OK
## And, from $T_sum_events
## we see all, except the second, are < simGenotypes$T_sampling[1]

## In what order? As given by the mutation time

order(simGenotypes$T_sum_events[1, ])
## Path is 6 1 3 5 4 

## Let's repeat this with the 3rd individual

simGenotypes$T_sum_events[3, ] < simGenotypes$T_sampling[3]
## only genes 3 and 6

## First 3 and then 6


## OJO! We would need to check that our output of paths is identical to the
## genotype given.



## yes, look at T_sum_events: T_events is time given ancestors satisfied. So total time is T_sum_events.

## Look at example for four. This is in Spanish

## > simGenotypes$T_events[4, ]
## [1] 0.8096 0.4213 0.2442 0.2207 0.1387 0.4097

## > simGenotypes$T_sum_events[4, ]
## [1] 0.8096 0.4213 0.2442 1.0303 0.9483 0.4097


## De hecho, el cuarto sujeto no tiene nada mutado y su tiempo es 

## > simGenotypes$T_sampling[4]
## [1] 0.1398


## Creo que el T_events es time to events given that parents are
## satisfied. Eso tiene sentido. Fíjate el tiempo en el T_sum_events, con el
## pose que hemos creado, donde los genes 4 y 5 dependen del 1:

## El T_sum_events del 4º gen = 0.8096 + .2207

## El T_sum_events del 5º gen = 0.8096 + .1387
