
## We cannot simulate a case where there both A and B strong effect, but joint
## increases fitness very little We want a case of negative epistasis (not of
## reciprocal sign epistasis). This simply shows that the simulation of data
## under posets, by how it is done (but what is the CBN model) does not allow
## for these effects.

## MHN can capture that, but would not warrant a causal interpretation (this is
## selection bias).


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
out <- evam(db2)
plot_DAG_fg(out, db2)
