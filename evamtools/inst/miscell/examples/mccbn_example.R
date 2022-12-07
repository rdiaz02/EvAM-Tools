library(evamtools)
library(mccbn)
data <- examples_csd$csd$AND$data

evamtools:::mccbn_hcbn_proc(data, max.iter.asa = 10,
    , ncores = 4)

posets <- mccbn::candidate_posets(data, rep(1, nrow(data)), 0.9)
poset0 <- posets[[length(posets)]]
fit <- mccbn::adaptive.simulated.annealing(poset0, data, L=100,
                                           max.iter.asa=10, seed=10L)
