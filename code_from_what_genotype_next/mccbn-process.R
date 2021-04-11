## Recall that for using newer versions of MC-CBN you need libboost < 1.74:
## https://github.com/cbg-ethz/MC-CBN/issues/5
## I install as follows
## sudo apt-get install libboost1.67-dev:amd64 libboost1.67-tools-dev    libboost-graph-parallel1.67-dev
## Which removes newer versions.
## A minimal wrapper. No bootstrap or anything like that for now.
library(mccbn)

mccbn_proc <- function(x) {
    stopifnot(!is.null(colnames(x)))
    fit <- mccbn::learn_network(x)
    mle_index <- which.max(fit$logliks)
    am <- fit$posets[[mle_index]]
    df1 <- igraph::as_data_frame(graph_from_adjacency_matrix(am))
    colnames(df1) <- c("From", "To")
    no_parent <- setdiff(colnames(x), df1[, 2])
    dfr <- rbind(
        data.frame(From = "Root", To = no_parent,
                   stringsAsFactors = FALSE),
        df1)
    dfr$edge = paste(dfr[, "From"],
                     dfr[, "To"],
                     sep = " -> ")
    ## of course, lambda is per child

    ## if using MC-CBN pre 2020-10-07 
    ## lambda <- fit$fits[[mle_index]]$par

    ## if using MC-CBN post 2020-10-07 
    lambda <- fit$fits[[mle_index]]$lambda

    names(lambda) <- colnames(am)
    dfr$lambda <- lambda[dfr$To]
    return(list(edges = dfr))
}


## ## Using the new functionality. Extremely slow.
## ## The only non-defaults are threads (thrds)
## ## and L.
## posets <- mccbn::candidate_posets(db2, rep(1, nrow(db2)), 0.9)
## poset0 <- posets[[length(posets)]]

## fit <- mccbn::adaptive.simulated.annealing(poset0,
##                                            db2, L = 100,
##                                            max.iter.asa = 10000L,
##                                            thrds = 8)
