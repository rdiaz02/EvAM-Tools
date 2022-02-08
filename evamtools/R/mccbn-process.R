## Copyright 2021, 2022 Pablo Herrera Nieto

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Not used now
if(FALSE) {
    ## This is using OT-CBN, not H-CBN2
    mccbn_proc <- function(x) {
        if (!requireNamespace("mccbn", quietly = TRUE))
            stop("MC-CBN (mccbn) no installed")
        
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



######################################################################
######################################################################
###
###  libboost issues
###
######################################################################
######################################################################


## You need to do first
## sudo apt-get install dh-autoreconf autoconf automake autotools-dev autoconf autoconf-archive


## Recall that for using newer versions of MC-CBN you need libboost < 1.74:
## https://github.com/cbg-ethz/MC-CBN/issues/5
## I install as follows
## sudo apt-get install libboost1.67-dev:amd64 libboost1.67-tools-dev    libboost-graph-parallel1.67-dev
## Which removes newer versions.

## Upgrading. These will be removed:

##   libboost-filesystem1.67-dev
##   libboost-graph-parallel1.67-dev libboost-iostreams1.67-dev
##   libboost-locale1.67-dev libboost-regex1.67-dev
##   libboost-serialization1.67-dev libboost-system1.67-dev
##   libboost-test1.67-dev libboost1.67-dev libboost1.67-tools-dev


## Going back to 1.67

## sudo apt-get remove --purge libboost1.74-tools-dev libboost1.74-dev libboost-wave1.74-dev libboost-type-erasure1.74-dev libboost-timer1.74-dev  libboost-test1.74-dev  libboost-stacktrace1.74-dev libboost-random1.74-dev  libboost-python1.74-dev libboost-program-options1.74-dev libboost-numpy1.74-dev  libboost-nowide1.74-dev  libboost-mpi1.74-dev libboost-mpi-python1.74.0  libboost-mpi-python1.74-dev  libboost-mpi-python1.67.0  libboost-math1.74-dev  libboost-log1.74-dev  libboost-locale1.74-dev libboost-graph1.74-dev libboost-iostreams1.74-dev  libboost-graph-parallel1.74-dev  libboost-filesystem1.74-dev libboost-fiber1.74-dev libboost-exception1.74-dev  libboost-date-time1.74-dev libboost-coroutine1.74-dev  libboost-context1.74-dev libboost-container1.74-dev  libboost-chrono1.74-dev libboost-atomic1.74-dev

## sudo apt-get --purge autoremove

## The following 1.74 are needed for other packages
## libboost-filesystem1.74.0:amd64 
## libboost-iostreams1.74.0:amd64  
## libboost-locale1.74.0:amd64 
## libboost-thread1.74.0:amd64

## sudo apt-get install libboost-filesystem1.67-dev libboost-graph-parallel1.67-dev libboost-iostreams1.67-dev libboost-locale1.67-dev libboost-regex1.67-dev libboost-serialization1.67-dev libboost-system1.67-dev libboost-test1.67-dev libboost1.67-dev libboost1.67-tools-dev libboost-graph1.67-dev libboost-graph1.67.0
