## Copyright 2021, 2022 Pablo Herrera Nieto, Ramon Diaz-Uriarte

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


## This is using OT-CBN, not H-CBN2
do_MCCBN_OT_CBN <- function(x) {
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

    lambda <- fit$fits[[mle_index]]$lambda

    names(lambda) <- colnames(am)
    dfr$lambda <- lambda[dfr$To]
    return(list(edges = dfr))
}


## #' @title Run MCCBN with H-CBN2
## #' 
## #' Using H-CBN2, as in example from https://github.com/cbg-ethz/MC-CBN
## #' 
## #' @param x Cross secitonal data. Matrix of genes (columns)
## #' and individuals (rows)
## #' @param max.iter.asa Int. Number of steps to run. Default: 100000
## #' @param ncores Int. Number of threads to use
## #' @param L Int. Number of samples to be drawn from the proposal in the E-step
## #' @param tmp_dir Directory name where the oput is located. 
## #' @param addname String to append to the temporary folder name. Default NULL
## #' @param silent Whether to run show message showing the folder name where MCCBN is run
## #' 
## #' @return A list with the adjacency matrix, the lambdas, the parent set
## #' and a data.frame with From-To edges and associated lambdas.

do_MCCBN_HCBN2 <- function(x 
                          , mccbn_hcbn2_opts = list(
                                tmp_dir = NULL,
                                addname = NULL,
                                silent = TRUE,
                                L = 100,
                                sampling = c("forward", "add-remove", "backward", "bernoulli", "pool"),
                                max.iter = 100L,
                                update.step.size = 20L,
                                tol = 0.001,
                                max.lambda.val = 1e+06,
                                T0 = 50,
                                adap.rate = 0.3,
                                acceptance.rate = NULL,
                                step.size = NULL,
                                max.iter.asa = 10000L,
                                neighborhood.dist = 1L,
                                adaptive = TRUE,
                                thrds = 1L,
                                verbose = FALSE,
                                seed = NULL
                           )
                           ) {
    if (!requireNamespace("mccbn", quietly = TRUE))
        stop("MC-CBN (mccbn) no installed")
    
    stopifnot(!is.null(colnames(x)))
    stopifnot(mccbn_hcbn2_opts$max.iter.asa >= 5)

    # Setting tmp folder
    # The functiona create asa.txt and poset.txt
    if (is.null(mccbn_hcbn2_opts$tmp_dir)) {
        tmp_dir <- tempfile()
        dirname0 <- NULL
        if (!is.null(mccbn_hcbn2_opts$addname)) {
            dirname0 <- tmp_dir
            tmp_dir <- paste0(tmp_dir, "/",
                              "_mccbn_", mccbn_hcbn2_opts$addname)
        }
        if (!mccbn_hcbn2_opts$silent)
            message(paste("\n Using dir", tmp_dir))
        if (dir.exists(tmp_dir)) {
            stop("dirname ", tmp_dir, "exists")
        }
        dir.create(tmp_dir, recursive = TRUE)
    } else {
        tmp_dir <- mccbn_hcbn2_opts$tmp_dir
    }
    posets <- mccbn::candidate_posets(x, rep(1, nrow(x)), 0.9)
    poset0 <- posets[[length(posets)]]
    fit <- mccbn::adaptive.simulated.annealing(poset = poset0,
                                               obs = x,
                                               outdir = tmp_dir,
                                               L = mccbn_hcbn2_opts$L,
                                               sampling = mccbn_hcbn2_opts$sampling,
                                               max.iter = mccbn_hcbn2_opts$max.iter,
                                               update.step.size = mccbn_hcbn2_opts$update.step.size,
                                               tol = mccbn_hcbn2_opts$tol,
                                               max.lambda.val = mccbn_hcbn2_opts$max.lambda.val,
                                               T0 = mccbn_hcbn2_opts$T0,
                                               adap.rate = mccbn_hcbn2_opts$adap.rate,
                                               acceptance.rate = mccbn_hcbn2_opts$acceptance.rate,
                                               step.size = mccbn_hcbn2_opts$step.size,
                                               max.iter.asa = mccbn_hcbn2_opts$max.iter.asa,
                                               neighborhood.dist = mccbn_hcbn2_opts$neighborhood.dist,
                                               adaptive = mccbn_hcbn2_opts$adaptive,
                                               thrds = mccbn_hcbn2_opts$thrds,
                                               verbose = mccbn_hcbn2_opts$verbose,
                                               seed = mccbn_hcbn2_opts$seed
                                               )

    ## Default iterations for asa is 10000 in original code, 100 for testing
    ## L: ; they used 100, but I do not know if this is too low

    am <- fit$poset
    colnames(am) <- rownames(am) <- colnames(x)
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
    lambda <- fit$lambda

    names(lambda) <- colnames(am)
    dfr$lambda <- lambda[dfr$To]
    return(list(edges = dfr))
}




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
