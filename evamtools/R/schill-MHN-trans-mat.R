## Copyright 2020, 2022 Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify it under
## the terms of the GNU Affero General Public License (AGPLv3.0) as published by
## the Free Software Foundation, either version 3 of the License, or (at your
## option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero General Public License for more details.

## You should have received a copy of the GNU Affero General Public License along
## with this program.  If not, see <http://www.gnu.org/licenses/>.


## Code to run Schill et al.'s MHN and obtain transition matrices between
## genotypes. Several implementations are provided, some using sparse
## matrices. In general, the best option is to use do_MHN2 for moderate to
## large numbers of features.

## Yes, I implemented several Theta -> transition rate, when Schill
## already provides Build.Q. Oh well, I don't remember what I was thinking.
## Anyway, this is another check. :-)

## observations (rows as patients, columns genes) -> transition matrix genotypes

##          Remember my transition matrices between genotypes have origin
##          in rows, destination in columns. Transposed w.r.t. to their
##          Figure 2, left.

## Identical to do_MHN, but with an argument for sparse, and corresponding
## additional code
do_MHN2 <- function(x,  lambda = 1/nrow(x), sparse = TRUE) {
    ## lambda 0.01 is what they use by default in their bioRxiv, p. 7 paper.
    ## In the paper it is 1/nrow(x). 
    mhnd <- Data.to.pD(x)
    theta <- Learn.MHN(mhnd, lambda = lambda)
    colnames(theta) <- rownames(theta) <- colnames(x)
    ## Reorder, so genotype names and transition rate matrix
    ## have names ordered
    oindex <- order(colnames(theta))
    theta <- theta[oindex, oindex]
    if(!sparse) {
        trm <- theta_to_trans_rate_3(theta,
                                     inner_transition = inner_transitionRate_3_1)
    } else {
        trm <- theta_to_trans_rate_3_SM(theta,
                                        inner_transition = inner_transitionRate_3_1)
    }
    
    return(list(
        theta = theta,
        transitionRateMatrix = trm,
        transitionMatrixTimeDiscretized =
            trans_rate_to_trans_mat(trm,
                                    method = "uniformization",
                                    paranoidCheck = TRUE),
        transitionMatrixCompExp =
            trans_rate_to_trans_mat(trm,
                                    method = "competingExponentials",
                                    paranoidCheck = TRUE)
                ))
}


## transition rate matrix -> transition matrix
##    if method == uniformization, we assume no diagonal entry
##    because we compute it
##   similar code in time-discretized-CBN.R, but w/o ability to handle
##   sparse matrices

trans_rate_to_trans_mat <- function(x,
                                    method = c("competingExponentials",
                                               "uniformization"),
                                    paranoidCheck = TRUE) {
    method <- match.arg(method)
    if(method == "competingExponentials") {
        ## Conditional on there being a transition, thus diagonals are
        ## zero
        sx <- rowSums(x)
        if ( inherits(x, "dgCMatrix")) { ## sparse matrices
            tm <- x
            ii <- which(sx > 0)
            for(i in ii) {
                tm[i, ] <- tm[i, ]/sx[i]
            }
        } else {
            tm <- sweep(x, 1, sx, "/")
            tm[nrow(tm), ] <- 0 ## last row is 0
        }
        if(paranoidCheck) {
            stopifnot(isTRUE(all.equal(c(rep(1, nrow(x) - 1), 0),
                                       rowSums(tm),
                                       check.attributes = FALSE)))  
        }
    } else if (method == "uniformization") {
        ## The time-discretized version. Diagonals can be non-zero
        ## Using uniformization method, as in Schill et al., 2020,
        ## "Modelling cancer progression using Mutual Hazard Networks",
        ## Bioinformatics
        ## Fig.6 legend, who cite Grassman, 1977
        ## (Wikipedia also has an entry:
        ## https://en.wikipedia.org/wiki/Uniformization_(probability_theory )
        dd <- diag(x)
        if(!isTRUE(all(dd == 0))) stop("Diagonal of x is not 0")
        ## In p. 243 they say
        ## "diagonal elements are defined as Qxx = -\sum y;x so that
        ## columns sum to zero". They say columns, because
        ## their Q is transposed relative to what I do
        diag(x) <- -1 * rowSums(x)
        gamma <- max(abs(diag(x)))
        tm <- diag(nrow(x)) + x/gamma
        if(paranoidCheck) {
            stopifnot(isTRUE(all.equal(rep(1, nrow(x)),
                                       rowSums(tm),
                                       check.attributes = FALSE)))
        }
    }
    return(tm)
}


## integer (number of genes) -> all genotypes as 0,1 vectors
allGenotypes_3 <- function(k) {
    ## From OncoSimulR
    f1 <- function(n) {
        ## combinations: from gtools
        lapply(seq.int(n), function(x) combinations(n = n, r = x))}

    list.of.vectors <- function(y) {
        ## there's got to be a simpler way
        lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}),
                      recursive = FALSE),
               function(m) m[[1]])
    }
   
    mutated <- list.of.vectors(f1(k))
    num_mutated <- lapply(mutated, length)
    
    ## number of genes, mutated positions -> binary genotype as vector of
    ## 0, 1
    binary_genotype <- function(x, k) {
        y <- rep(0L, k)
        y[x] <- 1L
        return(y)
    }
    bin_genot <- lapply(mutated, function(x) binary_genotype(x, k = k))
    
    return(list(num_mut = c(0, unlist(num_mutated)),
                mutated = c(list(NA), mutated),
                bin_genotype = c(list(rep(0, k)), bin_genot)))

}



## theta from Learn.MHN
##      function used to compute theta -> transition rate matrix
##      Note that the diagonal is not added
##         This computes the products of the Theta, as in Fig.2 right
theta_to_trans_rate_3 <- function(theta,
                                  inner_transition = inner_transitionRate_3_1) {

    Theta <- exp(theta)
    geneNames <- colnames(theta)
    
    k <- ncol(theta)
    genots <- allGenotypes_3(k)
    numGenots <- length(genots$num_mut)

    genotNames <- unlist(
        lapply(genots$bin_genotype,
               function(x)
                   paste(geneNames[which(x == 1L)], sep = "", collapse = ", "))
        )
    genotNames[genotNames == ""] <- "WT"

    ## This single call is the one that takes most time

    TRM <- vapply(seq.int(numGenots),
           function(x)
               transitionRateC3(x, genotypes = genots,
                               Theta = Theta, maxmut = k,
                               inner_transition = inner_transition),
           double(numGenots))
    TRM <- matrix(TRM, ncol = numGenots, byrow = TRUE)
    colnames(TRM) <- rownames(TRM) <- genotNames
    return(TRM)
}



## theta from Learn.MHN
##      theta -> transition rate matrix
theta_to_trans_rate_3_SM <- function(theta,
                                     inner_transition = inner_transitionRate_3_1) {
    
    ## t1 <- Sys.time()
    Theta <- exp(theta)
    geneNames <- colnames(theta)
    k <- ncol(theta)
    genots <- allGenotypes_3(k)
    numGenots <- length(genots$num_mut)

    genotNames <- unlist(
        lapply(genots$bin_genotype,
               function(x)
                   paste(geneNames[which(x == 1L)], sep = "", collapse = ", "))
        )
    genotNames[genotNames == ""] <- "WT"

    ## Initialize sparseMatrix in first call
    trmv <- transitionRateC3_SM(1, genotypes = genots,
                                Theta = Theta, maxmut = k,
                                inner_transition = inner_transition)
    TRM <- sparseMatrix(i = rep(1, length(trmv[, "j"])),
                        j = trmv[, "j"],
                        x = trmv[, "x"],
                        dims = c(numGenots, numGenots),
                        dimnames = list(genotNames, genotNames))
    ## Can skip last one
    for(i in seq.int(2, numGenots - 1)) {
        trmv <- transitionRateC3_SM(i, genotypes = genots,
                                    Theta = Theta, maxmut = k,
                                    inner_transition = inner_transition)
        TRM[i = i, j = trmv[, "j"]] <- trmv[, "x"]
    }
    return(TRM)
}


## genotype position (row), genotype position,
##    genotype data frame
##    Theta (as exp(theta)) -> transition rate x -> y
inner_transitionRate_3_1 <- function(i, j, genotypes, Theta) {
    
    if( genotypes$num_mut[j] != (genotypes$num_mut[i] + 1) ) return(0)

    x <- genotypes$bin_genotype[[i]]
    y <- genotypes$bin_genotype[[j]]

    if( length(posy <- which(y != x)) != 1 ) return(0)

    posx <- which(x == 1L)
    if(length(posx) == 0) return(Theta[posy, posy])

    return(Theta[posy, posy] * prod(Theta[posy, posx]))

}


## row number of genotype, genotypes data frame
##            Theta (as exp(theta))
##            maximum number of mutations ->
##      transition rate x -> y [only those that are not 0 by construction]
##   output to be used for populating a sparse matrix
transitionRateC3_SM <- function(i, genotypes,  Theta, maxmut,
                            inner_transition) {

    num_genots <- length(genotypes$num_mut)
    ## Last genotype
    ## But this is silly. This code is never reached.
    if(i == num_genots) return(cbind(j = num_genots, x = 0))

    ## t11 <- Sys.time()
    ## All genotypes we are sure we cannot transition to
    nmuts <- genotypes$num_mut[i]
    ## Necessarily same of fewer mutations
    ## Necessarily not reachable
    if( (nmuts + 2) <= maxmut ) {
        mi2 <- match(nmuts + 2, genotypes$num_mut)
        upper_gg <- (mi2 - 1)
    } else {
        upper_gg <- num_genots
    }

    ## All those with same number of mutations are also
    ## necessarily not reachable.
    mi1 <- match(nmuts + 1, genotypes$num_mut)

    ff <- function(x) inner_transition(i, x, genotypes, Theta)
    qs <- vapply(mi1:upper_gg,
                 ff,
                 double(1)
                 )
    return(cbind(j = mi1:upper_gg, 
                 x = qs))
}


## row number of genotype, genotypes data frame
##            Theta (as exp(theta))
##            maximum number of mutations -> transition rate x -> y
transitionRateC3 <- function(i, genotypes,  Theta, maxmut,
                            inner_transition) {

    num_genots <- length(genotypes$num_mut)
    ## Last genotype
    if(i == num_genots) return(rep(0, num_genots)) 

    ## All genotypes we are sure we cannot transition to
    nmuts <- genotypes$num_mut[i]
    ## Necessarily same of fewer mutations
    tmp1 <- rep(0, length.out = (i - 1))
    ## Necessarily not reachable
    ## t12 <- Sys.time()
    if( (nmuts + 2) <= maxmut ) {
        mi2 <- match(nmuts + 2, genotypes$num_mut)
        tmp3 <- rep(0, length.out = num_genots - mi2 + 1)
        upper_gg <- (mi2 - 1)
    } else {
        tmp3 <- double(0)
        upper_gg <- num_genots
    }

    ## All those with same number of mutations are also
    ## necessarily not reachable.
    mi1 <- match(nmuts + 1, genotypes$num_mut)
    tmp2 <- rep(0, length.out = mi1 - i)

    ff <- function(x) inner_transition(i, x, genotypes, Theta)
    qs <- vapply(mi1:upper_gg,
                 ff,
                 double(1)
                 )
    return(c(tmp1, tmp2, qs, tmp3))
}


######################################################################
######################################################################
######################################################################
######################################################################

###   The code below has other implementations of the same functionality
###   It is now not used for computations, but it is used for testing
###   as they are different implementations.




## ## Seems hard to make it any faster
## pd3 <- profileExpr(theta_to_trans_rate_3(Theta.BC, inner_transition = inner_transitionRate_3_1))
## ## pd32 <- profileExpr(theta_to_trans_rate_3(Theta.BC, inner_transition = inner_transitionRate_3_2))
## options(width = 150)
## funSummary(pd3)
## callSummary(pd3)

## ## funSummary(pd32)
## ## callSummary(pd32)


# library(codetools)
# checkUsageEnv(env = .GlobalEnv)



## observations (rows as patients, columns genes) -> transition matrix genotypes

##          Remember my transition matrices between genotypes have origin
##          in rows, destination in columns. Transposed w.r.t. to their
##          Figure 2, left.

## Since this does not allow for sparse matrices, left only for testing.
do_MHN <- function(x,  lambda = 1/nrow(x)) {
    ## lambda 0.01 is what they use by default in their bioRxiv, p. 7 paper.
    ## In the paper it is 1/nrow(x). See paper and emails.
  
    mhnd <- Data.to.pD(x)
    message("\n      MHN: done Data.to.pD ", date(), "\n")
    theta <- Learn.MHN(mhnd, lambda = lambda)
    message("\n      MHN: done Learn.MHN ", date(), "\n")
    colnames(theta) <- rownames(theta) <- colnames(x)
    ## Reorder, so genotype names and transition rate matrix
    ## have names ordered
    oindex <- order(colnames(theta))
    theta <- theta[oindex, oindex]
    trm <- theta_to_trans_rate_3(theta,
                                 inner_transition = inner_transitionRate_3_1)
    message("\n      MHN: done theta_to_trans_rate_3 ", date(), "\n")

    return(list(
        theta = theta,
        transitionRateMatrix = trm,
        transitionMatrixTimeDiscretized =
            trans_rate_to_trans_mat(trm,
                                    method = "uniformization",
                                    paranoidCheck = TRUE),
        transitionMatrixCompExp =
            trans_rate_to_trans_mat(trm,
                                    method = "competingExponentials",
                                    paranoidCheck = TRUE)
                ))
}

## theta from Learn.MHN -> transition rate matrix
##      Note that the diagonal is not added
##         This computes the products of the Theta, as in Fig.2 right
theta_to_trans_rate_1 <- function(theta) {
    Theta <- exp(theta)
    geneNames <- colnames(theta)

    k <- ncol(theta)
    genots <- c(list(rep(0, k)), allGenotypes_former(k))
    numGenots <- length(genots)
    genotNames <- unlist(
        lapply(genots,
               function(x)
                   paste(geneNames[which(x == 1L)], sep = "", collapse = ", "))
        )
    genotNames[genotNames == ""] <- "WT"
    ## this single call is the one that takes most time
    ## t1 <- Sys.time()
    ## These calls take tiny time
    dim1 <- rep(genots, numGenots)
    dim2 <- rep(genots, rep(numGenots, numGenots))
    ## t2 <- Sys.time()
    tr2 <- function(g1, g2) transitionRateB(g1, g2, Theta)

    TRM <- mapply(tr2,
                  dim1,
                  dim2)
    ## TRM <- mapply(function(g1, g2) transitionRate(g1, g2, Theta),
    ##               dim1,
    ##               dim2)

    ## TRM <- mapply(function(g1, g2) transitionRate(g1, g2, Theta),
    ##            rep(genots, numGenots),
    ##            rep(genots, rep(numGenots, numGenots)))

    TRM <- matrix(TRM, ncol = numGenots, byrow = FALSE)
    colnames(TRM) <- rownames(TRM) <- genotNames
    return(TRM)
}

## theta from Learn.MHN
##      theta -> transition rate matrix
theta_to_trans_rate <- function(theta,
                                inner_transition = inner_transitionRate_1) {

    ## t1 <- Sys.time()
    Theta <- exp(theta)
    geneNames <- colnames(theta)
    
    k <- ncol(theta)
    genots <- allGenotypes_2(k)
    numGenots <- nrow(genots)

    genotNames <- unlist(
        lapply(genots$bin_genotype,
               function(x)
                   paste(geneNames[which(x == 1L)], sep = "", collapse = ", "))
        )
    genotNames[genotNames == ""] <- "WT"


    TRM <- vapply(seq.int(numGenots),
           function(x)
               transitionRateC(x, genotypes = genots,
                               Theta = Theta, maxmut = k,
                               inner_transition = inner_transition),
           double(numGenots))
    TRM <- matrix(TRM, ncol = numGenots, byrow = TRUE)
    colnames(TRM) <- rownames(TRM) <- genotNames
    return(TRM)
}


## integer (number of genes) -> all genotypes as 0,1 vectors
allGenotypes_2 <- function(k) {
    ## From OncoSimulR
    f1 <- function(n) {
        lapply(seq.int(n), function(x) combinations(n = n, r = x))}

    list.of.vectors <- function(y) {
        ## there's got to be a simpler way
        lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}),
                      recursive = FALSE),
               function(m) m[[1]])
    }
   
    mutated <- list.of.vectors(f1(k))
    num_mutated <- lapply(mutated, length)
    
    ## number of genes, mutated positions -> binary genotype as vector of
    ## 0, 1
    binary_genotype <- function(x, k) {
        y <- rep(0L, k)
        y[x] <- 1L
        return(y)
    }
    bin_genot <- lapply(mutated, function(x) binary_genotype(x, k = k))
    ## Map(function(nm, m) list(num_mutated = nm, genot = m), num_mutated, bin_genot)
    df <- data.frame(num_mut = unlist(num_mutated),
                     mutated = I(mutated),  ## not needed?
                     bin_genotype = I(bin_genot))
    ## Add WT
    df <- rbind(data.frame(num_mut = 0,
                           mutated = NA,
                           bin_genotype = I(list(rep(0, k)))),
                df)
    return(df)
}



## integer (number of genes) -> all genotypes as 0,1 vectors
allGenotypes_former <- function(k) {
    ## From OncoSimulR
    f1 <- function(n) {
        lapply(seq.int(n), function(x) combinations(n = n, r = x))}


    list.of.vectors <- function(y) {
        ## there's got to be a simpler way
        lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}),
                      recursive = FALSE),
               function(m) m[[1]])
    }
   
    mutated <- list.of.vectors(f1(k))
    
    ## number of genes, mutated positions -> binary genotype as vector of
    ## 0, 1
    binary_genotype <- function(x, k) {
        y <- rep(0L, k)
        y[x] <- 1L
        return(y)
    }
    lapply(mutated, function(x) binary_genotype(x, k = k))
}


## this is faster
## genotype, genotype, Theta (as exp(theta)) -> transition rate x -> y
transitionRateB <- function(x, y, Theta) {
    if( (sum(y) != (sum(x) + 1)) || (sum(y != x) != 1) ) {        
        return(0)
    } else {
        posy <- which(y != x)
        posx <- which(x == 1L)
        if(length(posx) == 0) {
            return(Theta[posy, posy])
        } else {
            return(Theta[posy, posy] * prod(Theta[posy, posx]))
        }
    }
}


## this is faster
## row number of genotype, genotypes data frame
##            Theta (as exp(theta))
##            maximum number of mutations -> transition rate x -> y
transitionRateC <- function(i, genotypes,  Theta, maxmut,
                            inner_transition) {

    num_genots <- nrow(genotypes)
    ## Last genotype
    if(i == num_genots) return(rep(0, num_genots)) 

    ## t11 <- Sys.time()
    ## All genotypes we are sure we cannot transition to
    nmuts <- genotypes[i, "num_mut"]
    ## Necessarily same of fewer mutations
    tmp1 <- rep(0, length.out = (i - 1))
    ## Necessarily not reachable
    ## t12 <- Sys.time()
    if( (nmuts + 2) <= maxmut ) {
        mi2 <- match(nmuts + 2, genotypes$num_mut)
        tmp3 <- rep(0, length.out = num_genots - mi2 + 1)
        upper_gg <- (mi2 - 1)
    } else {
        tmp3 <- double(0)
        upper_gg <- num_genots
    }
    ## All those with same number of mutations are also
    ## necessarily not reachable.
    mi1 <- match(nmuts + 1, genotypes$num_mut)
    tmp2 <- rep(0, length.out = mi1 - i)

    ff <- function(x) inner_transition(i, x, genotypes, Theta)
    qs <- vapply(mi1:upper_gg,
                 ff,
                 double(1)
                 )
    return(c(tmp1, tmp2, qs, tmp3))
}



## genotype position (row), genotype position,
##    genotype data frame
##    Theta (as exp(theta)) -> transition rate x -> y
inner_transitionRate_1 <- function(i, j, genotypes, Theta) {

    x <- genotypes[[i, "bin_genotype"]]
    y <- genotypes[[j, "bin_genotype"]]
    
    if( (genotypes[j, "num_mut"] != (genotypes[i, "num_mut"] + 1) ) ||
        (length(posy <- which(y != x)) != 1) ){
        ## (sum(y != x) != 1) ) {
        return(0)
    } else {
        ## posy <- which(y != x)
        posx <- which(x == 1L)
        if(length(posx) == 0) {
            return(Theta[posy, posy])
        } else {
            return(Theta[posy, posy] * prod(Theta[posy, posx]))
        }
    }
}


## genotype position (row), genotype position,
##    genotype data frame
##    Theta (as exp(theta)) -> transition rate x -> y
inner_transitionRate_3_2 <- function(i, j, genotypes, Theta) {

    if( genotypes$num_mut[j] != (genotypes$num_mut[i] + 1) ) return(0)

    x <- genotypes$mutated[[i]]
    y <- genotypes$mutated[[j]]

    if( length(posy <- setdiff(y, x)) != 1) return(0)

    if(all(is.na(x))) return(Theta[posy, posy])

    return(Theta[posy, posy] * prod(Theta[posy, x]))

}


## genotype position (row), genotype position,
##    genotype data frame
##    Theta (as exp(theta)) -> transition rate x -> y
inner_transitionRate_2 <- function(i, j, genotypes, Theta) {

    if(genotypes[j, "num_mut"] != (genotypes[i, "num_mut"] + 1) ) return(0)

    x <- genotypes[[i, "mutated"]]
    y <- genotypes[[j, "mutated"]]

    if( length(posy <- setdiff(y, x)) != 1) return(0)

    if(all(is.na(x))) return(Theta[posy, posy])

    return(Theta[posy, posy] * prod(Theta[posy, x]))

}



## ## genotype position (row), genotype position,
## ##    genotype data frame
## ##    Theta (as exp(theta)) -> transition rate x -> y
## inner_transitionRate_2 <- function(i, j, genotypes, Theta) {

##     x <- genotypes[[i, "mutated"]]
##     y <- genotypes[[j, "mutated"]]
    
##     if( (genotypes[j, "num_mut"] != (genotypes[i, "num_mut"] + 1) ) ||
##         (length(posy <- setdiff(y, x)) != 1) ) {
##         ## (length(setdiff(y, x)) != 1) ) {
##         return(0)
##     } else {
##         ## posy <- setdiff(y, x)
##         if(all(is.na(x))) { ## x is WT
##             return(Theta[posy, posy])
##         } else {
##             return(Theta[posy, posy] * prod(Theta[posy, x]))
##         }
##     }
## }



## ## Not used anywhere anymore
## ## genotype, genotype, Theta (as exp(theta)) -> transition rate x -> y
## transitionRate <- function(x, y, Theta) {
##     if(sum(y) != (sum(x) + 1) ) {
##         return(0)
##     } else {
##         posy <- which(y != x)
##         if(length(posy) != 1) {
##             return(0)
##         } else {
##             posx <- which(x == 1L)
##             if(length(posx) == 0) {
##                 ret <- return(Theta[posy, posy])
##             } else {
##                 ret <- (Theta[posy, posy] * prod(Theta[posy, posx]))
##             }
##             if(length(ret) > 1) {
##                 cat("\n here")
##                 stop()
##             }
##             else (return(ret))
##             ## if(length(posx) == 0) return(Theta[posy, posy])
##             ## else return(Theta[posy, posy] * cumprod(Theta[posy, posx]))
##         }
##     }
## }
