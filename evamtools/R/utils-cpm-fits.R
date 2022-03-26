## Copyright 2022 Ramon Diaz-Uriarte

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



## convert data frame to a matrix with 0L and 1L
df_2_mat_integer <- function(x) {
    x1 <- as.matrix(x)
    if(max(abs(x - x1)) != 0) stop("failed conversion to matrix")
    x2 <- x1
    storage.mode(x2) <- "integer"
    if(max(abs(x2 - x1)) != 0) stop("Not in 0L, 1L")
    if(max(abs(x - x2)) != 0) stop("df not in 0L, 1L") ## paranoia
    return(x2)
}

## To make it explicit
## but do not set the last row to NaNs or similar.
## Following same logic as in trans_rate_to_trans_mat
rowScaleMatrix <- function(x) {
    tm <- x
    sx <- rowSums(x)
    ii <- which(sx > 0)
    for(i in ii) {
        tm[i, ] <- tm[i, ]/sx[i]
    }
    tm
}




boot_data_index <- function(x, boot) {
    ## Used by CBN and DiP
    ## boot is an integer. 0 means no boot
    ## that is because I reuse boot for two purposes
    boot <- as.logical(boot)
    if (boot) {
        ind <- sample(nrow(x), nrow(x), replace = TRUE)
        return(x[ind, , drop = FALSE])
    } else {
        return(x)
    }
}

add_pseudosamples <- function(x, n00 = "auto3") {
    if (n00 == "auto") {
        if(nrow(x) <= 500) {
            n00 <- round(nrow(x) * 0.10)
        } else {
            n00 <- round(nrow(x) * 0.05)
        }
    } else if(n00 == "auto2") {
        ## add only if max. freq. of any gene is > 95%
        fmax <- max(colSums(x)) / nrow(x)
        if(fmax > 0.95)
            n00 <- round(nrow(x) * 0.05)
        else
            n00 <- 0
    } else if(n00 == "auto3") {
        ## add only if any gene is 100%
        ## add just 1
        fmax <- max(colSums(x))/nrow(x)
        if(fmax == 1) {
            message("\n  Added one pseudosample \n ")
            n00 <- 1
        } else { 
            n00 <- 0
        }
    }
    return(rbind(x,
                 matrix(0L, nrow = n00, ncol = ncol(x))
                 ))
}


## Remove a fraction, frac, of the WT
## from a matrix of data. Used to examine
## the effect of having very few WT
remove_WT <- function(x, frac = 0.9) {
    which_wt <- which(rowSums(x) == 0)
    if(length(which_wt) > 0) {
        rmwt <- which_wt[1:(round(frac * length(which_wt)))]
    }
    return(x[-rmwt, ])
}

## Add N WT to the data. Used to examine the effect of having many WT.
add_WT <- function(x, N = 10000) {
    ncx <- ncol(x)
    x <- rbind(x, matrix(0, nrow = N, ncol = ncx))
    return(x)
}



## ## a simple check
## any_constant_col <- function(x) {
##     nr <- nrow(x)
##     mcs <- max(colSums(x))
##     any(mcs == nr)
## }
