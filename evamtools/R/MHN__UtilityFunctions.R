## Files MHN__*.R: MHN__UtilityFunctions.R, MHN__RegularizedOptimization.R,
## MHN__ModelConstruction.R, MHN__Likelihood.R, MHN__InlineFunctions.R,,
## MHN__ExampleApplications.R

##   Files obtained from https://github.com/RudiSchill/MHN
##   Commit 49a8cc0 from 2018-08-16
##   (we have added the "MHN__" and made minor modifications to conform to usage
##   within an R package).

##   License: MIT, which can be combined with software under the AGPL 3 as used
##   by the rest of this project. (Updated on 2022-04-04).
  
##   Author of code: from commit history, most likely Rudolf Schill.
  
##   Authors of paper/project: Schill, R., Solbrig, S., Wettig, T., & Spang, R.
  
##   Paper: Schill, R., Solbrig, S., Wettig, T., & Spang, R. (2020). Modelling
##   cancer progression using Mutual Hazard Networks. Bioinformatics, 36(1),
##   241–249. http://dx.doi.org/10.1093/bioinformatics/btz513


#Convert a state from a bit-vector to a natural number.
State.to.Int <- function(x){
  x <- as.logical(rev(x))
  packBits(rev(c(rep(FALSE, 32 - length(x)%%32), x)), type="integer") + 1
}

#Convert a data matrix, where each row is the bit-vector of a state,
#to a probability distribution as a vector.
Data.to.pD <- function(Data){
  Data <- as.matrix(Data)
  n <- ncol(Data)
  N <- 2^n
  
  Data <- apply(Data, 1, State.to.Int)

  pD <- tabulate(Data, nbins=N)
  pD <- pD/sum(pD)
  
  return(pD)
}

#Simulate an empirical sample from a probability distribution.
Finite.Sample <- function(pTh, k){
  N <- length(pTh)
  tabulate(sample(1:N, k, prob=pTh, replace=T), nbins=N) / k
}

#Kullback-Leibler divergence from model distribution q to true distribution p
KL.Div <- function(p,q){
  as.numeric(p%*%log(p) - p%*%log(q))
}

