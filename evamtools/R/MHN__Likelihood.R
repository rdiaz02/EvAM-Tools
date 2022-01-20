## Files MHN__*.R: MHN__UtilityFunctions.R, MHN__RegularizedOptimization.R,
## MHN__ModelConstruction.R, MHN__Likelihood.R, MHN__InlineFunctions.R,,
## MHN__ExampleApplications.R

##   Files obtained from https://github.com/RudiSchill/MHN
##   Commit 49a8cc0 from 2018-08-16
##   (we have added the "MHN__" and made minor modifications to conform to usage
##   within an R package).
  
##   License: no license information available in the repository nor the files.
  
##   Author of code: from commit history, most likely Rudolf Schill.
  
##   Authors of paper/project: Schill, R., Solbrig, S., Wettig, T., & Spang, R.
  
##   Paper: Schill, R., Solbrig, S., Wettig, T., & Spang, R. (2020). Modelling
##   cancer progression using Mutual Hazard Networks. Bioinformatics, 36(1),
##   241–249. http://dx.doi.org/10.1093/bioinformatics/btz513

## source("InlineFunctions.R")

#Multiplies Q (as a sum of Kroncker-products) with a vector x.
Q.vec <- function(Theta, x, diag=F, transp=F){
  n <- ncol(Theta)
  y <- rep(0, 2^n)
  
  for(i in 1:n){ #Should be parallelized with MPI
    y <- y + kronvec(exp(Theta[i,]), i, x, diag, transp)
  }    
  
  return(y)
}

#Solves [-Q+I]x = b using the Jacobi method. Convergence is guaranteed
#for n+1 iterations.
Jacobi <- function(Theta, b, transp=F, x=NULL){
  n <- ncol(Theta)
  if(is.null(x)) x <- rep(1,2^n)/(2^n)
  
  dg <- -Q.Diag(Theta) + 1
  
  for(i in 1:(n+1)){
    x <- b + Q.vec(Theta, x, diag=F, transp)
    x <- x/dg
  }
  
  return(x)
}

#Generate the probability distribution from a model Theta.
Generate.pTh <- function(Theta, p0 = NULL){
  n <- ncol(Theta)
  if(is.null(p0)) p0 <- c(1, rep(0, 2^n - 1))
  
  return(Jacobi(Theta,p0))
}


#Log-likelihood Score
Score <- function(Theta, pD){
  pTh <- Generate.pTh(Theta)
  as.numeric(pD %*% log(pTh)) 
}


#Gradient of the Score wrt Theta. Implements equation (7)
Grad <- function(Theta, pD){
  n <- sqrt(length(Theta))
  Theta <- matrix(Theta,nrow=n,ncol=n)
  
  p0    <- c(1, rep(0,2^n - 1))
  
  pTh <- Jacobi(Theta, p0)
  q   <- Jacobi(Theta, pD/pTh, transp=T)

  G <- matrix(0,nrow=n,ncol=n)
  for(i in 1:n){ #should be parallelized with MPI
    r     <- q * kronvec(exp(Theta[i,]), i, pTh, 1, 0)
    G[i,] <- grad_loop_j(i,n,r)
  }

  return(G)
}
