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


#Smooth approximation of the L1 penalty on Theta.
#(to be replaced with OWL-QN)
L1 <- function(Theta, eps=1e-05){
  diag(Theta) <- 0 
  sum(sqrt(Theta^2 + eps))
}

#Derivative of L1 penalty
L1_ <- function(Theta, eps=1e-05){
  diag(Theta) <- 0
  Theta / sqrt(Theta^2 + eps)
}


#Regularized Score 
Score.Reg <- function(Theta, pD, lambda){
  
  #Reshape parameters as matrix, after internal handling by BFGS as vectors.
  n <- sqrt(length(Theta))
  Theta <- matrix(Theta,nrow=n,ncol=n)  
  
  Score(Theta,pD) - lambda*L1(Theta)
} 

#Regularized Gradient
Grad.Reg  <- function(Theta, pD, lambda){
  n <- sqrt(length(Theta))
  Theta <- matrix(Theta,nrow=n,ncol=n)  
  
  Grad(Theta,pD) - lambda*L1_(Theta)
} 


#Learn an MHN from data.
Learn.MHN <- function(pD, init=NULL, lambda=0 ,maxit=5000, trace=0, reltol=1e-07, round=T){
  n <- log(length(pD), base=2)
  
  #Initialize the parameters from the independence model
  if(is.null(init)){
    init <- Learn.Indep(pD)
  } 
  
  opt <- optim(init, fn=Score.Reg, gr=Grad.Reg, pD, lambda, method="BFGS", 
               control=list(fnscale=-1,trace=trace,maxit=maxit,reltol=reltol))
  
  Theta <- matrix(opt$par,nrow=n,ncol=n)
  
  if(round){
    Theta <- round(Theta,2)
  } 
  
  return(Theta)
}


