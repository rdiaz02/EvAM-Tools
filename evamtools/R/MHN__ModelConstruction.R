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

## require("Matrix")

#Create a random MHN with (log-transformed) parameters Theta.
#Sparsity is given as percentage.
Random.Theta <- function(n, sparsity=0){
  Theta  <- matrix(0,nrow=n,ncol=n)
  
  diag(Theta)  <- NA
  nonZeros <- sample(which(!is.na(Theta)), size=(n^2 - n)*(1 - sparsity))
  
  Theta[nonZeros] <- rnorm(length(nonZeros))
  diag(Theta) <- rnorm(n)
  
  return(round(Theta,2))
} 

#Create a single subdiagonal of Q from the ith row in Theta.
Q.Subdiag <- function(Theta, i){
  row <- Theta[i,]
  n <- length(row)
  
  #start the subdiagonal with the base rate Theta_ii 
  s <- exp(row[i])
  
  #and duplicate it for each additional factor Theta_ij.
  for(j in 1:n){
    s <- c(s, s * exp(row[j]) * (i != j))
  }
  
  return(s)
}

#Build the transition rate matrix Q from its subdiagonals.
Build.Q <- function(Theta){
  n <- nrow(Theta)
  
  Subdiags <- c()
  for(i in 1:n){
    Subdiags <- cbind(Subdiags, Q.Subdiag(Theta, i))
  }
  
  Q <- Matrix::bandSparse(2^n, k = -2^(0 : (n-1)), diagonals=Subdiags)
  diag(Q) <- -colSums(Q)
  
  return(Q)
}

#Get the diagonal of Q. 
Q.Diag <- function(Theta){
  n <- ncol(Theta)
  dg <- rep(0, 2^n)
  
  for(i in 1:n){
    dg <- dg - Q.Subdiag(Theta, i)
  }
  
  return(dg)
}

#Learn an independence model from the data distribution, which assumes that no events interact. 
#Used to initialize the parameters of the actual model before optimization.
Learn.Indep <- function(pD){
  n <- log(length(pD), base=2)
  Theta <- matrix(0, nrow=n, ncol=n)
  
  for(i in 1:n){
    pD <- matrix(pD, nrow=2^(n-1), ncol=2, byrow=T)    
    
    perc <- sum(pD[,2])
    Theta[i,i] <- log(perc/(1-perc))
  }
  
  return(round(Theta,2))
}


