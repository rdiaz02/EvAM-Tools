source("InlineFunctions.R")

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
