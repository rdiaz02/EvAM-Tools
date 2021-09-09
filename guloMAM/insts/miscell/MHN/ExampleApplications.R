
source("UtilityFunctions.R")
source("ModelConstruction.R")
source("Likelihood.R")
source("RegularizedOptimization.R")


#Simulation-------------------------
set.seed(1)

#Create a true MHN with random parameters (in log-space)
Theta.true <- Random.Theta(n=8, sparsity=0.50)
pTh <- Generate.pTh(Theta.true)

#Estimate the model from an empirical sample
pD  <- Finite.Sample(pTh, 500)
Theta.hat <- Learn.MHN(pD, lambda=1/500)
KL.Div(pTh, Generate.pTh(Theta.hat))

#Given the true distribution, parameters can often be recovered exactly
Theta.rec <- Learn.MHN(pTh, lambda=0, reltol=1e-13)



#Cancer Progression Data----------------

Dat <- readRDS(file="data/BreastCancer.rds") 

pD <- Data.to.pD(Dat)
Theta.BC <- Learn.MHN(pD, lambda=0.01)

colnames(Theta.BC) <- colnames(Dat)
rownames(Theta.BC) <- colnames(Theta.BC)

View(exp(Theta.BC))







