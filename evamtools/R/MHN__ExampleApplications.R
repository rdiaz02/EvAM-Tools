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



## source("UtilityFunctions.R")
## source("ModelConstruction.R")
## source("Likelihood.R")
## source("RegularizedOptimization.R")


## #Simulation-------------------------
## set.seed(1)

## #Create a true MHN with random parameters (in log-space)
## Theta.true <- Random.Theta(n=8, sparsity=0.50)
## pTh <- Generate.pTh(Theta.true)

## #Estimate the model from an empirical sample
## pD  <- Finite.Sample(pTh, 500)
## Theta.hat <- Learn.MHN(pD, lambda=1/500)
## KL.Div(pTh, Generate.pTh(Theta.hat))

## #Given the true distribution, parameters can often be recovered exactly
## Theta.rec <- Learn.MHN(pTh, lambda=0, reltol=1e-13)



## #Cancer Progression Data----------------

## Dat <- readRDS(file="data/BreastCancer.rds") 

## pD <- Data.to.pD(Dat)
## Theta.BC <- Learn.MHN(pD, lambda=0.01)

## colnames(Theta.BC) <- colnames(Dat)
## rownames(Theta.BC) <- colnames(Theta.BC)

## View(exp(Theta.BC))







