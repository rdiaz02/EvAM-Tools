## Using an example from the PMCE repository to illustrate how to interpret the
## lambdas and how we process the output.

## lambdas_matrix shows entries where the lambdas are divided
## by the number of ancestors.

##  Using repo at https://github.com/BIMIB-DISCo/PMCE
## downloaded on 2022-02-7, commit 836575a

options(width = 150)

## The following example has downloaded the repository under ~/tmp.
## Adjust for your own case.
setwd("~/tmp/PMCE")
load("./Utilities/final_results/results.RData")


##########  Breast_Invasive_Carcinoma$hesbcn
results$Breast_Invasive_Carcinoma$hesbcn
## Compare with
readLines("./Utilities/results_files/Breast_Invasive_Carcinoma_aic.txt")

## For example, CDH1, that has a XOR relationship.
results$Breast_Invasive_Carcinoma$hesbcn$lambdas_matrix[, "CDH1"]
sum(results$Breast_Invasive_Carcinoma$hesbcn$lambdas_matrix[, "CDH1"])
## That is the seventh entry in Best Lambdas:
## CDH1 is the 7th element of "parent_set"


##########   Bladder_Urothelial_Carcinoma
##  KMT2D, with an OR relationship
results$Bladder_Urothelial_Carcinoma$hesbcn

readLines("./Utilities/results_files/Bladder_Urothelial_Carcinoma_aic.txt")
results$Bladder_Urothelial_Carcinoma$hesbcn$lambdas_matrix[, "KMT2D"]
sum(results$Bladder_Urothelial_Carcinoma$hesbcn$lambdas_matrix[, "KMT2D"])
## That is the first entry in Best Lambdas;
## KMT2D is the first element of "parent_set"


##  FGFR3, with an AND relationship
readLines("./Utilities/results_files/Bladder_Urothelial_Carcinoma_aic.txt")
results$Bladder_Urothelial_Carcinoma$hesbcn$lambdas_matrix[, "FGFR3"]
sum(results$Bladder_Urothelial_Carcinoma$hesbcn$lambdas_matrix[, "FGFR3"])
## That is the 8th entry in "Best Lambdas"

