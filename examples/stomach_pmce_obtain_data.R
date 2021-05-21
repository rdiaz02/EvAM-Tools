## Stomach data example
## https://github.com/BIMIB-DISCo/PMCE/blob/main/results_files/Stomach_Adenocarcinoma_aic.txt
## And it is results[[14]]
load("results.RData")
stomach <- results[[14]]$hesbcn$lambdas_matrix
pnoz <- which(stomach != 0, arr.ind = TRUE)

df1 <- data.frame(From = rownames(stomach)[pnoz[, 1]],
        To = colnames(stomach)[pnoz[, 2]], lambda = stomach[pnoz])

df1$Edge <- paste(df1$From, "->", df1$To)

relation <- results[[14]]$hesbcn$parent_set
df1$Relation <- relation[df1$To]

df1 <- df1[, c(1, 2, 4, 3, 5)]
stomach_pmce <- df1
