pwd0 <- getwd()
setwd("../code_from_what_genotype_next/")
source("code-all-methods-minimal.R")
setwd(pwd0)
rm(pwd0)

igraph_options(vertex.size = 14)
igraph_options(vertex.label.cex = 1.5)
igraph_options(edge.label.cex = 1.5)


## We need to have MCCBN installed for code below

###################################
## Simulating data under posets

## OT fails in most simple cases
## But note how the samples look in terms of genotype freqs.
## This emphasizes role of the WT, specially for CBN


true_p1 <- matrix(0, nrow = 3, ncol = 3)
true_p1[1, 2] <- 1
lambda_s <- 1
lambdas <- runif(3, 1/3*lambda_s, 3*lambda_s)
simGenotypes <- mccbn::sample_genotypes(700, true_p1,
                                        sampling_param = lambda_s,
                                        lambdas = lambdas)
db2 <- simGenotypes$obs_events
colnames(db2) <- LETTERS[1:3]
sampledGenotypes(db2)
out <- all_methods_2_trans_mat(db2)
plot_DAG_fg(out, db2)





## Simulating data under posets
true_p1 <- matrix(0, nrow = 4, ncol = 4)
true_p1[1, 2] <- 1
true_p1[3, 4] <- 1
lambda_s <- 1
lambdas <- runif(4, 1/3*lambda_s, 3*lambda_s)
simGenotypes <- mccbn::sample_genotypes(20000, true_p1,
                                        sampling_param = lambda_s,
                                        lambdas = lambdas)
db2 <- simGenotypes$obs_events
colnames(db2) <- LETTERS[1:4]
sampledGenotypes(db2)
out <- all_methods_2_trans_mat(db2)
plot_DAG_fg(out, db2)





## Simulating data under posets
true_p1 <- matrix(0, nrow = 6, ncol = 6)
true_p1[1, 2] <- 1
true_p1[2, 3] <- 1
true_p1[3, 4] <- 1
true_p1[5, 6] <- 1
lambda_s <- 1
lambdas <- runif(6, 1/6*lambda_s, 6*lambda_s)
lambdas <- runif(6, 1/2*lambda_s, lambda_s)
simGenotypes <- mccbn::sample_genotypes(20000, true_p1,
                                        sampling_param = lambda_s,
                                        lambdas = lambdas)
db2 <- simGenotypes$obs_events
colnames(db2) <- LETTERS[1:6]
sampledGenotypes(db2)
out <- all_methods_2_trans_mat(db2)
plot_DAG_fg(out, db2)

db3 <- remove_WT(db2, 1)
out3 <- all_methods_2_trans_mat(db3)
plot_DAG_fg(out3, db3)

db4 <- add_WT(db2, 100000)
out4 <- all_methods_2_trans_mat(db4)
plot_DAG_fg(out4, db4)



## Simulating data under posets
true_p1 <- matrix(0, nrow = 5, ncol = 5)
true_p1[1, 2] <- 1
true_p1[2, 3] <- 1
true_p1[4, 5] <- 1
lambda_s <- 1
lambdas <- runif(5, 1/5*lambda_s, 5*lambda_s)
simGenotypes <- mccbn::sample_genotypes(20000, true_p1,
                                        sampling_param = lambda_s,
                                        lambdas = lambdas)
db2 <- simGenotypes$obs_events
colnames(db2) <- LETTERS[1:5]
sampledGenotypes(db2)
out <- all_methods_2_trans_mat(db2)
plot_DAG_fg(out, db2)

db3 <- remove_WT(db2, 1  )
out3 <- all_methods_2_trans_mat(db3)
plot_DAG_fg(out3, db3)


