# pwd0 <- getwd()
# setwd("../../code_from_what_genotype_next/")
# source("code-all-methods-minimal.R")
# setwd(pwd0)
# rm(pwd0)

library(Matrix)
test1 <- list()
test1$edges <- data.frame(
    From = c("Root", "Root", "Root", "A", "B", "C", "D", "E"),
    To = c("A", "B", "C", "D", "D", "E", "F", "F"),
    Edge = c("Root->A", "Root->B", "Root->C", "A->D", "B->D", "C->E", "D->F", "E->F"),
    Lambdas = c(1, 2, 3, 4, 4, 5, 6, 6)
)
test1$parent_set <- c("Single", "Single", "Single", "AND", "Single", "XOR")
names(test1$parent_set) <- LETTERS[1:length(test1$parent_set)]


accesible_genotypes_t1 <- c("WT", "A", "B", "C", "A, B", "A, C", "B, C", "C, E"
    , "A, B, C", "A, B, D", "B, C, E", "C, E, F", "A, C, E", "A, B, C, D", "A, B, C, E"
    , "A, B, D, F", "A, C, E, F", "B, C, E, F",  "A, B, C, D, E"
    , "A, B, C, D, F", "A, B, C, E, F")

trm1 <- Matrix(0, nrow = length(accesible_genotypes_t1)
    , ncol = length(accesible_genotypes_t1), sparse = TRUE)

rownames(trm1) <- colnames(trm1) <- accesible_genotypes_t1

trm1["WT", "A"] <- 1
trm1["WT", "B"] <- 2
trm1["WT", "C"] <- 3
trm1["A", "A, B"] <- 2
trm1["A", "A, C"] <- 3
trm1["B", "A, B"] <- 1
trm1["B", "B, C"] <- 3
trm1["C", "A, C"] <- 1
trm1["C", "B, C"] <- 2
trm1["C", "C, E"] <- 5
trm1["A, B", "A, B, C"] <- 3 
trm1["A, B", "A, B, D"] <- 4
trm1["B, C", "A, B, C"] <- 1
trm1["B, C", "B, C, E"] <- 5
trm1["A, C", "A, B, C"] <- 2
trm1["A, C", "A, C, E"] <- 5
trm1["C, E", "A, C, E"] <- 1
trm1["C, E", "B, C, E"] <- 2
trm1["C, E", "C, E, F"] <- 6
trm1["A, B, D", "A, B, D, F"] <- 6 
trm1["A, B, D", "A, B, C, D"] <- 3
trm1["A, B, C", "A, B, C, D"] <- 4
trm1["A, B, C", "A, B, C, E"] <- 5
trm1["C, E, F", "A, C, E, F"] <- 1
trm1["C, E, F", "B, C, E, F"] <- 2
trm1["B, C, E", "A, B, C, E"] <- 1
trm1["B, C, E", "B, C, E, F"] <- 6
trm1["A, C, E", "A, B, C, E"] <- 2
trm1["A, C, E", "A, C, E, F"] <- 6
trm1["A, B, D, F", "A, B, C, D, F"] <- 3
trm1["A, B, C, D", "A, B, C, D, F"] <- 6
trm1["A, B, C, D", "A, B, C, D, E"] <- 5
trm1["A, B, C, E", "A, B, C, D, E"] <- 4
trm1["A, B, C, E", "A, B, C, E, F"] <- 6
trm1["A, C, E, F", "A, B, C, E, F"] <- 2
trm1["B, C, E, F", "A, B, C, E, F"] <- 1

test2 <- list()
test2$edges <- data.frame(From = c("Root", "Root", "A", "B", "C", "D")
                    , To = c("A", "B", "C", "D", "E", "E")
                    , Lambdas = c(1, 2, 3, 4, 5, 5))
test2$parent_set <- c("Single", "Single", "Single", "Single", "XOR")
names(test2$parent_set) <- LETTERS[1:length(test2$parent_set)]
accesible_genotypes_t2 <- c("WT", "A", "B"
    , "A, B", "A, C", "B, D"
    , "A, B, C", "A, B, D", "A, C, E", "B, D, E"
    , "A, B, C, D", "A, B, C, E", "A, B, D, E"
)

trm2 <- Matrix(0, nrow = length(accesible_genotypes_t2)
    , ncol = length(accesible_genotypes_t2), sparse = TRUE)

rownames(trm2) <- colnames(trm2) <- accesible_genotypes_t2

trm2["WT", "A"] <- 1
trm2["WT", "B"] <- 2
trm2["A", "A, B"] <- 2
trm2["A", "A, C"] <- 3
trm2["B", "A, B"] <- 1
trm2["B", "B, D"] <- 4
trm2["A, B", "A, B, C"] <- 3 
trm2["A, B", "A, B, D"] <- 4
trm2["A, C", "A, B, C"] <- 2
trm2["A, C", "A, C, E"] <- 5
trm2["B, D", "A, B, D"] <- 1
trm2["B, D", "B, D, E"] <- 5
trm2["A, B, D", "A, B, D, E"] <- 5 
trm2["A, B, D", "A, B, C, D"] <- 3
trm2["A, B, C", "A, B, C, D"] <- 4
trm2["A, B, C", "A, B, C, E"] <- 5
trm2["B, D, E", "A, B, D, E"] <- 1
trm2["A, C, E", "A, B, C, E"] <- 2

test3 <- test2
test3$parent_set <- c("Single", "Single", "Single", "Single", "OR")
names(test3$parent_set) <- LETTERS[1:length(test3$parent_set)]

accesible_genotypes_t3 <- c("WT", "A", "B"
    , "A, B", "A, C", "B, D"
    , "A, B, C", "A, B, D", "A, C, E", "B, D, E"
    , "A, B, C, D", "A, B, C, E", "A, B, D, E"
    , "A, B, C, D, E"
)

trm3 <- Matrix(0, nrow = length(accesible_genotypes_t3)
    , ncol = length(accesible_genotypes_t3), sparse = TRUE)

rownames(trm3) <- colnames(trm3) <- accesible_genotypes_t3

trm3["WT", "A"] <- 1
trm3["WT", "B"] <- 2
trm3["A", "A, B"] <- 2
trm3["A", "A, C"] <- 3
trm3["B", "A, B"] <- 1
trm3["B", "B, D"] <- 4
trm3["A, B", "A, B, C"] <- 3 
trm3["A, B", "A, B, D"] <- 4
trm3["A, C", "A, B, C"] <- 2
trm3["A, C", "A, C, E"] <- 5
trm3["B, D", "A, B, D"] <- 1
trm3["B, D", "B, D, E"] <- 5
trm3["A, B, D", "A, B, D, E"] <- 5 
trm3["A, B, D", "A, B, C, D"] <- 3
trm3["A, B, C", "A, B, C, D"] <- 4
trm3["A, B, C", "A, B, C, E"] <- 5
trm3["B, D, E", "A, B, D, E"] <- 1
trm3["A, C, E", "A, B, C, E"] <- 2
trm3["A, B, C, E", "A, B, C, D, E"] <- 4
trm3["A, B, C, D", "A, B, C, D, E"] <- 5
trm3["A, B, D, E", "A, B, C, D, E"] <- 3

test4 <- test2
test4$parent_set <- c("Single", "Single", "Single", "Single", "AND")
names(test4$parent_set) <- LETTERS[1:length(test4$parent_set)]

accesible_genotypes_t4 <- c("WT", "A", "B"
    , "A, B", "A, C", "B, D"
    , "A, B, C", "A, B, D"
    , "A, B, C, D"
    , "A, B, C, D, E"
)

trm4 <- Matrix(0, nrow = length(accesible_genotypes_t4)
    , ncol = length(accesible_genotypes_t4), sparse = TRUE)

rownames(trm4) <- colnames(trm4) <- accesible_genotypes_t4

trm4["WT", "A"] <- 1
trm4["WT", "B"] <- 2
trm4["A", "A, B"] <- 2
trm4["A", "A, C"] <- 3
trm4["B", "A, B"] <- 1
trm4["B", "B, D"] <- 4
trm4["A, B", "A, B, C"] <- 3 
trm4["A, B", "A, B, D"] <- 4
trm4["A, C", "A, B, C"] <- 2
trm4["B, D", "A, B, D"] <- 1
trm4["A, B, D", "A, B, C, D"] <- 3
trm4["A, B, C", "A, B, C, D"] <- 4
trm4["A, B, C, D", "A, B, C, D, E"] <- 5


compare_hesbcn_trms <- function(computed_trm, manual_trm){
    computed_trm <- computed_trm[Matrix::rowSums(computed_trm) > 0, ]
    ## Same genotypes
    computed_genotypes <- rownames(computed_trm)
    expect_equal(length(setdiff(rownames(computed_trm), rownames(manual_trm))), 0)

    
    ## By element comparison
    filled_elements <- Matrix::summary(computed_trm)[, c("i", "j")]

    for (i in (1:length(nrow(filled_elements)))){
        from_genotype <- computed_genotypes[filled_elements[i, 1]]
        to_genotype <- computed_genotypes[filled_elements[i, 2]]

        expect_equal(
            computed_trm[from_genotype, to_genotype] ,
            computed_trm[from_genotype, to_genotype] 
        )
    }
    ## Bulk comparison
    expect_equal(sum(computed_trm), sum(manual_trm))
}

computed_trm1 <- cpm_access_genots_paths_w_simplified_relationships(test1)$weighted_fgraph
computed_trm2 <- cpm_access_genots_paths_w_simplified_relationships(test2)$weighted_fgraph
computed_trm3 <- cpm_access_genots_paths_w_simplified_relationships(test3)$weighted_fgraph
computed_trm4 <- cpm_access_genots_paths_w_simplified_relationships(test4)$weighted_fgraph

compare_hesbcn_trms(computed_trm1, trm1)
compare_hesbcn_trms(computed_trm2, trm2)
compare_hesbcn_trms(computed_trm3, trm3)
compare_hesbcn_trms(computed_trm4, trm4)