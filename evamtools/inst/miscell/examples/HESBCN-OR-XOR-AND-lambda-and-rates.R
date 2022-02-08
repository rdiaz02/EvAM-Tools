data(examples_csd)
ex_hesbcn_or <- evamtools:::do_HESBCN(examples_csd$csd$OR$data, seed = 9)

## Output
## $adjacency_matrix
##      Root A B C D
## Root    0 1 1 0 0
## A       0 0 0 1 1
## B       0 0 0 1 1
## C       0 0 0 0 0
## D       0 0 0 0 0

## $lambdas_matrix
##      Root     A     B     C      D
## Root    0 7.364 2.367 0.000 0.0000
## A       0 0.000 0.000 6.569 0.1078
## B       0 0.000 0.000 6.569 0.1078
## C       0 0.000 0.000 0.000 0.0000
## D       0 0.000 0.000 0.000 0.0000

## $parent_set
##        A        B        C        D 
## "Single" "Single"    "XOR"     "OR" 

## $epsilon
## [1] 0.08162

## $lambdas
## [1]  7.3644  2.3673 13.1389  0.2156

## $edges
##   From To    Edge Lambdas Relation
## 1 Root  A Root->A  7.3644   Single
## 2 Root  B Root->B  2.3673   Single
## 3    A  C    A->C 13.1389      XOR
## 4    B  C    B->C 13.1389      XOR
## 5    A  D    A->D  0.2156       OR
## 6    B  D    B->D  0.2156       OR


## The transition rate matrix ("weighted_fgraph") and the transition matrix
## between genotypes (from the former, using competing exponentials)
evamtools:::cpm_access_genots_paths_w_simplified_relationships(ex_hesbcn_or)

## Get it annotated with column names


printSpMatrix(
    evamtools:::cpm_access_genots_paths_w_simplified_relationships(ex_hesbcn_or)$weighted_fgraph, col.names = TRUE)

## Output from OncoSimulR
oncosimul_out <- evamtools:::cpm2tm(ex_hesbcn_or$edges, max_f = NULL)

## 





## By hand
## Transition rate: WT -> A:
ex_hesbcn_or$edges[1, "Lambdas"]/sum(ex_hesbcn_or$edges[1:2, "Lambdas"])
trans_mat[1, 1:4]

## Fitness of the three genotypes A,B ;  A,D ; B,D are all different
oncosimul_out$accessible_genotypes["A, B"]
oncosimul_out$accessible_genotypes["A, D"]
oncosimul_out$accessible_genotypes["B, D"]





## A more complex example, showing AND, XOR, OR, depending on the run.

set.seed(2)
d1 <- data.frame(A = sample(c(1, 0), prob = c(0.7, 0.2), size = 200, replace = TRUE),
                 B = sample(c(1, 0), prob = c(0.85, 0.2), size = 200, replace = TRUE))
d1$C <- 0
d1$D <- 0
d1$E <- 0
d1$F <- 0
d1$C[(d1$A == 1) & (d1$B == 1)] <- 1
d1$D[(d1$A == 1) | (d1$B == 1)] <- 1
d1$D[100:200] <- 0
d1$E[xor((d1$A == 1), (d1$B == 1))] <- 1
d1$F[(d1$C == 1) & (d1$D == 1)] <- 1
d2 <- rbind(d1,
            data.frame(A = sample(c(1, 0), size = 25, prob = c(0.5, 0.2), replace = TRUE),
                       B = sample(c(1, 0), size = 25, prob = c(0.5, 0.2), replace = TRUE),
                       C = 0,
                       D = sample(c(1, 0), size = 25, replace = TRUE),
                       E = 0,
                       F = 0))
d3 <- d2
d3$C[(d3$A == 1) & (d3$B == 1)] <- 1

## Examples that mix output
d3_1 <- evamtools:::do_HESBCN(d3, seed = 26) ## AND, OR, Single
d3_2 <- evamtools:::do_HESBCN(d3, seed = 31)  ## AND, XOR, Single
