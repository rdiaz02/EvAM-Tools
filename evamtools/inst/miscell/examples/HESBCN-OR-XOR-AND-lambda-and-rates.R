data(examples_csd)
## The n = 10 is not a sensible value here. Done for the sake of speed.
ex_hesbcn_or <- evamtools:::do_HESBCN(examples_csd$csd$OR$data, seed = 9, n = 10)

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
## Root    0 8.083 2.585 0.000 0.0000
## A       0 0.000 0.000 8.914 0.2062
## B       0 0.000 0.000 8.914 0.2062
## C       0 0.000 0.000 0.000 0.0000
## D       0 0.000 0.000 0.000 0.0000


## $parent_set
##        A        B        C        D 
## "Single" "Single"    "XOR"     "OR" 

## $epsilon
## [1] 0.08181

## $lambdas
## [1]  8.0833  2.5854 17.8277  0.4124

## $edges
##   From To      Edge Lambdas Relation
## 1 Root  A Root -> A  8.0833   Single
## 2 Root  B Root -> B  2.5854   Single
## 3    A  C    A -> C 17.8277      XOR
## 4    B  C    B -> C 17.8277      XOR
## 5    A  D    A -> D  0.4124       OR
## 6    B  D    B -> D  0.4124       OR


## The transition rate matrix ("weighted_fgraph") and the transition matrix
## between genotypes (from the former, using competing exponentials)
evamtools:::cpm2tm(ex_hesbcn_or)

## Get it annotated with column names

require(Matrix)
printSpMatrix(
    evamtools:::cpm2tm(ex_hesbcn_or)$weighted_fgraph,
    col.names = TRUE)

printSpMatrix(
    trans_mat <- evamtools:::cpm2tm(ex_hesbcn_or)$trans_mat_genots,
    col.names = TRUE)

## Output from OncoSimulR
oncosimul_out <- evamtools:::cpm2F2tm(ex_hesbcn_or$edges, max_f = NULL)





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

## Example that mixes output
d3_1 <- evamtools:::do_HESBCN(d3, seed = 26, n = 100000) ## AND, OR, XOR, Single

