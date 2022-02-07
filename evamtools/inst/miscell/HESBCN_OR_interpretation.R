data(all_examples_csd)
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

## $edges
##   From To    Edge Lambdas Relation
## 1 Root  A Root->A  7.3644   Single
## 2 Root  B Root->B  2.3673   Single
## 3    A  C    A->C  6.5694      XOR
## 4    B  C    B->C  6.5694      XOR
## 5    A  D    A->D  0.1078       OR
## 6    B  D    B->D  0.1078       OR


trans_mat <- evamtools:::cpm_access_genots_paths_w_simplified_relationships(ex_hesbcn_or)$trans_mat_genots
oncosimul_out <- evamtools:::cpm2tm(ex_hesbcn_or$edges, max_f = NULL)


## By hand
## Transition rate: WT -> A:
ex_hesbcn_or$edges[1, "Lambdas"]/sum(ex_hesbcn_or$edges[1:2, "Lambdas"])
trans_mat[1, 1:4]

## Fitness of the three genotypes A,B ;  A,D ; B,D are all different
oncosimul_out$accessible_genotypes["A, B"]
oncosimul_out$accessible_genotypes["A, D"]
oncosimul_out$accessible_genotypes["B, D"]
