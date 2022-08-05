## This generates data/all_csd_examples.RData


dB_linear <- matrix(
  c(
    rep(c(1, 0, 0, 0), 100) #A
    , rep(c(1, 1, 0, 0), 105) #AB
    , rep(c(1, 1, 1, 0), 92) #ABC
    , rep(c(1, 1, 1, 1), 113) #ABCD
    , rep(c(0, 0, 0, 0), 7) #WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_linear) <- LETTERS[1:4]

## Add noise randomly in all the rest. Fix seed anyway
set.seed(22)

dB_AND <- matrix(
  c(
    rep(c(1, 0, 0, 0), round(200 * runif(1, 0.8, 1.2))) #A
    , rep(c(1, 0, 1, 0), round(100 * runif(1, 0.8, 1.2))) #AC
    , rep(c(1, 1, 0, 0), round(100 * runif(1, 0.8, 1.2))) #AB
    , rep(c(1, 1, 1, 0), round(50 * runif(1, 0.8, 1.2))) #ABC
    , rep(c(1, 1, 1, 1), round(10 * runif(1, 0.8, 1.2))) #ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.8, 1.2))) #WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_AND) <- LETTERS[1:4]

dB_OR <- matrix(
  c(
    rep(c(1, 0, 0, 0), round(400 * runif(1, 0.8, 1.2))) #A
    , rep(c(1, 0, 1, 0), round(200 * runif(1, 0.8, 1.2))) #AC
    , rep(c(1, 1, 0, 0), round(200 * runif(1, 0.8, 1.2))) #AB
    , rep(c(1, 1, 1, 0), round(100 * runif(1, 0.8, 1.2))) #ABC
    , rep(c(1, 1, 0, 1), round(150 * runif(1, 0.8, 1.2))) #ABD
    , rep(c(1, 0, 1, 1), round(150 * runif(1, 0.8, 1.2))) #ACD
    , rep(c(1, 1, 1, 1), round(50 * runif(1, 0.8, 1.2))) #ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.8, 1.2))) #WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_OR) <- LETTERS[1:4]

dB_XOR <- matrix(
  c(
    rep(c(1, 0, 0, 0), round(400 * runif(1, 0.8, 1.2))) #A
    , rep(c(1, 1, 0, 0), round(300 * runif(1, 0.8, 1.2))) #AB
    , rep(c(1, 0, 1, 0), round(300 * runif(1, 0.8, 1.2))) #AC
    , rep(c(1, 1, 1, 0), round(200 * runif(1, 0.8, 1.2))) #ABC
    , rep(c(1, 1, 0, 1), round(100 * runif(1, 0.8, 1.2))) #ABD
    , rep(c(1, 0, 1, 1), round(100 * runif(1, 0.8, 1.2))) #ACD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.8, 1.2))) #WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_XOR) <- LETTERS[1:4]

# A --> B ; C --> D; B XOR D for E 
dB_c1 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), round(300 * runif(1, 0.8, 1.2))) #A
    , rep(c(0, 0, 1, 0, 0), round(300 * runif(1, 0.8, 1.2))) #C
    , rep(c(1, 1, 0, 0, 0), round(200 * runif(1, 0.8, 1.2))) #AB
    , rep(c(0, 0, 1, 1, 0), round(200 * runif(1, 0.8, 1.2))) #CD
    , rep(c(1, 1, 1, 0, 0), round(100 * runif(1, 0.8, 1.2))) #ABC
    , rep(c(1, 0, 1, 1, 0), round(100 * runif(1, 0.8, 1.2))) #ACD
    , rep(c(1, 1, 0, 0, 1), round(100 * runif(1, 0.8, 1.2))) #ABE
    , rep(c(0, 0, 1, 1, 1), round(100 * runif(1, 0.8, 1.2))) #CDE
    , rep(c(1, 1, 1, 0, 1), round(100 * runif(1, 0.8, 1.2))) #ABCE
    , rep(c(1, 0, 1, 1, 1), round(100 * runif(1, 0.8, 1.2))) #ACDE
    , rep(c(1, 1, 1, 1, 0), round(50 * runif(1, 0.8, 1.2))) # ABCD
    , rep(c(0, 0, 0, 0, 0), round(10 * runif(1, 0.8, 1.2))) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c1) <- LETTERS[1:5]

# ((A AND B) or C) to reach D
# With missing genotypes ACD + BCD
dB_c2 <- matrix(
  c(
      rep(c(1, 0, 0, 0), round(300 * runif(1, 0.8, 1.2))) #A
    , rep(c(0, 1, 0, 0), round(300 * runif(1, 0.8, 1.2))) #B
    , rep(c(0, 0, 1, 0), round(300 * runif(1, 0.8, 1.2))) #C
    , rep(c(1, 1, 0, 0), round(200 * runif(1, 0.8, 1.2))) #AB
    , rep(c(1, 0, 1, 0), round(200 * runif(1, 0.8, 1.2))) #AC
    , rep(c(0, 1, 1, 0), round(200 * runif(1, 0.8, 1.2))) #BC
    , rep(c(0, 0, 1, 1), round(250 * runif(1, 0.8, 1.2))) #CD
    , rep(c(1, 1, 1, 0), round(100 * runif(1, 0.8, 1.2))) #ABC
    , rep(c(1, 1, 0, 1), round(200 * runif(1, 0.8, 1.2))) #ABD
    , rep(c(1, 1, 1, 1), round(200 * runif(1, 0.8, 1.2))) # ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.8, 1.2))) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c2) <- LETTERS[1:4]

# With ALL genotypes
dB_c2_2 <- matrix(
  c(
      rep(c(1, 0, 0, 0), round(300 * runif(1, 0.8, 1.2))) #A
    , rep(c(0, 1, 0, 0), round(300 * runif(1, 0.8, 1.2))) #B
    , rep(c(0, 0, 1, 0), round(300 * runif(1, 0.8, 1.2))) #C
    , rep(c(1, 1, 0, 0), round(200 * runif(1, 0.8, 1.2))) #AB
    , rep(c(1, 0, 1, 0), round(200 * runif(1, 0.8, 1.2))) #AC
    , rep(c(0, 1, 1, 0), round(200 * runif(1, 0.8, 1.2))) #BC
    , rep(c(0, 0, 1, 1), round(250 * runif(1, 0.8, 1.2))) #CD
    , rep(c(1, 1, 1, 0), round(100 * runif(1, 0.8, 1.2))) #ABC
    , rep(c(1, 1, 0, 1), round(200 * runif(1, 0.8, 1.2))) #ABD
    , rep(c(1, 0, 1, 1), round(200 * runif(1, 0.8, 1.2))) #ACD
    , rep(c(0, 1, 1, 1), round(200 * runif(1, 0.8, 1.2))) #BCD
    , rep(c(1, 1, 1, 1), round(200 * runif(1, 0.8, 1.2))) # ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.8, 1.2))) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c2_2) <- LETTERS[1:4]


## WT --> A AND B --> E
## XOR 
## WT --> C AND D --> E
## ((A AND B) XOR (C AND D) to reach E
dB_c3 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), round(300 * runif(1, 0.8, 1.2))) #A
    , rep(c(0, 1, 0, 0, 0), round(300 * runif(1, 0.8, 1.2))) #B
    , rep(c(0, 0, 1, 0, 0), round(300 * runif(1, 0.8, 1.2))) #C
    , rep(c(0, 0, 0, 1, 0), round(300 * runif(1, 0.8, 1.2))) #D
    , rep(c(1, 1, 0, 0, 0), round(200 * runif(1, 0.8, 1.2))) #AB
    , rep(c(0, 0, 1, 1, 0), round(200 * runif(1, 0.8, 1.2))) #CD
    , rep(c(1, 1, 0, 0, 1), round(100 * runif(1, 0.8, 1.2))) #ABE
    , rep(c(0, 0, 1, 1, 1), round(100 * runif(1, 0.8, 1.2))) #CDE
    , rep(c(0, 0, 0, 0, 0), round(10 * runif(1, 0.8, 1.2))) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c3) <- LETTERS[1:5]

## WT --> A --> B
## WT --> C --> D
dB_c4 <- matrix(
  c(
      rep(c(1, 0, 0, 0), round(300 * runif(1, 0.8, 1.2))) #A
    , rep(c(0, 0, 1, 0), round(300 * runif(1, 0.8, 1.2))) #C
    , rep(c(1, 1, 0, 0), round(200 * runif(1, 0.8, 1.2))) #AB
    , rep(c(0, 0, 1, 1), round(200 * runif(1, 0.8, 1.2))) #CD
    , rep(c(1, 1, 1, 0), round(100 * runif(1, 0.8, 1.2))) #ABC
    , rep(c(1, 0, 1, 1), round(100 * runif(1, 0.8, 1.2))) #ACD
    , rep(c(1, 1, 1, 1), round(50 * runif(1, 0.8, 1.2))) #ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.8, 1.2))) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c4) <- LETTERS[1:4]

## WT --> (A AND B AND C) --> D

dB_c5 <- matrix(
  c(
      rep(c(1, 0, 0, 0), round(500 * runif(1, 0.8, 1.2))) #A
    , rep(c(0, 1, 0, 0), round(500 * runif(1, 0.8, 1.2))) #B
    , rep(c(0, 0, 1, 0), round(500 * runif(1, 0.8, 1.2))) #C
    , rep(c(1, 1, 0, 0), round(400 * runif(1, 0.8, 1.2))) #AB
    , rep(c(1, 0, 1, 0), round(400 * runif(1, 0.8, 1.2))) #AC
    , rep(c(0, 1, 1, 0), round(400 * runif(1, 0.8, 1.2))) #BC
    , rep(c(1, 1, 1, 0), round(10 * runif(1, 0.8, 1.2))) #ABC
    , rep(c(1, 1, 1, 1), round(600 * runif(1, 0.8, 1.2)))  #ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.8, 1.2))) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c5) <- LETTERS[1:4]

## 
dB_two_ind <- matrix(
  c(
      rep(c(1, 0), round(50 * runif(1, 0.8, 1.2))) #A
    , rep(c(0, 1), round(300 * runif(1, 0.8, 1.2))) #B
    , rep(c(1, 1), round(400 * runif(1, 0.8, 1.2))) #AB
    , rep(c(0, 0), round(200 * runif(1, 0.8, 1.2))) # WT
  ), ncol = 2, byrow = TRUE
)
colnames(dB_two_ind) <- LETTERS[1:2]

## 
dB_two_ind_b <- matrix(
  c(
      rep(c(1, 0), round(50 * runif(1, 0.8, 1.2))) #A
    , rep(c(0, 1), round(50 * runif(1, 0.8, 1.2))) #B
    , rep(c(1, 1), round(400 * runif(1, 0.8, 1.2))) #AB
    , rep(c(0, 0), round(200 * runif(1, 0.8, 1.2))) # WT
  ), ncol = 2, byrow = TRUE
)
colnames(dB_two_ind_b) <- LETTERS[1:2]

## ((A AND B) OR (C AND D) to reach E
dB_c7 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), round(300 * runif(1, 0.8, 1.2))) #A
    , rep(c(0, 1, 0, 0, 0), round(300 * runif(1, 0.8, 1.2))) #B
    , rep(c(0, 0, 1, 0, 0), round(300 * runif(1, 0.8, 1.2))) #C
    , rep(c(0, 0, 0, 1, 0), round(300 * runif(1, 0.8, 1.2))) #D
    , rep(c(1, 1, 0, 0, 0), round(200 * runif(1, 0.8, 1.2))) #AB
    , rep(c(1, 0, 1, 0, 0), round(200 * runif(1, 0.8, 1.2))) #AC
    , rep(c(1, 0, 0, 1, 0), round(200 * runif(1, 0.8, 1.2))) #AD
    , rep(c(0, 1, 1, 0, 0), round(200 * runif(1, 0.8, 1.2))) #BC
    , rep(c(0, 1, 0, 1, 0), round(200 * runif(1, 0.8, 1.2))) #BD
    , rep(c(0, 0, 1, 1, 0), round(200 * runif(1, 0.8, 1.2))) #CD
    , rep(c(1, 1, 1, 0, 0), round(100 * runif(1, 0.8, 1.2))) #ABC
    , rep(c(1, 1, 0, 1, 0), round(100 * runif(1, 0.8, 1.2))) #ABD
    , rep(c(0, 1, 1, 1, 0), round(100 * runif(1, 0.8, 1.2))) #BCD
    , rep(c(1, 1, 0, 0, 1), round(200 * runif(1, 0.8, 1.2))) #ABE
    , rep(c(0, 0, 1, 1, 1), round(200 * runif(1, 0.8, 1.2))) #CDE
    , rep(c(1, 1, 1, 0, 1), round(200 * runif(1, 0.8, 1.2))) #ABCE
    , rep(c(1, 1, 0, 1, 1), round(200 * runif(1, 0.8, 1.2))) # ABDE
    , rep(c(1, 0, 1, 1, 1), round(200 * runif(1, 0.8, 1.2))) # ACDE
    , rep(c(0, 1, 1, 1, 1), round(200 * runif(1, 0.8, 1.2))) # BCDE
    , rep(c(1, 1, 1, 1, 1), round(200 * runif(1, 0.8, 1.2))) # ABCDE
    , rep(c(0, 0, 0, 0, 0), round(10 * runif(1, 0.8, 1.2))) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c7) <- LETTERS[1:5]

# Four genes, one example
dB_4g_1 <- matrix(
  c(
      rep(c(1, 0, 0, 0), round(100 * runif(1, 0.8, 1.2))) #A
    , rep(c(0, 1, 0, 0), round(150 * runif(1, 0.8, 1.2))) #B
    , rep(c(0, 0, 1, 0), round(200 * runif(1, 0.8, 1.2))) #C
    , rep(c(1, 1, 0, 0), round(200 * runif(1, 0.8, 1.2))) #AB
    , rep(c(1, 0, 1, 0), round(250 * runif(1, 0.8, 1.2))) #AC
    , rep(c(0, 1, 1, 0), round(300 * runif(1, 0.8, 1.2))) #BC
    , rep(c(1, 1, 1, 0), round(400 * runif(1, 0.8, 1.2))) #ABC
    , rep(c(0, 1, 1, 1), round(450 * runif(1, 0.8, 1.2))) #BCD
    , rep(c(1, 1, 1, 1), round(500 * runif(1, 0.8, 1.2))) #ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.8, 1.2))) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_4g_1) <- LETTERS[1:4]

# Four genes, another example
dB_4g_2 <- matrix(
  c(
      rep(c(1, 0, 0, 0), round(100 * runif(1, 0.8, 1.2))) #A
    , rep(c(0, 1, 0, 0), round(150 * runif(1, 0.8, 1.2))) #B
    , rep(c(0, 0, 1, 0), round(200 * runif(1, 0.8, 1.2))) #C
    , rep(c(1, 1, 0, 0), round(400 * runif(1, 0.8, 1.2))) #AB
    , rep(c(1, 0, 1, 0), round(320 * runif(1, 0.8, 1.2))) #AC
    , rep(c(0, 1, 1, 0), round(250 * runif(1, 0.8, 1.2))) #BC
    , rep(c(1, 1, 1, 0), round(250 * runif(1, 0.8, 1.2))) #ABC
    , rep(c(0, 1, 1, 1), round(450 * runif(1, 0.8, 1.2))) #BCD
    , rep(c(1, 1, 1, 1), round(500 * runif(1, 0.8, 1.2))) #ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.8, 1.2))) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_4g_2) <- LETTERS[1:4]

# Four genes, a third example with many unobserved genotypes.
dB_4g_3 <- matrix(
  c(
      rep(c(1, 0, 0, 0), round(250 * runif(1, 0.8, 1.2))) #A
    , rep(c(0, 1, 0, 0), round(150 * runif(1, 0.8, 1.2))) #B
    , rep(c(0, 0, 1, 0), round(200 * runif(1, 0.8, 1.2))) #C
    , rep(c(0, 1, 1, 0), round(250 * runif(1, 0.8, 1.2))) #BC
    , rep(c(1, 1, 1, 0), round(500 * runif(1, 0.8, 1.2))) #ABC
    , rep(c(0, 1, 1, 1), round(250 * runif(1, 0.8, 1.2))) #BCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.8, 1.2))) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_4g_3) <- LETTERS[1:4]

dag_linear <- matrix(0, ncol=5, nrow=5)
rownames(dag_linear) <- colnames(dag_linear) <- c("WT", "A", "B", "C", "D")
dag_linear["WT", "A"] <- 1
dag_linear["A" , "B"] <- 1
dag_linear["B", "C"] <- 1
dag_linear["C", "D"] <- 1

linear_parent_set <- c("Single", "Single", "Single", "Single")
names(linear_parent_set) <- c("A", "B", "C", "D")
linear_lambdas <- c(2 ,1 ,4 ,3)
names(linear_lambdas) <- c("A", "B", "C", "D")


dag_or <- matrix(0, ncol=5, nrow=5)
rownames(dag_or) <- colnames(dag_or) <- c("WT", "A", "B", "C", "D")
dag_or["WT", "A"] <- 1
dag_or["A" , "B"] <- 1
dag_or["A", "C"] <- 1
dag_or["C", "D"] <- 1
dag_or["B", "D"] <- 1

linear_parent_set <- c("Single", "Single", "Single", "Single")
names(linear_parent_set) <- c("A", "B", "C", "D")
and_parent_set <- c("Single", "Single", "Single", "AND")
names(and_parent_set) <- c("A", "B", "C", "D")
or_parent_set <- c("Single", "Single", "Single", "OR")
names(or_parent_set) <- c("A", "B", "C", "D")
xor_parent_set <- c("Single", "Single", "Single", "XOR")
names(xor_parent_set) <- c("A", "B", "C", "D")
and_lambdas <- c(1,2,3,4)
names(and_lambdas) <- c("A", "B", "C", "D")

or_parent_set <- c("Single", "Single", "Single", "OR")
names(or_parent_set) <- c("A", "B", "C", "D")

xor_parent_set <- c("Single", "Single", "Single", "XOR")
names(xor_parent_set) <- c("A", "B", "C", "D")




dag_and_or_xor <- matrix(0, ncol = 6, nrow = 6)
rownames(dag_and_or_xor) <- colnames(dag_and_or_xor) <- c("WT", "A", "B", "C", "D", "E")
dag_and_or_xor["WT", "A"] <- 1
dag_and_or_xor["WT", "B"] <- 1
dag_and_or_xor["A" , "C"] <- 1
dag_and_or_xor["B", "C"] <- 1
dag_and_or_xor["A", "D"] <- 1
dag_and_or_xor["B", "D"] <- 1
dag_and_or_xor["A", "E"] <- 1
dag_and_or_xor["B", "E"] <- 1
and_or_xor_parent_set <- c("Single", "Single", "AND", "OR", "XOR")
names(and_or_xor_parent_set) <- c("A", "B", "C", "D", "E")
and_or_xor_lambdas <- c(1, 2, 3, 4, 5)
names(and_or_xor_lambdas) <- c("A", "B", "C", "D", "E")



dag_c1 <- matrix(0, ncol=6, nrow=6)
rownames(dag_c1) <- colnames(dag_c1) <- c("WT", "A", "B", "C", "D", "E")
dag_c1["WT", "A"] <- 1
dag_c1["A" , "B"] <- 1
dag_c1["WT", "C"] <- 1
dag_c1["C", "D"] <- 1
dag_c1["B", "E"] <- 1
dag_c1["D", "E"] <- 1

c1_parent_set <- c("Single", "Single", "Single", "Sinlge", "XOR")
names(c1_parent_set) <- c("A", "B", "C", "D", "E")

dag_c2 <- matrix(0, ncol=5, nrow=5)
rownames(dag_c2) <- colnames(dag_c2) <- c("WT", "A", "B", "C", "D")
dag_c2["WT", "A"] <- 1
dag_c2["WT" , "B"] <- 1
dag_c2["WT", "C"] <- 1
dag_c2["A", "D"] <- 1
dag_c2["B", "D"] <- 1
dag_c2["C", "D"] <- 1

dag_c3 <- matrix(0, ncol=6, nrow=6)
rownames(dag_c3) <- colnames(dag_c3) <- c("WT", "A", "B", "C", "D", "E")
dag_c3["WT", "A"] <- 1
dag_c3["WT" , "B"] <- 1
dag_c3["WT", "C"] <- 1
dag_c3["WT", "D"] <- 1
dag_c3["A", "E"] <- 1
dag_c3["B", "E"] <- 1
dag_c3["C", "E"] <- 1
dag_c3["D", "E"] <- 1

dag_c4 <- matrix(0, ncol=5, nrow=5)
rownames(dag_c4) <- colnames(dag_c4) <- c("WT", "A", "B", "C", "D")
dag_c4["WT", "A"] <- 1
dag_c4["WT", "C"] <- 1
dag_c4["A", "B"] <- 1
dag_c4["C", "D"] <- 1

dag_2_ind <- matrix(0, ncol=3, nrow=3)
rownames(dag_2_ind) <- colnames(dag_2_ind) <- c("WT", "A", "B")
dag_2_ind["WT", "A"] <- 1
dag_2_ind["WT", "B"] <- 1

dag_4g <- matrix(0, ncol=5, nrow=5)
rownames(dag_4g) <- colnames(dag_4g) <- c("WT", "A", "B", "C", "D")
dag_4g["WT", "A"] <- 1
dag_4g["WT", "B"] <- 1
dag_4g["WT", "C"] <- 1
dag_4g["B", "D"] <- 1
dag_4g["C", "D"] <- 1

template_dag <- matrix(0, ncol = 2, nrow = 2)
rownames(template_dag) <- colnames(template_dag) <- c("WT", "A")
template_dag["WT", "A"] <- 1

template_thetas <- matrix(0, ncol = 10, nrow = 10)
    rownames(template_thetas) <- colnames(template_thetas) <- LETTERS[1:10]

## set.seed(1)
## mhn_example_lambdas <- matrix(Random.Theta(4, sparsity = 0.5))
mhn_example_lambdas <- matrix(
    c(-0.21, 0.00, 0.45, -1.24, 0, 0.41, 0, 0, 0, -1.27, -0.07, -1.22,
      0.62, 0, -0.01, -0.33), byrow = TRUE,
    ncol = 4, nrow = 4)

mhn_example_lambdas2 <- template_thetas
mhn_example_lambdas2[1:ncol(mhn_example_lambdas)
  , 1:ncol(mhn_example_lambdas)] <- mhn_example_lambdas
rownames(mhn_example_lambdas2) <- colnames(mhn_example_lambdas2) <-
    c("A", "B", "C", LETTERS[4:10])


examples_csd <- list(
  "upload" = list(
    # Empty = list(data = NULL, name = "Empty")
  ),
  "csd" = list(
    Empty = list(data = NULL, name = "Empty"),
    Linear = list(data = dB_linear,  name = "Linear", dag = dag_linear), 
    AND = list(data = dB_AND, name = "AND", dag = dag_or, dag_parent_set = and_parent_set), 
    OR = list(data = dB_OR, name = "OR", dag = dag_or, dag_parent_set = or_parent_set), 
    XOR = list(data = dB_XOR,name = "XOR", dag = dag_or, dag_parent_set = xor_parent_set), 
    two_ind = list(data = dB_two_ind, name = "Two indep. genes", dag = dag_2_ind), 
    two_ind_b = list(data = dB_two_ind_b, name = "Two indeo. genes, example b", dag = dag_2_ind), 
    c1 = list(data = dB_c1, name = "A --> B ; C --> D; B XOR D for E", dag = dag_c1, dag_parent_set = c1_parent_set), 
    c2 = list(ata = dB_c2, name = "Missing: ((A AND B) or C) to reach D", dag = dag_c2), 
    c2_2c2 = list(data = dB_c2_2, name = "All: ((A AND B) or C) to reach D", dag = dag_c2), 
    c4c2 = list(data = dB_c4, name = "Parallel", dag = dag_c4), 
    c3c2 = list(data = dB_c3, name = "Parallel XOR", dag = dag_c3), 
    c5c2 = list(data = dB_c5, name = "WT --> (A AND B AND C) --> D", dag = dag_c2), 
    c7c2 = list(data = dB_c7, name = "((A AND B) OR (C AND D) to reach E", dag = dag_c3), 
    d4gc1 = list(data = dB_4g_1, name = "Four genes, one example", dag = dag_4g), 
    d4gc2 = list(data = dB_4g_2, name = "Four genes, second example", dag = dag_4g), 
    d4gc3 = list(data = dB_4g_3, name = "Four genes, third example", dag = dag_4g)
  ),
  "dag" = list(
      DAG_Fork_3 = list(dag = template_dag, name = "DAG_Fork_3"),
      DAG_Linear = list(data = NULL,  name = "DAG_Linear", dag = dag_linear,
                        dag_parent_set = linear_parent_set,
                        lambdas = linear_lambdas
                        ),
      DAG_AND = list(data = NULL, name = "DAG_AND",
                     dag = dag_or, dag_parent_set = and_parent_set,
               lambdas = and_lambdas),
      DAG_OR = list(data = NULL, name = "DAG_OR", dag = dag_or,
                    dag_parent_set = or_parent_set), 
      DAG_XOR = list(data = NULL, name = "DAG_XOR", dag = dag_or,
                     dag_parent_set = xor_parent_set),
      DAG_A_O_X = list(data = NULL,
                       name = "DAG_A_O_X",
                       dag = dag_and_or_xor,
                      dag_parent_set = and_or_xor_parent_set,
                      lambdas = and_or_xor_lambdas)
  ),
  "matrix" = list(
     MHN_all_0 = list(thetas = template_thetas, name = "MHN_all_0"),
     MHN_Ex_1 = list(thetas = mhn_example_lambdas2, name = "MHN_Ex_1")      
  )
)

## If running from "something/EvAM-Tools/evamtools/inst/miscell/examples"
save(examples_csd, file = "../../../data/examples_csd.RData")


