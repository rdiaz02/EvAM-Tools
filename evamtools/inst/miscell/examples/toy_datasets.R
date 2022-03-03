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

dB_AND <- matrix(
  c(
    rep(c(1, 0, 0, 0), round(200 * runif(1, 0.85, 1.15))) #A
    , rep(c(1, 0, 1, 0), round(100 * runif(1, 0.85, 1.15))) #AC
    , rep(c(1, 1, 0, 0), round(100 * runif(1, 0.85, 1.15))) #AB
    , rep(c(1, 1, 1, 0), round(50 * runif(1, 0.85, 1.15))) #ABC
    , rep(c(1, 1, 1, 1), round(10 * runif(1, 0.85, 1.15))) #ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.85, 1.15))) #WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_AND) <- LETTERS[1:4]

dB_OR <- matrix(
  c(
    rep(c(1, 0, 0, 0), round(400 * runif(1, 0.85, 1.15))) #A
    , rep(c(1, 0, 1, 0), round(200 * runif(1, 0.85, 1.15))) #AC
    , rep(c(1, 1, 0, 0), round(200 * runif(1, 0.85, 1.15))) #AB
    , rep(c(1, 1, 1, 0), round(100 * runif(1, 0.85, 1.15))) #ABC
    , rep(c(1, 1, 0, 1), round(150 * runif(1, 0.85, 1.15))) #ABD
    , rep(c(1, 0, 1, 1), round(150 * runif(1, 0.85, 1.15))) #ACD
    , rep(c(1, 1, 1, 1), round(50 * runif(1, 0.85, 1.15))) #ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.85, 1.15))) #WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_OR) <- LETTERS[1:4]

dB_XOR <- matrix(
  c(
    rep(c(1, 0, 0, 0), round(400 * runif(1, 0.85, 1.15))) #A
    , rep(c(1, 1, 0, 0), round(300 * runif(1, 0.85, 1.15))) #AB
    , rep(c(1, 0, 1, 0), round(300 * runif(1, 0.85, 1.15))) #AC
    , rep(c(1, 1, 1, 0), round(200 * runif(1, 0.85, 1.15))) #ABC
    , rep(c(1, 1, 0, 1), round(100 * runif(1, 0.85, 1.15))) #ABD
    , rep(c(1, 0, 1, 1), round(100 * runif(1, 0.85, 1.15))) #ACD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.85, 1.15))) #WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_XOR) <- LETTERS[1:4]

# A --> B ; C --> D; B XOR D for E 
dB_c1 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), round(300 * runif(1, 0.85, 1.15))) #A
    , rep(c(0, 0, 1, 0, 0), round(300 * runif(1, 0.85, 1.15))) #C
    , rep(c(1, 1, 0, 0, 0), round(200 * runif(1, 0.85, 1.15))) #AB
    , rep(c(0, 0, 1, 1, 0), round(200 * runif(1, 0.85, 1.15))) #CD
    , rep(c(1, 1, 1, 0, 0), round(100 * runif(1, 0.85, 1.15))) #ABC
    , rep(c(1, 0, 1, 1, 0), round(100 * runif(1, 0.85, 1.15))) #ACD
    , rep(c(1, 1, 0, 0, 1), round(100 * runif(1, 0.85, 1.15))) #ABE
    , rep(c(0, 0, 1, 1, 1), round(100 * runif(1, 0.85, 1.15))) #CDE
    , rep(c(1, 1, 1, 0, 1), round(100 * runif(1, 0.85, 1.15))) #ABCE
    , rep(c(1, 0, 1, 1, 1), round(100 * runif(1, 0.85, 1.15))) #ACDE
    , rep(c(1, 1, 1, 1, 0), round(50 * runif(1, 0.85, 1.15))) # ABCD
    , rep(c(0, 0, 0, 0, 0), round(10 * runif(1, 0.85, 1.15))) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c1) <- LETTERS[1:5]

# ((A AND B) or C) to reach D
# With missing genotypes ACD + BCD
dB_c2 <- matrix(
  c(
      rep(c(1, 0, 0, 0), round(300 * runif(1, 0.85, 1.15))) #A
    , rep(c(0, 1, 0, 0), round(300 * runif(1, 0.85, 1.15))) #B
    , rep(c(0, 0, 1, 0), round(300 * runif(1, 0.85, 1.15))) #C
    , rep(c(1, 1, 0, 0), round(200 * runif(1, 0.85, 1.15))) #AB
    , rep(c(1, 0, 1, 0), round(200 * runif(1, 0.85, 1.15))) #AC
    , rep(c(0, 1, 1, 0), round(200 * runif(1, 0.85, 1.15))) #BC
    , rep(c(0, 0, 1, 1), round(250 * runif(1, 0.85, 1.15))) #CD
    , rep(c(1, 1, 1, 0), round(100 * runif(1, 0.85, 1.15))) #ABC
    , rep(c(1, 1, 0, 1), round(200 * runif(1, 0.85, 1.15))) #ABD
    , rep(c(1, 1, 1, 1), round(200 * runif(1, 0.85, 1.15))) # ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.85, 1.15))) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c2) <- LETTERS[1:4]

# With ALL genotypes
dB_c2_2 <- matrix(
  c(
      rep(c(1, 0, 0, 0), round(300 * runif(1, 0.85, 1.15))) #A
    , rep(c(0, 1, 0, 0), round(300 * runif(1, 0.85, 1.15))) #B
    , rep(c(0, 0, 1, 0), round(300 * runif(1, 0.85, 1.15))) #C
    , rep(c(1, 1, 0, 0), round(200 * runif(1, 0.85, 1.15))) #AB
    , rep(c(1, 0, 1, 0), round(200 * runif(1, 0.85, 1.15))) #AC
    , rep(c(0, 1, 1, 0), round(200 * runif(1, 0.85, 1.15))) #BC
    , rep(c(0, 0, 1, 1), round(250 * runif(1, 0.85, 1.15))) #CD
    , rep(c(1, 1, 1, 0), round(100 * runif(1, 0.85, 1.15))) #ABC
    , rep(c(1, 1, 0, 1), round(200 * runif(1, 0.85, 1.15))) #ABD
    , rep(c(1, 0, 1, 1), round(200 * runif(1, 0.85, 1.15))) #ACD
    , rep(c(0, 1, 1, 1), round(200 * runif(1, 0.85, 1.15))) #BCD
    , rep(c(1, 1, 1, 1), round(200 * runif(1, 0.85, 1.15))) # ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.85, 1.15))) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c2_2) <- LETTERS[1:4]


## WT --> A AND B --> E
## XOR 
## WT --> C AND D --> E
## ((A AND B) XOR (C AND D) to reach E
dB_c3 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), round(300 * runif(1, 0.85, 1.15))) #A
    , rep(c(0, 1, 0, 0, 0), round(300 * runif(1, 0.85, 1.15))) #B
    , rep(c(0, 0, 1, 0, 0), round(300 * runif(1, 0.85, 1.15))) #C
    , rep(c(0, 0, 0, 1, 0), round(300 * runif(1, 0.85, 1.15))) #D
    , rep(c(1, 1, 0, 0, 0), round(200 * runif(1, 0.85, 1.15))) #AB
    , rep(c(0, 0, 1, 1, 0), round(200 * runif(1, 0.85, 1.15))) #CD
    , rep(c(1, 1, 0, 0, 1), round(100 * runif(1, 0.85, 1.15))) #ABE
    , rep(c(0, 0, 1, 1, 1), round(100 * runif(1, 0.85, 1.15))) #CDE
    , rep(c(0, 0, 0, 0, 0), round(10 * runif(1, 0.85, 1.15))) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c3) <- LETTERS[1:5]

## WT --> A --> B
## WT --> C --> D
dB_c4 <- matrix(
  c(
      rep(c(1, 0, 0, 0), round(300 * runif(1, 0.85, 1.15))) #A
    , rep(c(0, 0, 1, 0), round(300 * runif(1, 0.85, 1.15))) #C
    , rep(c(1, 1, 0, 0), round(200 * runif(1, 0.85, 1.15))) #AB
    , rep(c(0, 0, 1, 1), round(200 * runif(1, 0.85, 1.15))) #CD
    , rep(c(1, 1, 1, 0), round(100 * runif(1, 0.85, 1.15))) #ABC
    , rep(c(1, 0, 1, 1), round(100 * runif(1, 0.85, 1.15))) #ACD
    , rep(c(1, 1, 1, 1), round(50 * runif(1, 0.85, 1.15))) #ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.85, 1.15))) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c4) <- LETTERS[1:4]

## WT --> (A AND B AND C) --> D

dB_c5 <- matrix(
  c(
      rep(c(1, 0, 0, 0), round(500 * runif(1, 0.85, 1.15))) #A
    , rep(c(0, 1, 0, 0), round(500 * runif(1, 0.85, 1.15))) #B
    , rep(c(0, 0, 1, 0), round(500 * runif(1, 0.85, 1.15))) #C
    , rep(c(1, 1, 0, 0), round(400 * runif(1, 0.85, 1.15))) #AB
    , rep(c(1, 0, 1, 0), round(400 * runif(1, 0.85, 1.15))) #AC
    , rep(c(0, 1, 1, 0), round(400 * runif(1, 0.85, 1.15))) #BC
    , rep(c(1, 1, 1, 0), round(10 * runif(1, 0.85, 1.15))) #ABC
    , rep(c(1, 1, 1, 1), round(600 * runif(1, 0.85, 1.15)))  #ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.85, 1.15))) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c5) <- LETTERS[1:4]

## Sign epistasis; why should this be sign epistasis?
dB_se <- matrix(
  c(
      rep(c(1, 0), round(50 * runif(1, 0.85, 1.15))) #A
    , rep(c(0, 1), round(300 * runif(1, 0.85, 1.15))) #B
    , rep(c(1, 1), round(400 * runif(1, 0.85, 1.15))) #AB
    , rep(c(0, 0), round(200 * runif(1, 0.85, 1.15))) # WT
  ), ncol = 2, byrow = TRUE
)
colnames(dB_se) <- LETTERS[1:2]

## Reciprocal sign epistasis: why is this RSE? It ain't
dB_rse <- matrix(
  c(
      rep(c(1, 0), round(50 * runif(1, 0.85, 1.15))) #A
    , rep(c(0, 1), round(50 * runif(1, 0.85, 1.15))) #B
    , rep(c(1, 1), round(400 * runif(1, 0.85, 1.15))) #AB
    , rep(c(0, 0), round(200 * runif(1, 0.85, 1.15))) # WT
  ), ncol = 2, byrow = TRUE
)
colnames(dB_rse) <- LETTERS[1:2]

## ((A AND B) OR (C AND D) to reach E
dB_c7 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), round(300 * runif(1, 0.85, 1.15))) #A
    , rep(c(0, 1, 0, 0, 0), round(300 * runif(1, 0.85, 1.15))) #B
    , rep(c(0, 0, 1, 0, 0), round(300 * runif(1, 0.85, 1.15))) #C
    , rep(c(0, 0, 0, 1, 0), round(300 * runif(1, 0.85, 1.15))) #D
    , rep(c(1, 1, 0, 0, 0), round(200 * runif(1, 0.85, 1.15))) #AB
    , rep(c(1, 0, 1, 0, 0), round(200 * runif(1, 0.85, 1.15))) #AC
    , rep(c(1, 0, 0, 1, 0), round(200 * runif(1, 0.85, 1.15))) #AD
    , rep(c(0, 1, 1, 0, 0), round(200 * runif(1, 0.85, 1.15))) #BC
    , rep(c(0, 1, 0, 1, 0), round(200 * runif(1, 0.85, 1.15))) #BD
    , rep(c(0, 0, 1, 1, 0), round(200 * runif(1, 0.85, 1.15))) #CD
    , rep(c(1, 1, 1, 0, 0), round(100 * runif(1, 0.85, 1.15))) #ABC
    , rep(c(1, 1, 0, 1, 0), round(100 * runif(1, 0.85, 1.15))) #ABD
    , rep(c(0, 1, 1, 1, 0), round(100 * runif(1, 0.85, 1.15))) #BCD
    , rep(c(1, 1, 0, 0, 1), round(200 * runif(1, 0.85, 1.15))) #ABE
    , rep(c(0, 0, 1, 1, 1), round(200 * runif(1, 0.85, 1.15))) #CDE
    , rep(c(1, 1, 1, 0, 1), round(200 * runif(1, 0.85, 1.15))) #ABCE
    , rep(c(1, 1, 0, 1, 1), round(200 * runif(1, 0.85, 1.15))) # ABDE
    , rep(c(1, 0, 1, 1, 1), round(200 * runif(1, 0.85, 1.15))) # ACDE
    , rep(c(0, 1, 1, 1, 1), round(200 * runif(1, 0.85, 1.15))) # BCDE
    , rep(c(1, 1, 1, 1, 1), round(200 * runif(1, 0.85, 1.15))) # ABCDE
    , rep(c(0, 0, 0, 0, 0), round(10 * runif(1, 0.85, 1.15))) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c7) <- LETTERS[1:5]

# Claudia Vasallo 1: representative
dB_cv1 <- matrix(
  c(
      rep(c(1, 0, 0, 0), round(100 * runif(1, 0.85, 1.15))) #A
    , rep(c(0, 1, 0, 0), round(150 * runif(1, 0.85, 1.15))) #B
    , rep(c(0, 0, 1, 0), round(200 * runif(1, 0.85, 1.15))) #C
    , rep(c(1, 1, 0, 0), round(200 * runif(1, 0.85, 1.15))) #AB
    , rep(c(1, 0, 1, 0), round(250 * runif(1, 0.85, 1.15))) #AC
    , rep(c(0, 1, 1, 0), round(300 * runif(1, 0.85, 1.15))) #BC
    , rep(c(1, 1, 1, 0), round(400 * runif(1, 0.85, 1.15))) #ABC
    , rep(c(0, 1, 1, 1), round(450 * runif(1, 0.85, 1.15))) #BCD
    , rep(c(1, 1, 1, 1), round(500 * runif(1, 0.85, 1.15))) #ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.85, 1.15))) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_cv1) <- LETTERS[1:4]

# Claudia Vasallo 2: Local maxima
dB_cv2 <- matrix(
  c(
      rep(c(1, 0, 0, 0), round(100 * runif(1, 0.85, 1.15))) #A
    , rep(c(0, 1, 0, 0), round(150 * runif(1, 0.85, 1.15))) #B
    , rep(c(0, 0, 1, 0), round(200 * runif(1, 0.85, 1.15))) #C
    , rep(c(1, 1, 0, 0), round(400 * runif(1, 0.85, 1.15))) #AB
    , rep(c(1, 0, 1, 0), round(320 * runif(1, 0.85, 1.15))) #AC
    , rep(c(0, 1, 1, 0), round(250 * runif(1, 0.85, 1.15))) #BC
    , rep(c(1, 1, 1, 0), round(250 * runif(1, 0.85, 1.15))) #ABC
    , rep(c(0, 1, 1, 1), round(450 * runif(1, 0.85, 1.15))) #BCD
    , rep(c(1, 1, 1, 1), round(500 * runif(1, 0.85, 1.15))) #ABCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.85, 1.15))) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_cv2) <- LETTERS[1:4]

# Claudia Vasallo 3: RMF
dB_cv3 <- matrix(
  c(
      rep(c(1, 0, 0, 0), round(250 * runif(1, 0.85, 1.15))) #A
    , rep(c(0, 1, 0, 0), round(150 * runif(1, 0.85, 1.15))) #B
    , rep(c(0, 0, 1, 0), round(200 * runif(1, 0.85, 1.15))) #C
    , rep(c(0, 1, 1, 0), round(250 * runif(1, 0.85, 1.15))) #BC
    , rep(c(1, 1, 1, 0), round(500 * runif(1, 0.85, 1.15))) #ABC
    , rep(c(0, 1, 1, 1), round(250 * runif(1, 0.85, 1.15))) #BCD
    , rep(c(0, 0, 0, 0), round(10 * runif(1, 0.85, 1.15))) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_cv3) <- LETTERS[1:4]

dag_linear <- matrix(0, ncol=5, nrow=5)
rownames(dag_linear) <- colnames(dag_linear) <- c("WT", "A", "B", "C", "D")
dag_linear["WT", "A"] <- 1
dag_linear["A" , "B"] <- 1
dag_linear["B", "C"] <- 1
dag_linear["C", "D"] <- 1

dag_or <- matrix(0, ncol=5, nrow=5)
rownames(dag_or) <- colnames(dag_or) <- c("WT", "A", "B", "C", "D")
dag_or["WT", "A"] <- 1
dag_or["A" , "B"] <- 1
dag_or["A", "C"] <- 1
dag_or["C", "D"] <- 1
dag_or["B", "D"] <- 1

and_parent_set <- c("Single", "Single", "Single", "AND")
names(and_parent_set) <- c("A", "B", "C", "D")
and_lambdas <- c(1,2,3,4)
names(and_lambdas) <- c("A", "B", "C", "D")

or_parent_set <- c("Single", "Single", "Single", "OR")
names(or_parent_set) <- c("A", "B", "C", "D")

xor_parent_set <- c("Single", "Single", "Single", "XOR")
names(xor_parent_set) <- c("A", "B", "C", "D")

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

dag_se <- matrix(0, ncol=3, nrow=3)
rownames(dag_se) <- colnames(dag_se) <- c("WT", "A", "B")
dag_se["WT", "A"] <- 1
dag_se["WT", "B"] <- 1

dag_cv <- matrix(0, ncol=5, nrow=5)
rownames(dag_cv) <- colnames(dag_cv) <- c("WT", "A", "B", "C", "D")
dag_cv["WT", "A"] <- 1
dag_cv["WT", "B"] <- 1
dag_cv["WT", "C"] <- 1
dag_cv["B", "D"] <- 1
dag_cv["C", "D"] <- 1

template_dag <- matrix(0, ncol = 2, nrow = 2)
rownames(template_dag) <- colnames(template_dag) <- c("WT", "A")
template_dag["WT", "A"] <- 1

template_thetas <- matrix(0, ncol = 10, nrow = 10)
    rownames(template_thetas) <- colnames(template_thetas) <- LETTERS[1:10]

mhn_example_lambdas <- matrix(c(
  c(-0.4, -1.04, -2.53),
  c(0.18, -0.92, -0.76),
  c(-1.05, -0.53, 0.23)
), ncol = 3, nrow = 3)

mhn_example_lambdas2 <- template_thetas
mhn_example_lambdas2[1:ncol(mhn_example_lambdas)
  , 1:ncol(mhn_example_lambdas)] <- mhn_example_lambdas
rownames(mhn_example_lambdas2) <- colnames(mhn_example_lambdas2) <- c("A1", "B2", "C3", LETTERS[4:10])


examples_csd <- list(
  "csd" = list(
    User = list(data = NULL, name = "User Data"),
    Linear = list(data = dB_linear,  name = "Linear", dag = dag_linear), 
    AND = list(data = dB_AND, name = "AND", dag = dag_or, dag_parent_set = and_parent_set), 
    OR = list(data = dB_OR, name = "OR", dag = dag_or, dag_parent_set = or_parent_set), 
    XOR = list(data = dB_XOR,name = "XOR", dag = dag_or, dag_parent_set = xor_parent_set), 
    se = list(data = dB_se, name = "Sign epistasis", dag = dag_se), 
    rse = list(data = dB_rse, name = "Reciprocal sign epistasis", dag = dag_se), 
    c1 = list(data = dB_c1, name = "A --> B ; C --> D; B XOR D for E", dag = dag_c1, dag_parent_set = c1_parent_set), 
    c2 = list(ata = dB_c2, name = "Missing: ((A AND B) or C) to reach D", dag = dag_c2), 
    c2_2c2 = list(data = dB_c2_2, name = "All: ((A AND B) or C) to reach D", dag = dag_c2), 
    c4c2 = list(data = dB_c4, name = "Parallel", dag = dag_c4), 
    c3c2 = list(data = dB_c3, name = "Parallel XOR", dag = dag_c3), 
    c5c2 = list(data = dB_c5, name = "WT --> (A AND B AND C) --> D", dag = dag_c2), 
    c7c2 = list(data = dB_c7, name = "((A AND B) OR (C AND D) to reach E", dag = dag_c3), 
    cv1c2 = list(data = dB_cv1, name = "CV#1 Representative", dag = dag_cv), 
    cv2c2 = list(data = dB_cv2, name = "CV#2 Local Maxima", dag = dag_cv), 
    cv3c2 = list(data = dB_cv3, name = "CV#3 RMF", dag = dag_cv)
  ),
  "dag" = list(
    User = list(dag = template_dag, name = "User Data"),
    AND = list(data = NULL, name = "AND", dag = dag_or, dag_parent_set = and_parent_set, lambdas = and_lambdas) 
  ),
  "matrix" = list(
    User = list(thetas = template_thetas, name = "User Data"),
    test1 = list(thetas = mhn_example_lambdas2, name = "Example 1")
  )
)

# save(examples_csd, file = "../../../data/examples_csd.RData")

# cpm_output <- list()
# for (i in names(all_examples_csd_2[["csd"]])[2:11]){
#   tmp <- all_examples_csd_2[["csd"]][[i]]$data
#   cpm_output[[i]] <- evam(tmp)
# }

# save(cpm_output, file = "../../data/cpm_output.RData")
