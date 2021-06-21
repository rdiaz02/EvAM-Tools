
dB_AND <- matrix(
  c(
    rep(c(1, 0, 0, 0), 100) #A
    , rep(c(1, 1, 0, 0), 100) #AB
    , rep(c(1, 1, 1, 0), 100) #ABC
    , rep(c(1, 1, 1, 1), 100) #ABCD
    , rep(c(0, 0, 0, 0), 10) #WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_AND) <- LETTERS[1:4]

dB_OR <- matrix(
  c(
    rep(c(1, 0, 0, 0), 200) #A
    , rep(c(1, 0, 1, 0), 100) #AC
    , rep(c(1, 1, 0, 0), 100) #AB
    , rep(c(1, 1, 1, 0), 50) #ABC
    , rep(c(1, 1, 0, 1), 50) #ABD
    , rep(c(1, 0, 1, 1), 50) #ACD
    , rep(c(1, 1, 1, 1), 10) #ABCD
    , rep(c(0, 0, 0, 0), 10) #WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_OR) <- LETTERS[1:4]

dB_XOR <- matrix(
  c(
    rep(c(1, 0, 0, 0), 200) #A
    , rep(c(1, 1, 0, 0), 100) #AB
    , rep(c(1, 0, 1, 0), 100) #AC
    , rep(c(1, 1, 0, 1), 50) #ABD
    , rep(c(1, 0, 1, 1), 50) #ACD
    , rep(c(0, 0, 0, 0), 10) #WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_XOR) <- LETTERS[1:4]

# A --> B ; C --> D; B XOR D for E 
dB_c1 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), 300) #A
    , rep(c(0, 0, 1, 0, 0), 300) #C
    , rep(c(1, 1, 0, 0, 0), 200) #AB
    , rep(c(0, 0, 1, 1, 0), 200) #CD
    , rep(c(1, 1, 1, 0, 0), 100) #ABC
    , rep(c(1, 0, 1, 1, 0), 100) #ACD
    , rep(c(1, 1, 0, 0, 1), 100) #ABE
    , rep(c(0, 0, 1, 1, 1), 100) #CDE
    , rep(c(1, 1, 1, 0, 1), 100) #ABCE
    , rep(c(1, 0, 1, 1, 1), 100) #ACDE
    , rep(c(1, 1, 1, 1, 0), 50) # ABCD
    , rep(c(0, 0, 0, 0, 0), 10) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c1) <- LETTERS[1:5]

# ((A AND B) or C) to reach D
# With missing genotypes ACD + BCD
dB_c2 <- matrix(
  c(
      rep(c(1, 0, 0, 0), 300) #A
    , rep(c(0, 1, 0, 0), 300) #B
    , rep(c(0, 0, 1, 0), 300) #C
    , rep(c(1, 1, 0, 0), 200) #AB
    , rep(c(1, 0, 1, 0), 200) #AC
    , rep(c(0, 1, 1, 0), 200) #BC
    , rep(c(0, 0, 1, 1), 250) #CD
    , rep(c(1, 1, 1, 0), 100) #ABC
    , rep(c(1, 1, 0, 1), 200) #ABD
    , rep(c(1, 1, 1, 1), 200) # ABCD
    , rep(c(0, 0, 0, 0), 10) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c2) <- LETTERS[1:4]

# With ALL genotypes
dB_c2_2 <- matrix(
  c(
      rep(c(1, 0, 0, 0), 300) #A
    , rep(c(0, 1, 0, 0), 300) #B
    , rep(c(0, 0, 1, 0), 300) #C
    , rep(c(1, 1, 0, 0), 200) #AB
    , rep(c(1, 0, 1, 0), 200) #AC
    , rep(c(0, 1, 1, 0), 200) #BC
    , rep(c(0, 0, 1, 1), 250) #CD
    , rep(c(1, 1, 1, 0), 100) #ABC
    , rep(c(1, 1, 0, 1), 200) #ABD
    , rep(c(1, 0, 1, 1), 200) #ACD
    , rep(c(0, 1, 1, 1), 200) #BCD
    , rep(c(1, 1, 1, 1), 200) # ABCD
    , rep(c(0, 0, 0, 0), 10) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c2) <- LETTERS[1:4]


## WT --> A AND B --> E
## XOR 
## WT --> C AND D --> E
## ((A AND B) XOR (C AND D) to reach E
dB_c3 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), 300) #A
    , rep(c(0, 1, 0, 0, 0), 300) #B
    , rep(c(0, 0, 1, 0, 0), 300) #C
    , rep(c(0, 0, 0, 1, 0), 300) #D
    , rep(c(1, 1, 0, 0, 0), 200) #AB
    , rep(c(0, 0, 1, 1, 0), 200) #CD
    , rep(c(1, 1, 0, 0, 1), 100) #ABE
    , rep(c(0, 0, 1, 1, 1), 100) #CDE
    , rep(c(0, 0, 0, 0, 0), 10) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c3) <- LETTERS[1:5]

## WT --> A --> B
## WT --> C --> D
dB_c4 <- matrix(
  c(
      rep(c(1, 0, 0, 0), 300) #A
    , rep(c(0, 0, 1, 0), 300) #C
    , rep(c(1, 1, 0, 0), 200) #AB
    , rep(c(0, 0, 1, 1), 200) #CD
    , rep(c(1, 1, 1, 0), 100) #ABC
    , rep(c(1, 0, 1, 1), 100) #ACD
    , rep(c(1, 1, 1, 1), 50) #ABCD
    , rep(c(0, 0, 0, 0), 10) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c4) <- LETTERS[1:4]

## WT --> (A AND B AND C) --> D

dB_c5 <- matrix(
  c(
      rep(c(1, 0, 0, 0), 300) #A
    , rep(c(0, 1, 0, 0), 300) #B
    , rep(c(0, 0, 1, 0), 300) #C
    , rep(c(1, 1, 0, 0), 200) #AB
    , rep(c(1, 0, 1, 0), 200) #AC
    , rep(c(0, 1, 1, 0), 200) #BC
    , rep(c(1, 1, 1, 0), 100) #ABC
    , rep(c(1, 1, 1, 1), 50)  #ABCD
    , rep(c(0, 0, 0, 0), 10) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c5) <- LETTERS[1:4]

## Sign epistasis
dB_se <- matrix(
  c(
      rep(c(1, 0), 50) #A
    , rep(c(0, 1), 300) #B
    , rep(c(1, 1), 400) #AB
    , rep(c(0, 0), 200) # WT
  ), ncol = 2, byrow = TRUE
)
colnames(dB_se) <- LETTERS[1:2]

## Reciprocal sign epistasis
dB_rse <- matrix(
  c(
      rep(c(1, 0), 50) #A
    , rep(c(0, 1), 50) #B
    , rep(c(1, 1), 400) #AB
    , rep(c(0, 0), 200) # WT
  ), ncol = 2, byrow = TRUE
)
colnames(dB_rse) <- LETTERS[1:2]

## ((A AND B) OR (C AND D) to reach E
dB_c7 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), 300) #A
    , rep(c(0, 1, 0, 0, 0), 300) #B
    , rep(c(0, 0, 1, 0, 0), 300) #C
    , rep(c(0, 0, 0, 1, 0), 300) #D
    , rep(c(1, 1, 0, 0, 0), 200) #AB
    , rep(c(1, 0, 1, 0, 0), 200) #AC
    , rep(c(1, 0, 0, 1, 0), 200) #AD
    , rep(c(0, 1, 1, 0, 0), 200) #BC
    , rep(c(0, 1, 0, 1, 0), 200) #BD
    , rep(c(0, 0, 1, 1, 0), 200) #CD
    , rep(c(1, 1, 1, 0, 0), 100) #ABC
    , rep(c(1, 1, 0, 1, 0), 100) #ABD
    , rep(c(0, 1, 1, 1, 0), 100) #BCD
    , rep(c(1, 1, 0, 0, 1), 200) #ABE
    , rep(c(0, 0, 1, 1, 1), 200) #CDE
    , rep(c(1, 1, 1, 0, 1), 200) #ABCE
    , rep(c(1, 1, 0, 1, 1), 200) # ABDE
    , rep(c(1, 0, 1, 1, 1), 200) # ACDE
    , rep(c(0, 1, 1, 1, 1), 200) # BCDE
    , rep(c(1, 1, 1, 1, 1), 200) # ABCDE
    , rep(c(0, 0, 0, 0, 0), 10) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c7) <- LETTERS[1:5]

# Claudia Vasallo 1: representative
dB_cv1 <- matrix(
  c(
      rep(c(1, 0, 0, 0), 100) #A
    , rep(c(0, 1, 0, 0), 150) #B
    , rep(c(0, 0, 1, 0), 200) #C
    , rep(c(1, 1, 0, 0), 200) #AB
    , rep(c(1, 0, 1, 0), 250) #AC
    , rep(c(0, 1, 1, 0), 300) #BC
    , rep(c(1, 1, 1, 0), 400) #ABC
    , rep(c(0, 1, 1, 1), 450) #BCD
    , rep(c(1, 1, 1, 1), 500) #ABCD
    , rep(c(0, 0, 0, 0), 10) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_cv1) <- LETTERS[1:4]

# Claudia Vasallo 2: Local maxima
dB_cv2 <- matrix(
  c(
      rep(c(1, 0, 0, 0), 100) #A
    , rep(c(0, 1, 0, 0), 150) #B
    , rep(c(0, 0, 1, 0), 200) #C
    , rep(c(1, 1, 0, 0), 400) #AB
    , rep(c(1, 0, 1, 0), 320) #AC
    , rep(c(0, 1, 1, 0), 250) #BC
    , rep(c(1, 1, 1, 0), 250) #ABC
    , rep(c(0, 1, 1, 1), 450) #BCD
    , rep(c(1, 1, 1, 1), 500) #ABCD
    , rep(c(0, 0, 0, 0), 10) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_cv2) <- LETTERS[1:4]

# Claudia Vasallo 3: RMF
dB_cv3 <- matrix(
  c(
      rep(c(1, 0, 0, 0), 250) #A
    , rep(c(0, 1, 0, 0), 150) #B
    , rep(c(0, 0, 1, 0), 200) #C
    , rep(c(0, 1, 1, 0), 250) #BC
    , rep(c(1, 1, 1, 0), 500) #ABC
    , rep(c(0, 1, 1, 1), 250) #BCD
    , rep(c(0, 0, 0, 0), 10) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_cv3) <- LETTERS[1:4]

## Make a test with the real world data of the pmce paper

# do_HyperTraPS("Bladder_Urothelial_Carcinoma.csv", "HyperTraPS_examples/HP_pmce", runs = 1000, bi=500)
# sample_freqs(dB_OR, "HyperTraPS_examples/HP_OR/freqs.jpg")

all_examples <- list(
    AND = dB_AND 
  , OR = dB_OR 
  , XOR = dB_XOR 
  , rse = dB_rse 
  , se = dB_se 
  , c1 = dB_c1 
  , c2 = dB_c2 
  , c2_2 = dB_c2_2 
  , c3 = dB_c3 
  , c4 = dB_c4 
  , c5 = dB_c5 
  # , c6 = dB_c6 
  , c7 = dB_c7 
  , cv1 = dB_cv1
  , cv2 = dB_cv2
  , cv3 = dB_cv3
)